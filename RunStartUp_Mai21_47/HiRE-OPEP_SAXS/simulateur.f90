!>
!! @mainpage Hire-RNA Force Field
!!
!! Project Main Page.


program simulator
  use DEFS
  implicit none


  
#ifdef MPI
  ! Initialise the MPI interface and get the number of tasks and their id
  call MPI_Init( Error )
  call MPI_Comm_Size( MPI_COMM_WORLD,ntasks,Error )
  call MPI_Comm_Rank( MPI_COMM_WORLD,taskid,Error )
#else
  ntasks = 1
  taskid = 0
#endif 

 
  ! First read the definition file
  call definitions()

  ! Call the appropriate method
  if (SIMULATION_TYPE .eq. 'serial') then
     call serial()
  else if (SIMULATION_TYPE .eq. 'replica'  .and. REPLICA_TYPE .eq. 'T_Exchange') then
     call replica()
  else if (SIMULATION_TYPE .eq. 'tempering') then
     call site()
  else if (SIMULATION_TYPE .eq. 'replica'  .and. REPLICA_TYPE .eq. 'E_Scale') then
     call Hamiltonian_replica()
  else if (SIMULATION_TYPE .eq. 'saxs')  then
     call SAXS_guided_replica()
  else 
    write(6,*) 'Problem: Wrong simulation type'
  endif
  
#ifdef MPI
  call MPI_Finalize( Error )
#endif
  stop
end program

!******************************************************************************
!*
!> @brief Routine that launches MD if a serial mode is selected
!*
!******************************************************************************
subroutine serial()
  use DEFS
  use md_defs
  use SAXS_scoring, only: finalise_SAXS
  implicit none
  
  integer(8)             :: current_time
  type (t_conformations) :: conformation
  integer                :: imd_connected

  ! First create the structure
  conformation%temperature = target_temperature
  conformation%energyscale = 1.0         ! scales ( = 1.0)
  conformation%logfile = MASTER_LOGFILE
  conformation%path = ''
  allocate(conformation%pos(VECSIZE))
  allocate(conformation%posref(VECSIZE))
  allocate(conformation%vel(VECSIZE))

  allocate(imd_forces(VECSIZE))
  imd_forces = 0

  if (SIMULATION_METHOD .eq. 'MD') then
    call initialise_md(current_time, conformation)


#ifdef IMD_SIMULATION
    if(is_imd) then
    ! Interactive Molecular Dynamics
      call interactor_start(NATOMS, imd_forcescale, imd_wait, imd_port, imd_debug)
      call sleep(1)
      call interactor_poll(imd_connected)
      do while (imd_connected .eq. 0)
        call interactor_poll(imd_connected)
      end do
    endif
#endif

    if(restart .ne. 'restart') then  
      call thermalize_md(conformation)
    endif

    call run_md(current_time, n_production, conformation)

    ! lm759 > Finalise SAXS module
    if (compute_saxs_serial) call finalise_SAXS()

  else 
    write(*,*) 'Wrong simulation_method'
    stop
  endif

  return
end subroutine serial

!******************************************************************************
!*
!> @brief Routine that handles the replica for MD if this option is selected
!*
!******************************************************************************
subroutine replica()
  use DEFS
  use RANDOM  
  use md_defs
  use restart_module
  use SAXS_scoring, only: finalise_SAXS
  implicit none
  
  integer :: i, j, itemp, ierror, vec_length
  integer(8) :: istep, master_istep, n_steps, current_step, start_step, end_step, current_istep
  integer :: ioresult
  real(8) :: en1, en2, T1, T2, ran_number
  real(8) :: ran3

  logical file_exists
  character(len=2000) :: line
  character(len=100) :: commande, commande1, commande2, fname, fname2, checkname
  character(len=30) :: word,chaine, lpath
  type (t_conformations),dimension(:), allocatable :: conformations
  real(8), dimension(:), allocatable :: vecteur
  integer :: REMDex_idx

  ! Prepare the vector to be used in the data exchange between the nodes. Contains the information
  ! of the structure t_conformations
  vec_length = 3*vecsize+5
  allocate(vecteur(vec_length))

  ! Generate the conformation information
  allocate(conformations(N_REPLICA))
  do i =1, N_REPLICA
    conformations(i)%id = i
    conformations(i)%temperature = T_replica(i)
    conformations(i)%energyscale = 1.0         ! scales ( = 1.0)
    allocate(conformations(i)%pos(VECSIZE))
    allocate(conformations(i)%posref(VECSIZE))
    allocate(conformations(i)%vel(VECSIZE))
  enddo
  
  ! If we have a restart, then we read the main status of the run. 
  ! All the positions and velocities are read in initial_md
  if(restart .eq. 'restart') call read_master_restart(master_istep)
    
  ! Loop for each node, if there is less nodes than replicas, then some will
  ! do two  or more temperatures. The number of replicas divided by the number
  ! of nodes must be an integer. 
  do i=1, N_REPLICA/ntasks
    T_id = (i-1) *ntasks + taskid +1

    ! Now, check that the subdirectories are there
    call convert_to_chain(T_id,chaine,"0")
    lpath = 'p' // chaine(29:30)
    
    ! Regrettably, inquire cannot tell whether or not a directory exists,
    ! we must therefore base our information on a file test.existence
    checkname = trim(lpath) //"/test.existence"
    open(unit=83,file=checkname,status='unknown',iostat=ioresult)
    close(83)
    
    if (ioresult .eq. 0 ) then
      file_exists = .true.
    else
      file_exists = .false.
    endif

    ! If they do not exist, then create them
    if (.not.file_exists) then 
      commande = 'mkdir ' // trim(lpath) // char(0)
#ifdef MP
      call lsystem(commande)  ! Bug on MP requires a different call
#else
      call system(commande)
#endif
    endif
    ! Store the name of the logfile and of the paths
    conformations(T_id)%logfile = trim(lpath)// '/' // trim(MASTER_LOGFILE)  // '.' // chaine(29:30)
    conformations(T_id)%path = trim(lpath) // '/'
  enddo
    
  ! We first initalise and thermalize the various simulations
  do i=1, N_REPLICA/ntasks
    T_id = (i-1)*ntasks + taskid + 1  
    if (SIMULATION_METHOD .eq. 'MD') then
      call initialise_md(current_istep,conformations(T_id))
      if(restart .ne. 'restart') call thermalize_md(conformations(T_id))
    endif
  end do

  ! Then we prepare the header for the switches
  if (taskid .eq. 0) call switch_header()

#ifdef MPI
    ! We now wait here until MD simulations are finished
    call MPI_Barrier(MPI_COMM_WORLD,ierror)
    ! And then, we collect back the information to the central node
    call  collect_info_to_main_node(n_replica,conformations)
#endif

    ! We then start the simulation per se
    if (SIMULATION_METHOD .eq. 'MD') then
      n_steps = n_production / n_step_exchange
    else
      STOP 'Bad simulation method. Aborting.'
    endif

    start_step = 1
    if (restart .eq. 'restart') start_step = master_istep+1
    
    do istep = start_step, n_steps
      current_step = (istep-1) * n_step_exchange
      end_step = istep * n_step_exchange
      do j=1, N_REPLICA/ntasks
        T_id = (j-1) * ntasks + taskid + 1
        
        ! If we restart, then if the simulation was waiting for exchange, we do not 
        ! launch a new simulation. Otherwise, we do that.
        if (restart .eq. 'restart' .and. current_istep .gt. master_istep) exit  

        if (SIMULATION_METHOD .eq. 'MD') then
          call run_md(current_step, end_step, conformations(T_id))
        endif
        
        ! We now save the current simulation
        call save_restart(istep,conformations(T_id))
      end do
      
      ! We can reset the restart as we are now continuing our simulation normally
      restart = 'new'
           
#ifdef MPI
      ! We now wait here until MD simulations are finished
      call MPI_Barrier(MPI_COMM_WORLD,ierror)
      ! And then, we collect back the information to the central node
      call  collect_info_to_main_node(n_replica,conformations)
      call MPI_Barrier(MPI_COMM_WORLD,ierror)
#endif

      if (taskid .eq. 0) then
        do REMDex_idx = 1, 100
        ! if istep is even, j=2, if odd, j=3
!        do j=modulo(istep,2)+2, N_REPLICA,2
         ! ran3 is in [0;1), j is in [2;N] ( or [2;N+1) )
          j = int(ran3()*(N_REPLICA-1)+2)
          en1 = conformations(j-1)%energy
          en2 = conformations(j)%energy

          T1  = conformations(j-1)%temperature
          T2  = conformations(j)%temperature
      
          ran_number = ran3()
      
          ! p(T2->T1) = min(1, exp{-[(beta2-beta1)(en1-en2)]}
          if (ran_number .le. exp((1.0d0/T1 - 1.0d0/T2)* (en1-en2))) then
            itemp = conformations(j-1)%id
            conformations(j-1)%id = conformations(j)%id
            conformations(j)%id = itemp
            call switch(vecsize,conformations(j-1)%pos,   conformations(j)%pos)
            call switch(vecsize,conformations(j-1)%posref,conformations(j)%posref)
            call switch(vecsize,conformations(j-1)%vel,   conformations(j)%vel)
          
            ! Renormalise the velocities by the ratio of temperatures after the exchange
            if (SIMULATION_METHOD .eq. 'MD') then    
              conformations(j-1)%vel(:) = sqrt(T1/T2)*conformations(j-1)%vel(:)
              conformations(j  )%vel(:) = sqrt(T2/T1)*conformations(j)%vel(:)
            endif
          endif
!        end do
        end do


        ! Write the output of the exchange in replica.dat 
        open(unit=FREP,file=REPLICAFILE,status='unknown',action='write',position='append',iostat=ierror)
        open(unit=FLOG,file=MASTER_LOGFILE,status='unknown',action='write',position='append',iostat=ierror)
        write(word,'(i8,i10)') istep,conformations(1)%counter
        line = '  ' // word
        
        do j=1, N_REPLICA
          write(word,'(i8)') conformations(j)%id
          line = trim(line) //word
        end do
        write(FLOG,"(A)") trim(line)
        write(FREP,"(A)") trim(line)
        close(FREP)
        close(FLOG)
      endif

#ifdef MPI
      ! We now broadcast the information to the other nodes
      call MPI_Barrier(MPI_COMM_WORLD,ierror)
      do i=1, N_REPLICA
        call conformation_to_vector(vec_length,vecteur,conformations(i))
        vecteur(vec_length) = i
        call MPI_Bcast(vecteur,vec_length,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
        T_id = nint(vecteur(vec_length))
        call vector_to_conformation(vec_length,vecteur,conformations(T_id))
      end do
#endif

      ! The following makes sure that we always have the backup restart files all
      ! at the same time - 
      do j=1, N_REPLICA/ntasks
        T_id = (j-1) * ntasks + taskid + 1
        fname = trim(conformations(T_id)%path) // RESTARTFILE
        commande = 'cp ' // trim(fname)
        fname = trim(conformations(T_id)%path) // 'restart.backup'
        commande1 = trim(commande) // ' ' // trim(fname) // char(0)
        commande = 'cp ' // trim(fname)
        fname2 = trim(fname) // '.1'
        commande2 = trim(commande) // ' ' // trim(fname2) // char(0)

#ifdef MP
        call lsystem(trim(commande2))  ! Bug on MP requires a different call
#else
        call system(commande2)
#endif

#ifdef MP
        call lsystem(trim(commande1))  ! Bug on MP requires a different call
#else
        call system(commande1)

#endif

      end do

      ! And we save the master restart information and update the backups
      do j=1, N_REPLICA/ntasks
        T_id = (j-1) * ntasks + taskid + 1
        ! We now save the current simulation
        call save_restart(istep,conformations(T_id))
      end do
        
      call save_master_restart(istep)
    end do

  ! lm759 > Finalise SAXS module
  if (compute_saxs_serial) call finalise_SAXS()

  deallocate(vecteur)
  
  return
end subroutine replica


!******************************************************************************
!*
!> @brief Routine that handles the replica for MD if this option is selected
!*
!******************************************************************************
subroutine SAXS_guided_replica()
  use DEFS
  use RANDOM
  use md_defs
  use restart_module
  use SAXS_scoring
  implicit none

  integer :: i, j, itemp, ierror, vec_length, current_istep
  integer(8) :: istep, master_istep, n_steps, current_step, start_step, end_step
  integer :: n_rate = 0
  integer :: ioresult
  real(8) :: en1, en2, T1, T2, ran_number
  real(8) :: delta_SAXS
  real(8) :: ran3
#ifdef DEBUG_EXCHANGE_RATE
  integer :: attempts = 0
  integer :: exchange = 0
#endif

  logical file_exists
  character(len=2000) :: line
  character(len=100) :: commande, commande1, commande2, fname, fname2, checkname
  character(len=30) :: word,chaine, lpath
  type (t_conformations),dimension(:), allocatable :: conformations
  real(8), dimension(:), allocatable :: vecteur
  integer :: REMDex_idx

  ! Prepare the vector to be used in the data exchange between the nodes. Contains the information
  ! of the structure t_conformations
  vec_length = 3*vecsize+6
  allocate(vecteur(vec_length))

  ! Generate the conformation information
  allocate(conformations(N_REPLICA))
  do i =1, N_REPLICA
    conformations(i)%id = i
    conformations(i)%temperature = T_replica(i)
    conformations(i)%energyscale = 1.0         ! scales ( = 1.0)
    conformations(i)%score = 0.0
    allocate(conformations(i)%pos(VECSIZE))
    allocate(conformations(i)%posref(VECSIZE))
    allocate(conformations(i)%vel(VECSIZE))
  enddo

  ! If we have a restart, then we read the main status of the run.
  ! All the positions and velocities are read in initial_md
  if(restart .eq. 'restart') call read_master_restart(master_istep)

  ! Loop for each node, if there is less nodes than replicas, then some will
  ! do two  or more temperatures. The number of replicas divided by the number
  ! of nodes must be an integer.
  do i=1, N_REPLICA/ntasks
    T_id = (i-1) *ntasks + taskid +1

    ! Now, check that the subdirectories are there
    call convert_to_chain(T_id,chaine,"0")
    lpath = 'p' // chaine(29:30)

    ! Regrettably, inquire cannot tell whether or not a directory exists,
    ! we must therefore base our information on a file test.existence
    checkname = trim(lpath) //"/test.existence"
    open(unit=83,file=checkname,status='unknown',iostat=ioresult)
    close(83)

    if (ioresult .eq. 0 ) then
      file_exists = .true.
    else
      file_exists = .false.
    endif

    ! If they do not exist, then create them
    if (.not.file_exists) then
      commande = 'mkdir ' // trim(lpath) // char(0)
#ifdef MP
      call lsystem(commande)  ! Bug on MP requires a different call
#else
      call system(commande)
#endif
    endif
    ! Store the name of the logfile and of the paths
    conformations(T_id)%logfile = trim(lpath)// '/' // trim(MASTER_LOGFILE)  // '.' // chaine(29:30)
    conformations(T_id)%path = trim(lpath) // '/'
  enddo

  ! We first initalise and thermalize the various simulations
  do i=1, N_REPLICA/ntasks
    T_id = (i-1)*ntasks + taskid + 1
    if (SIMULATION_METHOD .eq. 'MD') then
      call initialise_md(current_istep,conformations(T_id))
      if(restart .ne. 'restart') call thermalize_md(conformations(T_id))
    endif
  end do

  if (taskid .eq. 0) then
    ! Then we prepare the header for the switches
    call switch_header()

    ! Finally, ensure consistency of system with regard to the target curve
    call check_SAXS_consistency(conformations(1)%pos)
  end if

#ifdef MPI
    ! We now wait here until MD simulations are finished
    call MPI_Barrier(MPI_COMM_WORLD,ierror)
    ! And then, we collect back the information to the central node
    call  collect_info_to_main_node(n_replica,conformations)
#endif

    ! We then start the simulation per se
    if (SIMULATION_METHOD .eq. 'MD') then
      n_steps = n_production / n_step_exchange
    else
      STOP 'Bad simulation method. Aborting.'
    endif

    start_step = 1
    if (restart .eq. 'restart') start_step = master_istep+1

    do istep = start_step, n_steps
      current_step = (istep-1) * n_step_exchange
      end_step = istep * n_step_exchange
      do j=1, N_REPLICA/ntasks
        T_id = (j-1) * ntasks + taskid + 1

        ! If we restart, then if the simulation was waiting for exchange, we do not
        ! launch a new simulation. Otherwise, we do that.
        if (restart .eq. 'restart' .and. current_istep .gt. master_istep) exit

        if (SIMULATION_METHOD .eq. 'MD') then
          call run_md(current_step, end_step, conformations(T_id))
        endif

        ! We now save the current simulation
        call save_restart(istep,conformations(T_id))

        ! And finally, we update the SAXS score and log it for each replica
        update_saxs_score: if ( modulo(n_rate,n_rate_saxs) .eq. 0 ) then
          fname = trim(conformations(T_id)%path) // 'saxs_curve_current.dat'
          open(unit=50, file=fname)
          conformations(T_id)%score = SAXS_energy(conformations(T_id)%pos, write_to_unit=50)
          close(50)
          fname = trim(conformations(T_id)%path) // 'saxs.log'
          open (unit=50, file=fname, status='unknown',action='write',position='append',iostat=ierror)
          write(50,*) conformations(T_id)%score
          close(50)
        end if update_saxs_score
      end do
      n_rate = n_rate+1

      ! We can reset the restart as we are now continuing our simulation normally
      restart = 'new'

#ifdef MPI
      ! We now wait here until MD simulations are finished
      call MPI_Barrier(MPI_COMM_WORLD,ierror)
      ! And then, we collect back the information to the central node
      call  collect_info_to_main_node(n_replica,conformations)
      call MPI_Barrier(MPI_COMM_WORLD,ierror)
#endif

      if (taskid .eq. 0) then

        do REMDex_idx = 1, 100
        ! if istep is even, j=2, if odd, j=3
!        do j=modulo(istep,2)+2, N_REPLICA,2
         ! ran3 is in [0;1), j is in [2;N] ( or [2;N+1) )
          j = int(ran3()*(N_REPLICA-1)+2)
          en1 = conformations(j-1)%energy
          en2 = conformations(j)%energy

          T1  = conformations(j-1)%temperature
          T2  = conformations(j)%temperature

          delta_SAXS = log10(conformations(j-1)%score / conformations(j)%score)

          ran_number = ran3()

          ! p(T2->T1) = min(1, exp{-[(beta2-beta1)(en1-en2)]}
          if (ran_number .le. exp((1.0d0/T1 - 1.0d0/T2)* (en1-en2 + saxs_alpha * delta_SAXS))) then
#ifdef DEBUG_EXCHANGE_RATE
            attempts = attempts + 1
            exchange = exchange + 1
            print *, exchange, '/', attempts
#endif
            itemp = conformations(j-1)%id
            conformations(j-1)%id = conformations(j)%id
            conformations(j)%id = itemp
            call switch(vecsize,conformations(j-1)%pos,   conformations(j)%pos)
            call switch(vecsize,conformations(j-1)%posref,conformations(j)%posref)
            call switch(vecsize,conformations(j-1)%vel,   conformations(j)%vel)

            ! Renormalise the velocities by the ratio of temperatures after the exchange
            if (SIMULATION_METHOD .eq. 'MD') then
              conformations(j-1)%vel(:) = sqrt(T1/T2)*conformations(j-1)%vel(:)
              conformations(j  )%vel(:) = sqrt(T2/T1)*conformations(j)%vel(:)
            endif
#ifdef DEBUG_EXCHANGE_RATE
          else
            attempts = attempts + 1
            print *, exchange, '/', attempts
#endif
          endif
!        end do
        end do


        ! Write the output of the exchange in replica.dat
        open(unit=FREP,file=REPLICAFILE,status='unknown',action='write',position='append',iostat=ierror)
        open(unit=FSAXS,file=SAXSFILE,status='unknown',action='write',position='append',iostat=ierror)
        open(unit=FLOG,file=MASTER_LOGFILE,status='unknown',action='write',position='append',iostat=ierror)
        write(word,'(i8,i10)') istep,conformations(1)%counter
        write(FSAXS,'(i8,i10)') istep,conformations(1)%counter
        line = '  ' // word

        do j=1, N_REPLICA
          write(word,'(i8)') conformations(j)%id
          line = trim(line) //word
          write(FSAXS,'(i4,ES15.7)') j, conformations(j)%score
        end do
        write(FLOG,"(A)") trim(line)
        write(FREP,"(A)") trim(line)
        close(FREP)
        close(FLOG)
        close(FSAXS)
      endif

#ifdef MPI
      ! We now broadcast the information to the other nodes
      call MPI_Barrier(MPI_COMM_WORLD,ierror)
      do i=1, N_REPLICA
        call conformation_to_vector(vec_length,vecteur,conformations(i))
        vecteur(vec_length) = i
        call MPI_Bcast(vecteur,vec_length,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
        T_id = nint(vecteur(vec_length))
        call vector_to_conformation(vec_length,vecteur,conformations(T_id))
      end do
#endif

      ! The following makes sure that we always have the backup restart files all
      ! at the same time -
      do j=1, N_REPLICA/ntasks
        T_id = (j-1) * ntasks + taskid + 1
        fname = trim(conformations(T_id)%path) // RESTARTFILE
        commande = 'cp ' // trim(fname)
        fname = trim(conformations(T_id)%path) // 'restart.backup'
        commande1 = trim(commande) // ' ' // trim(fname) // char(0)
        commande = 'cp ' // trim(fname)
        fname2 = trim(fname) // '.1'
        commande2 = trim(commande) // ' ' // trim(fname2) // char(0)

#ifdef MP
        call lsystem(trim(commande2))  ! Bug on MP requires a different call
#else
        call system(commande2)
#endif

#ifdef MP
        call lsystem(trim(commande1))  ! Bug on MP requires a different call
#else
        call system(commande1)

#endif

      end do

      ! And we save the master restart information and update the backups
      do j=1, N_REPLICA/ntasks
        T_id = (j-1) * ntasks + taskid + 1
        ! We now save the current simulation
        call save_restart(istep,conformations(T_id))
      end do

      call save_master_restart(istep)
    end do

  ! lm759 > Finalise SAXS module
  call finalise_SAXS()

  deallocate(vecteur)

  return
end subroutine SAXS_guided_replica

!########################################
!******************************************************************************
!*
!> @brief Routine that handles the Hamiltonian replica exchange if this option is selected
!*           (Rozita 30 Oct 08)
!*
!******************************************************************************
subroutine Hamiltonian_replica()
  use DEFS
  use RANDOM  
  use md_defs
  use restart_module
  use SAXS_scoring, only: finalise_SAXS
  implicit none
  
  integer :: i, j, itemp, ierror, vec_length
  integer(8) :: istep, master_istep, n_steps, current_istep, current_step, start_step, end_step
  integer :: ioresult
  real(8) :: en1, en2, T1, T2, ran_number, T, beta, Hi_conf1, Hi_conf2, Hj_conf1, Hj_conf2
  real(8) :: ran3, lambda1, lambda2
  integer :: tN_REPLICA

  logical file_exists
  character(len=2000) :: line
  character(len=100) :: commande, commande1, commande2, fname, fname2, checkname
  character(len=30) :: word,chaine, lpath
  type (t_conformations),dimension(:), allocatable :: conformations
  real(8), dimension(:), allocatable :: vecteur
  real(8), dimension(VECSIZE) :: force1, force2

  ! Prepare the vector to be used in the data exchange between the nodes. Contains the information
  ! of the structure t_conformations
  vec_length = 3*vecsize+6
  allocate(vecteur(vec_length))

  ! Generate the conformation information
  tN_REPLICA = N_E_REPLICA + N_REPLICA

  allocate(conformations(tN_REPLICA))

  do i =1, N_REPLICA
    conformations(i)%id = i
    conformations(i)%temperature = T_replica(i)
    conformations(i)%energyscale = 1.0         ! scales ( = 1.0)
    conformations(i)%score= 0.0
    allocate(conformations(i)%pos(VECSIZE))
    allocate(conformations(i)%posref(VECSIZE))
    allocate(conformations(i)%vel(VECSIZE))
  enddo


  do i =N_REPLICA+1, tN_REPLICA
    conformations(i)%id = i
    conformations(i)%temperature = T_replica(N_REPLICA)    !target_temperature
    conformations(i)%energyscale = E_replica(i-N_REPLICA)    ! scales (from simulator.sh)
    allocate(conformations(i)%pos(VECSIZE))
    allocate(conformations(i)%posref(VECSIZE))
    allocate(conformations(i)%vel(VECSIZE))
  enddo
  
  ! If we have a restart, then we read the main status of the run. 
  ! All the positions and velocities are read in initial_md
  if(restart .eq. 'restart') call read_master_restart(master_istep)
    
  ! Loop for each node, if there is less nodes than replicas, then some will
  ! do two  or more temperatures. The number of replicas divided by the number
  ! of nodes must be an integer. 
  do i=1, tN_REPLICA/ntasks
    T_id = (i-1) *ntasks + taskid +1
    ! Now, check that the subdirectories are there

    call convert_to_chain(T_id,chaine,"0")
    lpath = 'p' // chaine(29:30)

    ! Regrettably, inquire cannot tell whether or not a directory exists,
    ! we must therefore base our information on a file test.existence
    checkname = trim(lpath) //"/test.existence"
    open(unit=83,file=checkname,status='unknown',iostat=ioresult)
    close(83)
    
    if (ioresult .eq. 0 ) then
      file_exists = .true.
    else  
      file_exists = .false.
    endif

    ! If they do not exist, then create them
    if (.not.file_exists) then 
      commande = 'mkdir ' // trim(lpath) // char(0)
#ifdef MP
      call lsystem(commande)  ! Bug on MP requires a different call
#else
      call system(commande)
#endif
    endif
    ! Store the name of the logfile and of the paths
    conformations(T_id)%logfile = trim(lpath)// '/' // trim(MASTER_LOGFILE)  // '.' // chaine(29:30)
    conformations(T_id)%path = trim(lpath) // '/'
  enddo
    
  ! We first initalise and thermalize the various simulations
  do i=1, tN_REPLICA/ntasks
    T_id = (i-1)*ntasks + taskid + 1  
    if (SIMULATION_METHOD .eq. 'MD') then
      call initialise_md(current_istep,conformations(T_id))
      if(restart .ne. 'restart') call thermalize_md(conformations(T_id))
    endif
  end do


  ! Then we prepare the header for the switches
  if (taskid .eq.0) call switch_header()

#ifdef MPI
    ! We now wait here until MD simulations are finished
    call MPI_Barrier(MPI_COMM_WORLD,ierror)
    ! And then, we collect back the information to the central node
    call  collect_info_to_main_node(tN_REPLICA,conformations)
#endif

    ! We then start the simulation per se
    if (SIMULATION_METHOD .eq. 'MD' ) then
      n_steps = n_production / n_step_exchange
    else
      n_steps = 0
    endif

    start_step = 1
    if (restart .eq. 'restart') start_step = master_istep+1
    do istep = start_step, n_steps
      current_step = (istep-1) * n_step_exchange
      end_step = istep * n_step_exchange

      do j=1, tN_REPLICA/ntasks
        T_id = (j-1) * ntasks + taskid + 1
        
        ! If we restart, then if the simulation was waiting for exchange, we do not 
        ! launch a new simulation. Otherwise, we do that.
        if (restart .eq. 'restart' .and. current_istep .gt. master_istep) exit  

        if (SIMULATION_METHOD .eq. 'MD') then
          call run_md(current_step, end_step, conformations(T_id))
        endif
        
        ! We now save the current simulation
        call save_restart(istep,conformations(T_id))
      end do
      
      ! We can reset the restart as we are now continuing our simulation normally
      restart = 'new'

        
#ifdef MPI
      ! We now wait here until MD simulations are finished
      call MPI_Barrier(MPI_COMM_WORLD,ierror)
      ! And then, we collect back the information to the central node
      call  collect_info_to_main_node(tN_REPLICA,conformations)
      call MPI_Barrier(MPI_COMM_WORLD,ierror)
#endif

!################     temperature and hamiltonian replica Exchange      ############
!                            every  n_step_exchange step
!###################################################################################


      if (taskid .eq. 0) then 
        do j=modulo(istep,2)+2, tN_REPLICA,2
          en1 = conformations(j-1)%energy
          en2 = conformations(j)%energy

          T1  = conformations(j-1)%temperature
          T2  = conformations(j)%temperature

          ran_number = ran3()
         
          !==============     temperature replica exchange     ================ 

          if (j .le. N_REPLICA) then
            ! p(T2->T1) = min(1, exp{-[(beta2-beta1)(en1-en2)]}

            if (ran_number .le. exp((1.0d0/T1 - 1.0d0/T2)* (en1-en2))) then
              itemp = conformations(j-1)%id
              conformations(j-1)%id = conformations(j)%id
              conformations(j)%id = itemp
              call switch(vecsize,conformations(j-1)%pos,   conformations(j)%pos)
              call switch(vecsize,conformations(j-1)%posref,conformations(j)%posref)
              call switch(vecsize,conformations(j-1)%vel,   conformations(j)%vel)
          
              ! Renormalise the velocities by the ratio of temperatures after the exchange
              if (SIMULATION_METHOD .eq. 'MD') then    
                conformations(j-1)%vel(:) = sqrt(T1/T2)*conformations(j-1)%vel(:)
                conformations(j  )%vel(:) = sqrt(T2/T1)*conformations(j)%vel(:)
              endif
            endif

          endif
         !=============     Hamiltonian replica exchange      ================== 


          if (j .gt. N_REPLICA) then
            ! p(scale2->scale1) = min(1, exp{-beta [(H1(X1)+H2(X2)) - (H1(X2)+H2(X1)) ]}


            T = conformations(j)%temperature
            beta = 1.0d0/T
            lambda1 = conformations(j-1)%energyscale
            lambda2 = conformations(j)%energyscale

            Hi_conf1 = en1
            Hj_conf2 = en2

            call calcforce(lambda1, conformations(j)%pos,force1,Hi_conf2)
            call calcforce(lambda2, conformations(j-1)%pos,force2,Hj_conf1)

            if (ran_number .le. exp((-beta)*( (Hi_conf1+Hj_conf2)-(Hi_conf2+Hj_conf1) ) )) then
              itemp = conformations(j-1)%id
              conformations(j-1)%id = conformations(j)%id
              conformations(j)%id = itemp
              call switch(vecsize,conformations(j-1)%pos,   conformations(j)%pos)
              call switch(vecsize,conformations(j-1)%posref,conformations(j)%posref)
              call switch(vecsize,conformations(j-1)%vel,   conformations(j)%vel)
        
            endif

          endif

        end do
 

        ! Write the output of the exchange in replica.dat 
        open(unit=FREP,file=REPLICAFILE,status='unknown',action='write',position='append',iostat=ierror)
        open(unit=FLOG,file=MASTER_LOGFILE,status='unknown',action='write',position='append',iostat=ierror)
        write(word,'(i8,i10)') istep,conformations(1)%counter
        line = '  ' // word
        
        do j=1, tN_REPLICA
          write(word,'(i8)') conformations(j)%id
          line = trim(line) //word
        end do
        write(FLOG,"(A)") trim(line)
        write(FREP,"(A)") trim(line)
        close(FREP)
        close(FLOG)
      endif

#ifdef MPI
      ! We now broadcast the information to the other nodes
      call MPI_Barrier(MPI_COMM_WORLD,ierror)
      do i=1, tN_REPLICA
        call conformation_to_vector(vec_length,vecteur,conformations(i))
        vecteur(vec_length) = i
        call MPI_Bcast(vecteur,vec_length,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
        T_id = nint(vecteur(vec_length))
        call vector_to_conformation(vec_length,vecteur,conformations(T_id))
      end do
#endif

      ! The following makes sure that we always have the backup restart files all
      ! at the same time -
      do j=1, tN_REPLICA/ntasks
        T_id = (j-1) * ntasks + taskid + 1
        fname = trim(conformations(T_id)%path) // RESTARTFILE
        commande = 'cp ' // trim(fname)
        fname = trim(conformations(T_id)%path) // 'restart.backup'
        commande1 = trim(commande) // ' ' // trim(fname) // char(0)
        commande = 'cp ' // trim(fname)
        fname2 = trim(fname) // '.1'
        commande2 = trim(commande) // ' ' // trim(fname2) // char(0)

#ifdef MP
        call lsystem(trim(commande2))  ! Bug on MP requires a different call
#else   
        call system(commande2)
#endif

#ifdef MP
        call lsystem(trim(commande1))  ! Bug on MP requires a different call
#else
        call system(commande1)
#endif
 
      end do

      ! And we save the master restart information and update the backups
      do j=1, tN_REPLICA/ntasks
        T_id = (j-1) * ntasks + taskid + 1
        ! We now save the current simulation
        call save_restart(istep,conformations(T_id))
      end do
        
      call save_master_restart(istep)
    end do

  ! lm759 > Finalise SAXS module
  if (compute_saxs_serial) call finalise_SAXS()

  deallocate(vecteur)
  
  return
end subroutine Hamiltonian_replica




!***************************************************************************************
!> @brief Switches two vectors : a -> b and b -> a
!***************************************************************************************
subroutine switch(length,a,b)
  use DEFS
  implicit none
  
  integer, intent(in) :: length
  real(8), dimension(length), intent(inout) :: a, b
  real(8), dimension(:), allocatable :: c
  
  allocate(c(length))
  
  c = a  ! VECTORIAL OPERATION
  a = b  ! VECTORIAL OPERATION
  b = c  ! VECTORIAL OPERATION
  
  deallocate(c)
  
  return
end subroutine switch


subroutine collect_info_to_main_node(nreplica,conformations)
  use DEFS
  implicit none
  
  integer, intent(in) :: nreplica
  type(t_conformations), dimension(nreplica), intent(inout) :: conformations

#ifdef MPI
  integer :: i, isteps, vec_length,ierror
  real(8), dimension(:), allocatable :: vecteur
  real(8), dimension(:,:), allocatable ::  vecteur_in

  integer :: number_replica

  ! Prepare the vector
  vec_length = 3*vecsize+5
  allocate(vecteur(vec_length))

  if (REPLICA_TYPE .eq. 'T_Exchange') then
    number_replica = N_REPLICA
  else
    number_replica = (N_E_REPLICA + N_REPLICA)
  endif


  allocate(vecteur_in(vec_length,number_replica))
  
  ! And then, we collect back the information to the central node
  do isteps=1, (number_replica)/ntasks
  
    ! If this is not the central node, then we prepare the return
    T_id = (isteps-1)*ntasks + taskid+1   ! Define the temperature id
    ! Generates a vector to transmit (faster)
    call conformation_to_vector(vec_length,vecteur,conformations(T_id))
    vecteur(vec_length) = T_id 

    call MPI_Barrier(MPI_COMM_WORLD,ierror)
    call MPI_Gather(vecteur,vec_length,MPI_REAL8,vecteur_in,vec_length,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
    call MPI_Barrier(MPI_COMM_WORLD,ierror)
    if (taskid .eq. 0) then
      ! If this is the central node, then collect from all nodes
      do i=1, ntasks
        T_id = nint(vecteur_in(vec_length,i))  ! This identifies the right set to copy to
        ! Convert the vector to into the rigth information
        call vector_to_conformation(vec_length,vecteur_in(:,i),conformations(T_id))
      end do
    
    endif
  end do
  deallocate(vecteur_in)  
  deallocate(vecteur)
#endif
  return
end subroutine collect_info_to_main_node

subroutine conformation_to_vector(vec_length,vecteur,conformation)
  use DEFS
  implicit none

  integer, intent(in) :: vec_length
  type(t_conformations), intent(in)  :: conformation
  real(8), dimension(vec_length), intent(out) :: vecteur

  vecteur(          1:  vecsize) = conformation%pos(1:vecsize)
  vecteur(  vecsize+1:2*vecsize) = conformation%vel(1:vecsize)
  vecteur(2*vecsize+1:3*vecsize) = conformation%posref(1:vecsize)
  vecteur(3*vecsize+1) = conformation%energy
  vecteur(3*vecsize+2) = conformation%temperature
  vecteur(3*vecsize+3) = conformation%id
  vecteur(3*vecsize+4) = conformation%energyscale

  return
end subroutine conformation_to_vector

subroutine vector_to_conformation(vec_length,vecteur,conformation)
  use DEFS
  implicit none

  integer, intent(in) :: vec_length
  type(t_conformations), intent(inout) :: conformation
  real(8), dimension(vec_length), intent(in)    :: vecteur

  conformation%pos         = vecteur(          1:  vecsize)
  conformation%vel         = vecteur(  vecsize+1:2*vecsize)
  conformation%posref      = vecteur(2*vecsize+1:3*vecsize)
  conformation%energy      = vecteur(3*vecsize+1)
  conformation%temperature = vecteur(3*vecsize+2) 
  conformation%id          = nint(vecteur(3*vecsize+3))
  conformation%energyscale = vecteur(3*vecsize+4) 

  return
end subroutine vector_to_conformation

!******************************************************************************
!*
!> @brief After the thermalization (or restart), we write the headers in the various
!> file that will be used
!*
!******************************************************************************
subroutine switch_header()
  use DEFS
  implicit none

  integer :: i,ierror, number_replica
  character(len=2000) :: line
  character(len=30) :: word,chaine, lpath

  ! Open the necessary files. If replicafile already exists (it is a restart), then we do not 
  ! had a header since it has already one. 
  open(unit=FLOG,file=MASTER_LOGFILE,status='unknown',action='write',position='append',iostat=ierror)
  if (restart .ne. 'restart') then
    open(unit=FREP,file=REPLICAFILE,status='unknown',action='write',position='rewind',iostat=ierror)
  endif
  
  ! Create the files and add the header for the replicas.
  line = '      istep' //'   counter' 


  if (REPLICA_TYPE .eq. 'T_Exchange') then
    number_replica = N_REPLICA
  else
    number_replica = (N_E_REPLICA + N_REPLICA)
  endif

  do i=1, number_replica
    call convert_to_chain(i,chaine,"0")
    lpath = 'R' // chaine(29:30)
    write(word,'(A8)') trim(lpath)
    line = trim(line) // word
  end do
  write(flog,'(A)') trim(line)
  if (restart .ne. 'restart') write(frep,'(A)') trim(line)

  line = '     ******  ********'
  do i=1, number_replica
    line = trim(line) // '  ******'
  end do
  write(flog,'(A)') trim(line)
  if (restart .ne. 'restart') write(frep,'(A)') trim(line)


  do i=1, N_TEMP
    call convert_to_chain(i,chaine,"0")
    lpath = 'R' // chaine(29:30)
    write(word,'(A8)') trim(lpath)
    line = trim(line) // word
  end do
  write(flog,'(A)') trim(line)
  if (restart .ne. 'restart') write(frep,'(A)') trim(line)

  line = '     ******  ********'
  do i=1, N_TEMP
    line = trim(line) // '  ******'
  end do
  write(flog,'(A)') trim(line)
  if (restart .ne. 'restart') write(frep,'(A)') trim(line)
  close(flog)
  if (restart .ne. 'restart') close(frep)

end subroutine switch_header

!***************************************************************
! > @brief      Simulated Tempering                                    !
!***************************************************************

subroutine site()
  use DEFS
  use RANDOM
  use md_defs
  use restart_module
  use SAXS_scoring, only: finalise_SAXS
  implicit none

  integer :: i, j, ierror, vec_length, current_istep
  integer(8) :: istep, master_istep, n_steps, current_step, start_step, end_step
  integer :: ioresult
  real(8) :: en1, T1, T2, ran_number,F1,F2,F0,T0,B1,B2,E1,E2
  real(8) :: ran3,DELTAP,DELTAM,ran_numberp
  integer :: exch(N_TEMP,N_TEMP)
  double precision::energy(N_TEMP),norm(N_TEMP), pdiv

  logical file_exists
  character(len=2000) :: line
  character(len=100) :: commande, checkname
  character(len=30) :: word,chaine, lpath
  type (t_conformations),dimension(:), allocatable :: conformations
  real(8), dimension(:), allocatable :: vecteur

  ! Prepare the vector to be used in the data exchange between the nodes. Contains the information
  ! of the structure t_conformations
  vec_length = 3*vecsize+5
  allocate(vecteur(vec_length))

  ! Generate the conformation information
  allocate(conformations(N_TEMP))
  do i =1, N_TEMP
    conformations(i)%id = i
    conformations(i)%temperature = T_tempering(i)
    conformations(i)%free_energy = F_tempering(i)
    conformations(i)%E_energy    = E_tempering(i)
    conformations(i)%N_energy    = N_tempering(i)
    conformations(i)%energyscale = 1.0         ! scales ( = 1.0)
    allocate(conformations(i)%pos(VECSIZE))
    allocate(conformations(i)%posref(VECSIZE))
    allocate(conformations(i)%vel(VECSIZE))
  enddo

  ! If we have a restart, then we read the main status of the run.
  ! All the positions and velocities are read in initial_md
  if(restart .eq. 'restart') call read_master_restart(master_istep)

  ! Loop for each node, if there is less nodes than replicas, then some will
  ! do two  or more temperatures. The number of replicas divided by the number
  ! of nodes must be an integer.

  do i=1, N_TEMP
    T_id = i

    ! Now, check that the subdirectories are there
    call convert_to_chain(T_id,chaine,"0")
    lpath = 'p' // chaine(29:30)

    ! Regrettably, inquire cannot tell whether or not a directory exists,
    ! we must therefore base our information on a file test.existence
    checkname = trim(lpath) //"/test.existence"
    open(unit=83,file=checkname,status='unknown',iostat=ioresult)
    close(83)

    if (ioresult .eq. 0 ) then
      file_exists = .true.
    else
      file_exists = .false.
    endif

    ! If they do not exist, then create them
    if (.not.file_exists) then
      commande = 'mkdir ' // trim(lpath) // char(0)
#ifdef MP
      call lsystem(commande)  ! Bug on MP requires a different call
#else
      call system(commande)
#endif
    endif
    ! Store the name of the logfile and of the paths
    conformations(T_id)%logfile = MASTER_LOGFILE
    conformations(T_id)%path = trim(lpath) // '/'
  enddo

!!!!! END PREPARATION: distribute temperatures, prepare files !!!!!!!!!!
!!!!! We then initalise and thermalize the simulation at temperature T_id = 1 on node 0

!    if (taskid .eq. 0) T_id = taskid + 1

!    do i=1,N_TEMP

    T_id = last_id

    if (SIMULATION_METHOD .eq. 'MD') then
      call initialise_md(current_istep,conformations(T_id))
      if(restart .ne. 'restart') call thermalize_md(conformations(T_id))
    endif

!    end do

!    do i=2,N_REPLICA
!    conformations(i)%vel(:) = 0.0
!    conformations(i)%pos(:) = 0.0
!    conformations(i)%posref(:) = 0.0
!    conformations(i)%energy    = 0.0
!    enddo

! Then we prepare the header for the switches

 if (taskid .eq. 0) call switch_header()

   
!!!!!!!!!!!!!!Finishing equilibration process at all Temperature !!!!!!!!!

! We then start the simulation per se

    if (SIMULATION_METHOD .eq. 'MD') then
      n_steps = n_production / n_step_exchange
    else
      n_steps = 0
    endif

    start_step = 1
!    if (restart .eq. 'restart') start_step = master_istep+1


!!!!!!!!!!!!!!!! Enter MAIN LOOP !!!!!!!!!!!!!!!!!!!!!!!
!! We start with the lowest temperature!!!!!!!!!!!!!!!!!

!     norm(1)               = 1.0

    do i=1,N_TEMP

        energy(i)               =  E_tempering(i)
        norm(i)                 =  N_tempering(i)

    enddo


    T_id = last_id

!  do i =1, N_TEMP
!  write(*,*) T_tempering(i), F_tempering(i), E_tempering(i),N_tempering(i),T_id
!  enddo
!  write(*,*) start_step,n_steps
 
    open(unit=FREP,file=REPLICAFILE,status='unknown',action='write',position='append',iostat=ierror)
    open (unit=6, file='temp_traj.txt', status='unknown',action='write',position='append',iostat=ierror)

    start_step = 1
    if (restart .eq. 'restart') start_step = master_istep+1
    do istep = start_step, n_steps
        current_step = (istep-1) * n_step_exchange
        end_step = istep * n_step_exchange


        write(word,'(i8,i10)') istep, T_id
        line = '  ' // word

        do j=1, N_TEMP
          write(word,'(f12.3)') conformations(j)%free_energy
          line = trim(line) //word
        end do
        write(FREP,"(A)") trim(line)


        if(mod(istep,10).eq.0) then
            write(6,*) istep,T_id,conformations(T_id)%energy,conformations(T_id)%temperature
        endif


        if(taskid.eq.0) then   !!!! Run on node 0

        ! If we restart, then if the simulation was waiting for exchange, we do not
        ! launch a new simulation. Otherwise, we do that.

! PHUONG

        if (restart .eq. 'restart' .and. current_istep .gt. master_istep) exit

        if (SIMULATION_METHOD .eq. 'MD') then
          call run_md(current_step, end_step, conformations(T_id))
        endif


        energy(T_id) = energy(T_id) + conformations(T_id)%energy
        norm(T_id)   = norm(T_id)   + 1.0

        if (mod(int(sum(norm)), N_TEMP*500) .eq. 0) then
          pdiv = 0
          do i = 1,N_TEMP
            pdiv = pdiv + abs(norm(i)/sum(norm) - 1.0d0/float(N_temp))
          enddo
          pdiv = pdiv/2000.0
          energy = energy*(1-pdiv)
          norm = norm*(1-pdiv)
        endif

        ! We now save the current simulation
        call save_restart(istep,conformations(T_id))

        endif

!! test here
        
!        if(istep.eq.start_step) then
!        open (unit=10, file='start_rerun.dat', status='unknown',action='write',position='append',iostat=ierror)
!        write(10,*) T_id,conformations(T_id)%energy,conformations(T_id)%temperature
!        endif
         if (mod(int(sum(norm)), N_TEMP*50) .eq. 0) then
         open (unit=9, file='input_rerun.txt', status='replace',action='write',position='rewind',iostat=ierror)
         write(9,*)'setenv Free_energy_tempering "',(conformations(j)%free_energy,j=1,N_TEMP), '"'
         write(9,*)'setenv Last_accum_energy "',(energy(j),j=1,N_TEMP), '"'
         write(9,*)'setenv Last_accum_norm "',(norm(j),j=1,N_TEMP), '"'
         write(9,*)'setenv Last_id',T_id
         close(9)
         endif
!!!!!!!!Calculate free energy!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!        conformations(1)%free_energy  = 0.0d0

        

        do i=2,N_TEMP

        T1 = conformations(i-1)%temperature        
        T2 = conformations(i)%temperature
        B1 = 1.0d0/T1
        B2 = 1.0d0/T2
        if(energy(i-1).eq.0.0) then
            E1 = 0.0
        else
            E1 = energy(i-1)/norm(i-1)
        endif

        if(energy(i).eq.0.0) then
            E2 = 0.0
        else
            E2 = energy(i)/norm(i)
        endif

        F1 = conformations(i-1)%free_energy
        F2 = F1 + (B2 - B1)*(E2 + E1)/2.0d0

        conformations(i)%free_energy = F2

        enddo

!        do i=2,N_TEMP
!        write(*,*) conformations(i)%free_energy
!        enddo

!         write(*,*) T_id,conformations(T_id)%energy,conformations(T_id)%temperature

!!!!!!!!! END n_step_exchange steps at temperature T_id = 1!!!!!!

      ! We can reset the restart as we are now continuing our simulation normally

      restart = 'new'

!      if (taskid .eq. T_id) then

        j = T_id

        en1 = conformations(j)%energy

        T1  = conformations(j)%temperature   !! This is the one we ran MD
        F1  = conformations(j)%free_energy

        if(T_id.lt.N_TEMP) then
          T2  = conformations(j+1)%temperature
          F2  = conformations(j+1)%free_energy
        endif

        if(T_id.gt.1) then       
          T0  = conformations(j-1)%temperature
          F0  = conformations(j-1)%free_energy
        endif

        ran_number  = ran3()
        ran_numberp = ran3()


        DELTAP = (1.0d0/T2 - 1.0d0/T1)*en1 - (F2 - F1)
        DELTAM = (1.0d0/T0 - 1.0d0/T1)*en1 - (F0 - F1)
  
!        write(*,*) DELTAP,DELTAM
!        write(*,*) ran_number, ran_numberp
!        write(*,*) ""

          if( (ran_numberp.lt.0.50).and.(ran_number .le. exp(-DELTAM)).and.(j.gt.1) ) then

            conformations(j-1)%id         = conformations(j)%id       ! tested: not important
            conformations(j-1)%energy     = conformations(j)%energy   ! tested: not important
            conformations(j-1)%logfile    = conformations(j)%logfile  ! tested: not important

            conformations(j-1)%pos(:)     = conformations(j)%pos(:)
            conformations(j-1)%posref(:)  = conformations(j)%posref(:)

            conformations(j-1)%vel(:)     = sqrt(T0/T1)*conformations(j)%vel(:)

            T_id = T_id - 1
            exch(j,j-1) = 1
          else if(j.gt.1) then
            T_id = T_id
            exch(j,j-1) = 0
          endif  ! Accept down move

          if( (ran_numberp.gt.0.50).and.(ran_number .le. exp(-DELTAP)).and.(j.lt.N_TEMP) ) then

               conformations(j+1)%id         = conformations(j)%id
               conformations(j+1)%energy     = conformations(j)%energy
               conformations(j+1)%logfile    = conformations(j)%logfile


               conformations(j+1)%pos(:)     = conformations(j)%pos(:)
               conformations(j+1)%posref(:)  = conformations(j)%posref(:)

               conformations(j+1)%vel(:)     = sqrt(T2/T1)*conformations(j)%vel(:)


              T_id = T_id + 1
              exch(j,j+1) = 1            
          else if(j.lt.N_TEMP) then
            write(*,*) T_id, N_TEMP, j, ran_numberp, ran_number
              T_id = T_id
              exch(j,j+1) = 0
           endif  ! Accept up move


          ! Write the output of the exchange in replica.dat
!        open(unit=FREP,file=REPLICAFILE,status='unknown',action='write',position='append',iostat=ierror)
!        open(unit=FLOG,file=MASTER_LOGFILE,status='unknown',action='write',position='append',iostat=ierror)
!        write(word,'(i8,i10)') istep,conformations(1)%counter
!        line = '  ' // word

!        do j=1, N_REPLICA
!          write(word,'(i8)') conformations(j)%id
!          line = trim(line) //word
!        end do
!        write(FLOG,"(A)") trim(line)
!        write(FREP,"(A)") trim(line)
!        close(FREP)
!        close(FLOG)


!        endif ! taskid = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! And we save the master restart information and update the backups
        ! We now save the current simulation

        call save_restart(istep,conformations(T_id))
  
        call save_master_restart(istep)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    end do   !    Production run

    close(FREP)
    close(6)


    ! lm759 > Finalise SAXS module
    if (compute_saxs_serial) call finalise_SAXS()

    deallocate(vecteur)

    return
end subroutine site


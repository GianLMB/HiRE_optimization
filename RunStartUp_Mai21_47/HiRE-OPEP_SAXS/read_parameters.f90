!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Reads the parameter file 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine definitions()
  use defs
  use RANDOM
  use md_defs
  use SAXS_scoring
  implicit none
  integer :: i,ierror, ntemps, nenergies
  character(8) :: date
  character(10) :: time
  character(5) :: zone
  integer, dimension(8) :: value
  real(8), dimension(:), allocatable :: dum_T_replica
  real(8), dimension(:), allocatable :: dum_E_replica
  real(8), dimension(:), allocatable :: dum_T_tempering
  real(8), dimension(:), allocatable :: dum_F_tempering
  real(8), dimension(:), allocatable :: dum_E_tempering
  real(8), dimension(:), allocatable :: dum_N_tempering

  logical :: exists_already

  character(len=100) :: fname
  character(len=150) :: commande
  character(len=9)   :: digit = "123456789"
  character(len=12)  :: restrict_motion
  character(len=50) dummy
  character(len=2000) dummy2
  character(len=2000) dummy3
  character(len=2000) dummy4
  character(len=2000) dummy5
  character(len=2000) dummy6

  integer, parameter :: fchain = 17

  character(len=20) :: molecule_type
  COMMON /MOLECULETYPE/ molecule_type 

  call GET_ENVIRONMENT_VARIABLE('MOLECULE_TYPE', dummy)
  if (dummy .eq. '') then
  write(*,*) 'Error: variable MOLECULE_TYPE is not defined'
  stop
  else
     read(dummy,*) molecule_type
     molecule_type = trim(molecule_type)
     if(molecule_type .NE. "RNA" .AND. molecule_type .NE. "DNA" .and. molecule_type .ne. "protein") then
     write(*,*) 'Error: variable MOLECULE_TYPE has to be RNA, DNA, or protein'
     stop
     endif
  endif
 
  call GET_ENVIRONMENT_VARIABLE('Simulation_Type', dummy)

  if (dummy .eq. '') then
     SIMULATION_TYPE = 'serial'
  else
     read(dummy,*) SIMULATION_TYPE
  endif

  if (SIMULATION_TYPE .eq. 'serial') then
     call GET_ENVIRONMENT_VARIABLE('Temperature', dummy)
     if (dummy .eq. '') then   
        write(*,*) 'Error: variable TEMPERATURE is not defined'
        stop
     else
        read(dummy,*) target_temperature
        ! Conversion from Kelvin to kcal/mol
        target_temperature = target_temperature * 0.59227 / 298
     endif
  endif

  call GET_ENVIRONMENT_VARIABLE('debug_status',dummy)
  if (dummy .eq. '') then
     debug_status = 'normal'
  else
     read(dummy,*) debug_status
  endif

  call GET_ENVIRONMENT_VARIABLE('use_qbug',dummy)
  if (dummy .eq. '') then
     use_qbug = .false.
  else
     read(dummy,*) dummy2
     use_qbug = dummy2 == 'yes'
  endif


!!!=============Modification for INMA


  call GET_ENVIRONMENT_VARIABLE('Internal_Coordinates',dummy)
  
  if (dummy .eq. '') then
     use_ics = .false.
  else
      read(dummy,*) use_ics
  endif
  
  

  if(use_ics)then
     MINIMIZATION_TYPE='MIN_steep'
     
     call GET_ENVIRONMENT_VARIABLE('Backbone_Variables',dummy)
     if (dummy .eq. '') then
        use_back = .false.
     else
        read(dummy,*) use_back
     endif
     call GET_ENVIRONMENT_VARIABLE('Base_Sugar_Variables',dummy)
     if (dummy .eq. '') then
        use_base = .false.
     else
        read(dummy,*) use_base
     endif
     call GET_ENVIRONMENT_VARIABLE('Valence_Angle_Variables',dummy)
     if (dummy .eq. '') then
        use_val = .false.
     else
        read(dummy,*) use_val
     endif
  else
     call GET_ENVIRONMENT_VARIABLE('MINIMIZATION_TYPE',dummy)
     if (dummy .eq. '') then
        MINIMIZATION_TYPE='MIN_damped'
     else
        read(dummy,*) MINIMIZATION_TYPE
     endif
  endif
  
  

  call GET_ENVIRONMENT_VARIABLE('Normal_Modes',dummy)
  if (dummy .eq. '') then
     use_nma = .false.
  else
     read(dummy,*) use_nma
  endif

 !!==============
  

  call GET_ENVIRONMENT_VARIABLE('Titration',dummy)
  if (dummy .eq. '') then
     use_tit = .false.
  else
     read(dummy,*) use_tit
  endif
  
    call GET_ENVIRONMENT_VARIABLE('n_steps_tit',dummy)
  if (dummy .eq. '') then
     n_steps_tit = 1
  else
     read(dummy,*) n_steps_tit
  endif
  

  call GET_ENVIRONMENT_VARIABLE('Restart_Run',dummy)
  if (dummy .eq. '') then
     restart = 'new'
  else
     read(dummy,*) restart
  endif

  call GET_ENVIRONMENT_VARIABLE('Potential_Scaling_Factor',dummy)
  if (dummy .eq. '') then
     force_scaling_factor = 1.0d0;
  else
     read(dummy,*) force_scaling_factor
  endif

  ! periodic boundary condition
  call GET_ENVIRONMENT_VARIABLE('Periodic_Boundary_Condition',dummy)
  if (dummy .eq. '') then
     PBC = .false.
  else
     read(dummy,*) PBC
  endif

  ! If perdiodic boundary conditions, we must define a Box Length, using for
  ! periodic boundary condition
  if (PBC) then 
     call GET_ENVIRONMENT_VARIABLE('Box_Length',dummy)
     if (dummy .eq. '') then
        write(*,*) 'Error: variable Box_Length is not defined but ', &
             'periodic boundary conditions is true.'
        stop
     else
        read(dummy,*) BL
     endif
  endif

  ! The center of mass as a reference for periodic boundary condition
  ! it is only applied while writing the pdb file, to prevent having
  ! broken chains.

  call GET_ENVIRONMENT_VARIABLE('PDB_center_of_mass',dummy)
  if (dummy .eq. '') then
     C_M = .false.
  else
     read(dummy,*) C_M
  endif

  ! Read the number of atoms. 
  call GET_ENVIRONMENT_VARIABLE('NATOMS', dummy)
  if (dummy .eq. '') then   
     write(*,*) 'Error: variable NATOMS is not defined'
  else
     read(dummy,*) NATOMS
  endif

  ! Set the random number generator
  call GET_ENVIRONMENT_VARIABLE('RANDOM_SEED', dummy)
  if (dummy .eq. '') then
     idum=0
  else
     read(dummy,*) idum
  endif

  ! If idum equal to  zero, then we use the clock for the random number
  if (idum .eq. 0 ) then
     call date_and_time(date,time,zone,value)
     idum = -1 * mod( (1000 * value(7) + value(8)), 1024)
  endif

  ! Read the name of the log file which keeps tracks of what is happening
  call GET_ENVIRONMENT_VARIABLE('LOGFILE', dummy)
  if (dummy .eq. '') then
     MASTER_LOGFILE   = 'log.file'
  else
     read(dummy,*) MASTER_LOGFILE
  endif

  ! Name of the file where the replica events are stored
  call GET_ENVIRONMENT_VARIABLE('REPLICAFILE',dummy)
  if (dummy .eq. '') then
     REPLICAFILE = 'replicas.dat'
  else
     read(dummy,*) REPLICAFILE
  endif

  ! Number of digits used for naming the configuration files when each 
  ! configuration is stored separately. 
  call GET_ENVIRONMENT_VARIABLE('NDIGITS',dummy)
  if (dummy .eq. '') then
     ndigits_filename = 5
  else
     read(dummy,*) ndigits_filename
  endif

  ! Name of the file where the chain information is stored
  call GET_ENVIRONMENT_VARIABLE('chain_file', dummy)
  if (dummy .eq. '') then
     chain_file = 'ichain.dat'
  else
     read(dummy,*) chain_file
  endif

  ! Use xtc (default is pdb
  call GET_ENVIRONMENT_VARIABLE('use_xtc',dummy)
  if (dummy .eq. '') then
     usextc = .false.
  else 
     read(dummy,*) usextc
  endif

  ! Store the configurations file in seperate or a single file per temperature
  call GET_ENVIRONMENT_VARIABLE('single_file',dummy)
  if (dummy .eq. '') then
     singlefile = .false.
  else 
     read(dummy,*) singlefile
  endif

  ! if we need to different initial configurations for different temperatures (in replica exchange)
  ! we have to set init_single_file to .false. in simulateur.sh (read changes.txt) for more information

  call GET_ENVIRONMENT_VARIABLE('init_singlefile',dummy)
  if (dummy .eq. '') then
     init_single_file = .true.
  else 
     read(dummy,*) init_single_file
  endif

  ! Restrict the motion - if yes, it is possible to add springs only during thermalization (restraint)
  ! or during the whole simulation (constraint). The reference position is defined on the Ca after 
  ! the initial minimization. 
  call GET_ENVIRONMENT_VARIABLE('Restrict_Motion',dummy)
  if (dummy .eq. '') then
     restrict_motion = 'none'
     constrained_fragments = .false.
     restrained_fragments = .false.
  else
     read(dummy, *) restrict_motion
     if (restrict_motion .eq. 'constraint') then
        constrained_fragments = .true.
        restrained_fragments = .false.
     else if (restrict_motion .eq. 'restraint') then
        constrained_fragments = .false.
        restrained_fragments = .true.
     else
        constrained_fragments = .false.
        restrained_fragments = .false.
     end if
  endif

  ! Defines the spring constant used if motion is restricted. 
  if (constrained_fragments.or.restrained_fragments) then
     call GET_ENVIRONMENT_VARIABLE('Spring_Constant', dummy)
     if (dummy .eq. '' ) then
        k_spring = 10.0
     else
        read(dummy,*) k_spring
     endif
  endif


  ! We now read the number of fragments
  inquire(file=chain_file,exist=exists_already)
  if (exists_already) then
     open(unit=fchain,file=chain_file,status='unknown',action='read')
     read(fchain,*) nfrag

     if (constrained_fragments.or.restrained_fragments) then
        allocate(list_fragments(nfrag,3))  
     else
        allocate(list_fragments(nfrag,2))  
     endif

     ! The constrained fragments are indicated in the ichain.dat file. 
     do i=1, nfrag
        if (constrained_fragments.or.restrained_fragments) then
           read(fchain,*) list_fragments(i,1),list_fragments(i,2),list_fragments(i,3)
        else
           read(fchain,*) list_fragments(i,1),list_fragments(i,2)
        endif
     enddo
     close(fchain)
  else
     write(6,"('ERROR: File ', A, ' does not exist.')") chain_file
     stop
  endif

  VECSIZE = 3 * NATOMS
  VECSIZE1 =  VECSIZE / NFRAG

  ! We first check whether the file "LOGFILE" exists. If so, we copy it before
  ! we start.
  if (taskid .eq. 0 ) then
     fname = MASTER_LOGFILE
     do i=1, 9
        inquire(file=fname,exist=exists_already)
        if (exists_already)  then
           fname = trim(MASTER_LOGFILE) // "." // digit(i:i)
        else
           if (i .gt. 1 ) then
              commande = "mv " // MASTER_LOGFILE // "  " // fname
#ifdef MP
              call lsystem(commande)  ! Bug on MP requires a different call
#else
              call system(commande)
#endif
           endif
           exit
        endif
     end do
  endif

!!!!!!!!!!!!!!!!!!!!!!READ TEMPERATURE 

  ! IF the simulation is a simulated tempering

if (SIMULATION_TYPE .eq. 'tempering') then
     call GET_ENVIRONMENT_VARIABLE('Number_Tempering', dummy)
     if (dummy .eq. '') then
        write(6,"('Error: must give the number of temperatures - Number_Tempering')")
     else
        read(dummy,*) N_TEMP
     endif

   call GET_ENVIRONMENT_VARIABLE('Temperatures_tempering',dummy3)
     if (dummy3 .eq. '') then
        write(6,"('Error: Must indicate the temperatures for the tempering')")
        stop
     else
        ! Because of the length of the input, the environment variable might be
        ! defined on many lines. Some compilers have problems with that so we
        ! have to first read the data, including errors which are transformed into
        ! zeros then remove them.
        allocate(T_tempering(N_TEMP))
        allocate(dum_T_tempering(2*N_TEMP))
        dum_T_tempering = 0.0d0
        read(dummy3,*,end=444) dum_T_tempering
444     ntemps = 0
        do i=1, 2*N_TEMP
           if ( (dum_T_tempering(i) .gt. 0.0d0)  .and. (ntemps .lt. N_TEMP) ) then
              ntemps = ntemps + 1
              T_tempering(ntemps) = dum_T_tempering(i)
           endif
        end do
        do i=1, N_TEMP
           T_tempering(i) = T_tempering(i) * 0.59227 / 298
        end do
     endif

!!!!!!!!!!!!! READ FREE ENERGY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   call GET_ENVIRONMENT_VARIABLE('Free_energy_tempering',dummy4)
     if (dummy4 .eq. '') then
        write(6,"('Error: Must indicate the free energy for the tempering')")
        stop
     else
        ! Because of the length of the input, the environment variable might be
        ! defined on many lines. Some compilers have problems with that so we
        ! have to first read the data, including errors which are transformed into
        ! zeros then remove them.
        allocate(F_tempering(N_TEMP))
        allocate(dum_F_tempering(2*N_TEMP))
        dum_F_tempering = -10000.0d0
        read(dummy4,*,end=555) dum_F_tempering
555     ntemps = 0
        do i=1, 2*N_TEMP   
           if ( (ntemps .lt. N_TEMP) ) then
              ntemps = ntemps + 1
              F_tempering(ntemps) = dum_F_tempering(i)
           endif
        end do
        do i=1, N_TEMP   
           F_tempering(i) = F_tempering(i)
        end do
     endif

!!!! read files for restart ST

   call GET_ENVIRONMENT_VARIABLE('Last_accum_energy',dummy5)
     if (dummy5 .eq. '') then
        write(6,"('Error: Must indicate the last accumulated energy')")
        stop
     else
        ! Because of the length of the input, the environment variable might be
        ! defined on many lines. Some compilers have problems with that so we
        ! have to first read the data, including errors which are transformed into
        ! zeros then remove them.
        allocate(E_tempering(N_TEMP))
        allocate(dum_E_tempering(2*N_TEMP))
        dum_E_tempering = -10000.0d0
        read(dummy5,*,end=666) dum_E_tempering
666     ntemps = 0
        do i=1, 2*N_TEMP
           if ( (ntemps .lt. N_TEMP) ) then
              ntemps = ntemps + 1
              E_tempering(ntemps) = dum_E_tempering(i)
           endif
        end do
        do i=1, N_TEMP
           E_tempering(i) = E_tempering(i)
        end do
     endif


   call GET_ENVIRONMENT_VARIABLE('Last_accum_norm',dummy6)
     if (dummy6 .eq. '') then
        write(6,"('Error: Must indicate the last norm')")
        stop
     else
        ! Because of the length of the input, the environment variable might be
        ! defined on many lines. Some compilers have problems with that so we
        ! have to first read the data, including errors which are transformed into
        ! zeros then remove them.
        allocate(N_tempering(N_TEMP))
        allocate(dum_N_tempering(2*N_TEMP))
        dum_N_tempering = -10000.0d0
        read(dummy6,*,end=777) dum_N_tempering
777     ntemps = 0
        do i=1, 2*N_TEMP
           if ( (ntemps .lt. N_TEMP) ) then
              ntemps = ntemps + 1
              N_tempering(ntemps) = dum_N_tempering(i)
           endif
        end do
        do i=1, N_TEMP
           N_tempering(i) = N_tempering(i)
        end do
     endif

     call GET_ENVIRONMENT_VARIABLE('Last_id', dummy)
     if (dummy .eq. '') then
        write(6,"('ERROR: Must indicate the last id')")
        stop
     else
        read(dummy,*) last_id
     endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     call GET_ENVIRONMENT_VARIABLE('Exchange_every_n_steps', dummy)
     if (dummy .eq. '') then
        write(6,"('ERROR: Must indicate the number of steps between exchanges')")
        stop
     else
        read(dummy,*) n_step_exchange
     endif

     endif

    


  ! IF the simulation is a replica exchange, then we set up the replicas
  if (SIMULATION_TYPE .eq. 'replica' .or. SIMULATION_TYPE .eq. 'saxs') then
     call GET_ENVIRONMENT_VARIABLE('Number_Replica', dummy)
     if (dummy .eq. '') then
        write(6,"('Error: must give the number of replica - Number_Replica')")
     else
        read(dummy,*) N_REPLICA
     endif

     call GET_ENVIRONMENT_VARIABLE('Number_E_Replica', dummy)
     if (dummy .eq. '') then
        REPLICA_TYPE = 'T_Exchange'
     else
        read(dummy,*) N_E_REPLICA
        if (N_E_REPLICA .eq. 0) then 
           REPLICA_TYPE = 'T_Exchange'  !     if there is no N_E_REPLICA or the number is set to "0" 
                                        !     we have T_exchange (normal replica exchange)
        else
           REPLICA_TYPE = 'E_Scale'     !     else, we have hamiltonian replica exchange
        endif 
     endif

     call GET_ENVIRONMENT_VARIABLE('Temperatures_Replica',dummy2)
     if (dummy2 .eq. '') then
        write(6,"('Error: Must indicate the temperatures of the replicas')")
        stop
     else 
        ! Because of the length of the input, the environment variable might be
        ! defined on many lines. Some compilers have problems with that so we
        ! have to first read the data, including errors which are transformed into
        ! zeros then remove them. 
        allocate(T_replica(N_REPLICA))
        allocate(dum_T_replica(2*N_REPLICA))
        dum_T_replica = 0.0d0
        read(dummy2,*,end=222) dum_T_replica
222     ntemps = 0
        do i=1, 2*N_REPLICA
           if ( (dum_T_replica(i) .gt. 0.0d0)  .and. (ntemps .lt. N_REPLICA) ) then
              ntemps = ntemps + 1
              T_replica(ntemps) = dum_T_replica(i)
           endif
        end do
        do i=1, N_REPLICA
           T_replica(i) = T_replica(i) * 0.59227 / 298
        end do
     endif

     call GET_ENVIRONMENT_VARIABLE('Exchange_every_n_steps', dummy)
     if (dummy .eq. '') then
        write(6,"('ERROR: Must indicate the number of steps between exchanges')")
        stop
     else
        read(dummy,*) n_step_exchange
     endif

  endif

  ! IF the simulation is a Hamiltonian replica exchange, then we set up the replicas
  if (SIMULATION_TYPE .eq. 'replica' .and. REPLICA_TYPE .eq. 'E_Scale') then

     call GET_ENVIRONMENT_VARIABLE('Scales_Replica',dummy2)
     if (dummy2 .eq. '') then
        write(6,"('Error: Must indicate the scales of the Energy replica')")
        stop
     else 
        ! Because of the length of the input, the environment variable might be
        ! defined on many lines. Some compilers have problems with that so we
        ! have to first read the data, including errors which are transformed into
        ! zeros then remove them. 
        allocate(E_replica(N_E_REPLICA))
        allocate(dum_E_replica(2*N_E_REPLICA))
        dum_E_replica = 0.0d0
        read(dummy2,*,end=333) dum_E_replica
333     nenergies = 0

        do i=1, 2*N_E_REPLICA
           if ( (dum_E_replica(i) .gt. 0.0d0)  .and. (nenergies .lt. N_E_REPLICA) ) then
              nenergies = nenergies + 1
              E_replica(nenergies) = dum_E_replica(i)
           endif
        end do

     endif

  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!  SAXS  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (SIMULATION_TYPE .eq. 'saxs') then

    call getenv('Number_Replica', dummy)
    if (dummy .eq. '') then
      write(6,"('Error: must give the number of replica - Number_Replica')")
    else
      read(dummy,*) N_REPLICA
    end if

    call getenv('SAXS_Calculation_Type', dummy)
    if (dummy .eq. '') then
      in_solution_curve = .true.
    else if (dummy .eq. 'vacuum') then
      in_solution_curve = .false.
    else if (dummy .eq. 'solution') then
      in_solution_curve = .true.
    else
      STOP 'Unrecognized SAXS Calculation Type'
    end if

    if (in_solution_curve) then

      call getenv('SAXS_Refine_Hydration', dummy)
      if (dummy .eq. '') then
        refine_hydration_layer = .false.
      else
        read(dummy,*) refine_hydration_layer
      end if

      if (refine_hydration_layer) then

        call getenv('SAXS_Dummy_Water_Radius', dummy)
        if (dummy .eq. '') then
          DX = 7.5
        else
          read(dummy,*) DX
        end if

        call getenv('SAXS_Dummy_Water_Layers', dummy)
        if (dummy .eq. '') then
          n_shells = 1
        else
          read(dummy,*) n_shells
        end if

        call getenv('SAXS_w_shell', dummy)
        if (dummy .eq. '') then
          SAXS_W_shell = 0.3d0
        else
          read(dummy,*) SAXS_w_shell
        end if

      end if

    end if

    call getenv('SAXS_Norm_Type', dummy)
    if (dummy .eq. '') then
      SAXS_norm_type = 2
    else
      read(dummy,*) SAXS_norm_type
    endif

    call getenv('SAXS_Vect_Max', dummy)
    if (dummy .eq. '') then
      SAXS_vect_max = 1.0d0
    else
      read(dummy,*) SAXS_vect_max
    endif

    call getenv('SAXS_Profiling_Rate', dummy)
    if (dummy .eq. '') then
      n_rate_saxs = 2
    else
      read(dummy,*) n_rate_saxs
    endif

    call getenv('SAXS_Energy_Coefficient', dummy)
    if (dummy .eq. '') then
      saxs_alpha = 200.0
    else
      read(dummy,*) saxs_alpha
    endif
    
    call getenv('SAXSPRINT',dummy)
    if (dummy .eq. '') then
        saxs_print = .false.
        print *, 'reading sim_f', saxs_print
    else
        read(dummy,*) saxs_print
    endif

  endif   
    


!=========================================================================================================

  if (taskid .eq. 0 ) then
     open(unit=FLOG,file=MASTER_LOGFILE,status='unknown',action='write',position='rewind',iostat=ierror)
     write(FLOG,*) '**************************************************************************'
     write(FLOG,'(A15,F8.3)') ' Version     : ', version
     write(FLOG,*) '**************************************************************************'
     write(FLOG,*) ' '
     write(FLOG,'(A39,I12)')   ' Number of atoms                     : ', NATOMS
     write(FLOG,'(A39,F12.6)') ' Temperature      (kelvin)           : ', target_temperature/0.59227*298
     write(FLOG,'(A39,F12.6)') ' Temperature      (kcal/mol)         : ', target_temperature
     write(FLOG,'(A39,I12)')   ' Number of fragments                 : ', NFRAG
     write(FLOG,'(A39,A12)')   ' Restriction on motion        : ', restrict_motion
     if (constrained_fragments.or.restrained_fragments) then
        write(FLOG,'(A39,F12.6)') '     Value of the spring constant    : ', k_spring
     endif
     write(FLOG,'(A39,I12)')   ' Random seed                         : ', idum
     write(FLOG,'(A39,I12)')   ' Number of digits in file name       : ', ndigits_filename
     write(FLOG,'(A39,A12)')   ' Debug status                        : ', debug_status
     write(FLOG,'(A39,F12.6)') ' Potential scaling factor            : ', force_scaling_factor

     write(FLOG,'(A39,L12)')   ' Periodic Boundary Condition         : ', PBC
     write(FLOG,'(A39,F12.6)') ' Box Length                          : ', BL
     write(FLOG,'(A39,L12)')   ' center of mass for pdb              : ', C_M
     write(FLOG,'(A39,A12)')   ' Minimization Type                   : ', MINIMIZATION_TYPE 
     write(FLOG,'(A39,A10)')   ' Simulation type                     : ', SIMULATION_TYPE

     if (SIMULATION_TYPE .eq. 'replica') then
        write(FLOG,'(A39,L12)')   ' Save PDB as single file         : ', singlefile
        write(FLOG,'(A39,I12)')   ' Number of replicas              : ', N_REPLICA
        write(FLOG,"         ('     Temperatures                    :   ',6F8.4)") T_replica
        write(FLOG,'(A39,I12)')   ' Exchange every n steps          : ', n_step_exchange
        write(FLOG,'(A39,L12)')   ' single file as initial for REMD : ', init_single_file
        write(FLOG,'(A39,A12)')   ' Replica type                    : ', REPLICA_TYPE
     endif
     if (SIMULATION_TYPE .eq. 'replica'  .and. REPLICA_TYPE .eq. 'E_Scale') then
        write(FLOG,'(A39,I12)')   ' Number of energy scales         : ', N_E_REPLICA
        write(FLOG,"         ('     energy scales                   :   ',6F8.4)") E_replica
     endif
     if (SIMULATION_TYPE .eq. 'saxs') then
        write(FLOG,'(A39,I12)')   ' SAXS Profiling Rate             : ', n_rate_saxs
        write(FLOG,'(A39,F12.6)') ' SAXS Energy Coefficient         : ', saxs_alpha
        write(FLOG,'(A39,I12)')   ' SAXS Norm Type                  : ', SAXS_norm_type
        write(FLOG,'(A39,F5.3)')  ' SAXS Vector Max                 : ', SAXS_Vect_Max
        write(FLOG,'(A39,I5)')    ' SAXS Vector Max Point           : ', max_q_point
        write(FLOG,'(A39,A12)')   ' SAXS Calculation Type           : ', merge(' In Solution', '    In Vacuo', in_solution_curve)
        write(FLOG,'(A39,L12)')   ' SAXS Refine hydration layer     : ', refine_hydration_layer
        write(FLOG,'(A39,F12.6)') ' SAXS Dummy Water Radius         : ', DX
        write(FLOG,'(A39,I12)')   ' SAXS Dummy Water Layers         : ', n_shells
     endif
     close(FLOG)
     open(FTIT, file='titration_MD.out', status='unknown',action='write',position='append')
     open(FTITpc, file='titration_PC.out', status='unknown',action='write',position='append')
     open(HBFILE, file='hb-on-the-fly.out', status='unknown',action='write',position='append')
     open(SAXSs, file='SAXS_score.dat', status='unknown',action='write',position='append')
     open(SAXSc, file='SAXS_curve.dat', status='unknown',action='write',position='append')
  endif

  call GET_ENVIRONMENT_VARIABLE('Simulation_Method', SIMULATION_METHOD)

  if (SIMULATION_METHOD .eq. "MD") then
     call read_parameters_md()
  else 
     write(*,*) 'Error : only MD accepted as simulation method'
     stop
  endif

  if(taskid .eq. 0 ) then
     open(unit=FLOG,file=MASTER_LOGFILE,status='unknown',action='write',position='append',iostat=ierror)
     write(FLOG,*)
     write(FLOG,'(A39,A12)')      ' Simulation method                   : ', SIMULATION_METHOD
     write(FLOG,*) ' '
     write(FLOG,*) '**************************************************************************'
     close(FLOG)
  endif
  allocate(pos(VECSIZE))       
  allocate(posref(VECSIZE))       
  allocate(force(VECSIZE))       
  allocate(atomic_type(NATOMS))
  allocate(mass(vecsize))

  ! We first set-up pointers for the x, y, z components in the position and
  ! forces

  x    => pos(1:3*natoms:3)
  y    => pos(2:3*natoms:3)
  z    => pos(3:3*natoms:3)

  xref => posref(1:3*natoms:3)
  yref => posref(2:3*natoms:3)
  zref => posref(3:3*natoms:3)

  fx   => force(1:3*natoms:3)
  fy   => force(2:3*natoms:3)
  fz   => force(3:3*natoms:3)

  return
end subroutine definitions

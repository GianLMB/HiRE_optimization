module SAXS_scoring

  use geometric_corrections
  implicit none

  ! A (somewhat) protected UNIT for file I/O
  integer :: QFILE = 123
  character*1, parameter :: comment_pattern='#'
  ! FA / CG parameters
  character*30 :: parameter_file = 'SAXS_grains.dat' 
  character*180 :: saxs_input_file = "SAXS"
  integer :: num_atoms ! A raw copy of NATOMS (defs)
  integer :: grain_number = 0
  character*5, dimension(:), allocatable :: atom_names
  character*5, dimension(:), allocatable :: Grains
  real*8, dimension(:,:), allocatable :: F_grains
  real*8, dimension(:), allocatable :: Grains_cutoff
  ! SAXS Target Curve Parameters
  logical :: use_target = .true.
  real*8, dimension(:), allocatable :: target_curve
  real*8 :: saxs_match_threshold  = 0.5 ! Threshold for check_SAXS_consistency
  ! SAXS Curve Parameters
  real*8 :: max_q = 1.0d0
  real*8 :: SAXS_vect_max = 1.0d0 
  integer :: max_q_point
  integer :: num_points = 200
  real*8 :: delta_q
  real*8 :: delta_w
  ! SAXS conformation arrays and parameters
  character*5, dimension(:), allocatable :: AtomNames_TOT
  real*8, dimension(:), allocatable :: pos_SOL, pos_TOT
  real*8, dimension(:,:), allocatable :: DistanceMatrix_TOT, F_CG_TOT
  integer :: num_SOL, num_TOT
  ! Computation Schemes
  real*8 :: SAXS_alpha = 200.0 ! SAXS coupling constant for Replicas
  integer :: SAXS_norm_type = 2
  logical :: use_mean_correction = .true.
  real*8, dimension(:), allocatable :: mean_correction_curve
  character*100 :: mean_correction_file
  logical :: in_solution_curve = .false. 
  logical :: explicit_sol_contribution = .false.
  logical :: refine_hydration_layer = .false. 
  logical :: add_hydration_shell = .false. 
  logical :: linear_intensity = .false.
  ! Solvent Parameters
  integer :: MAXSOL = 5000 ! Slightly dynamic, increased if necessary
  character*4 :: solvent_name = 'HOH '
  real*8 :: solvent_electron_density = 0.334
  real*8 :: SAXS_w_shell = 0.3d0
  real*8 :: dx
  integer :: n_shells

  contains

    ! lm759> Initialisation for SAXS
    subroutine initialise_SAXS()

      use defs
      implicit none
      integer :: point_q
      logical :: found_mean_correction = .false.
      real*8 :: q

      ! Basic definitions
      if(.not.allocated(atom_names)) allocate(atom_names(NATOMS))
      atom_names(:) = ATOMIC_TYPE(:)
      num_atoms = NATOMS
      delta_q = max_q / num_points
      max_q_point = int(SAXS_Vect_Max / delta_q)
      add_hydration_shell = in_solution_curve .and. refine_hydration_layer

      ! Allocate and fill SAXS target
      if (use_target .and. .not. allocated(target_curve)) then
        allocate(target_curve(0: 100 - 1))   !lal
        open(unit=QFILE, file='saxs_target.dat', status='old')
        q_loop: do point_q = 0, 100 - 1     !lal
          read(QFILE,*) q, target_curve(point_q)
        enddo q_loop
        close(QFILE)
      endif

      ! Allocate and fill SAXS mean correction
      if (in_solution_curve) then
        mean_correction_file = "saxs_mean_correction_Solution.dat"
      else
        mean_correction_file = "saxs_mean_correction_Vacuo.dat"
      endif
      inquire(file=mean_correction_file, exist=found_mean_correction)
      if (use_mean_correction .and. .not. found_mean_correction) then
        print *, 'Could not find mean correction file ', trim(mean_correction_file)
        STOP 5
      end if
      use_mean_correction = use_mean_correction .and. found_mean_correction
      if (use_mean_correction .and. .not. allocated(mean_correction_curve)) then
        allocate(mean_correction_curve(0:num_points-1))
        open(unit=QFILE, file=mean_correction_file, status='old')
        do point_q = 0, num_points-1
          read(QFILE,*) q, mean_correction_curve(point_q)
        enddo
        close(QFILE)
      endif
  
      ! Global initialization from parameter file
      call init_hash_CG()

    end subroutine initialise_SAXS

    ! lm759> Finalization for SAXS
    subroutine finalise_SAXS()
      if(use_target) deallocate(target_curve)
      if (in_solution_curve) deallocate(mean_correction_curve)
      deallocate(atom_names)
      deallocate(Grains, F_Grains, Grains_cutoff)
    end subroutine finalise_SAXS

    ! lm759> Initialization for the locally used conformation
    subroutine initialise_SAXS_conformation(pos)

      use defs, only: use_qbug
      implicit none
      real*8, dimension(1:3*num_atoms) :: pos
      real*8 :: start, finish

      ! Fill the shell and compute num_TOT
      if (use_qbug) call cpu_time(start)
      if_add_hydration_shell: if (add_hydration_shell) then
        allocate(pos_SOL(1:3*MAXSOL))
        pos_SOL = 0.0
        call hydrate_pos_outer_shell(pos)
        num_TOT = num_atoms + num_SOL
      else
        num_TOT=num_atoms
      endif if_add_hydration_shell

      ! Allocate stuff
      allocate(pos_TOT(1:3*num_TOT))
      allocate(AtomNames_TOT(num_TOT))
      allocate(F_CG_TOT(0:num_points-1,1:num_TOT))
      allocate(DistanceMatrix_TOT(1:num_TOT,1:num_TOT))

      ! Fill conformation arrays
      pos_TOT(1:3*num_atoms) = pos(1:3*num_atoms)
      AtomNames_TOT(1:num_atoms) = atom_names(1:num_atoms)
      if (add_hydration_shell) then
        pos_TOT(3*num_atoms+1:3*num_TOT) = pos_SOL(1:3*num_SOL)
        AtomNames_TOT(num_atoms+1:num_TOT) = 'HOH  '
      endif
      call fill_structureFactor_array(num_TOT, AtomNames_TOT, F_CG_TOT)
      call fill_half_distance_matrix(num_TOT, pos_TOT, DistanceMatrix_TOT)

      ! Print parameters
      if (use_qbug) then
        call cpu_time(finish)
        PRINT *, " "
        PRINT *, " "
        PRINT '("CONFORMATION TIME: ", f6.3)', FINISH-START
        PRINT *, "  Number of atoms: ", num_atoms, num_sol, num_tot
        PRINT *, " "
        PRINT *, " "
        close(222)
      endif

    end subroutine initialise_SAXS_conformation

    ! lm759> SAXS logaritmic intensity
    function SAXS_intensity(pos) result (logI)

      implicit none
      real*8, dimension(1:3*num_atoms) :: pos
      real*8, dimension(0:max_q_point-1) :: logI, I1, I0
      integer :: grain_i, grain_j, q

      ! Compute intensity
      I1 = 0.0d0
      do grain_i = 1, num_TOT
        I1(0) = I1(0) + F_CG_TOT(0,grain_i)**2
        q_loop: do q = 1, max_q_point-1
          I1(q)  = I1(q) + F_CG_TOT(q,grain_i)**2
        enddo q_loop
        do grain_j = grain_i+1, num_TOT
          I1(0) = I1(0) + 2 * F_CG_TOT(0,grain_i) * F_CG_TOT(0,grain_j)
        q_loop1: do q = 1, max_q_point-1
             I1(q)  = I1(q) + 2 * F_CG_TOT(q,grain_i) * F_CG_TOT(q,grain_j) * &
                 & sin(DistanceMatrix_TOT(grain_j,grain_i) * q * delta_q) / &
                 & (q * delta_q * DistanceMatrix_TOT(grain_j,grain_i)) 
        enddo q_loop1
        enddo
      enddo

      ! Log10 everything
        logI(:) = log10(I1(:)) 

      ! Add mean correction 
      if (use_mean_correction) then
        logI(:) = logI(:) + mean_correction_curve(:)
      endif

      ! UnLog10 everything lal
!        I1(:)=10**logI(:)

      ! Rescale lal
!      if (in_solution_curve .and. refine_hydration_layer) then
!        I1(:) = I1(:) * ( 10**target_curve(0) / I1(0) )
!      endif 

      ! ReLog10 everything lal
!        logI(:) = log10(I1(:)) 

      return
      
    end function SAXS_intensity

    ! lm759> SAXS cruve, energy and force implementation (as 2-norm)
    subroutine SAXS_energy_and_force(pos, logI, Esaxs, F_saxs)

      use defs, only: use_qbug
      implicit none
      real*8, dimension(1:3*num_atoms), intent(in) :: pos
      real*8, dimension(0:max_q_point-1), intent(out) :: logI
      real*8, dimension(1:3*num_atoms), intent(out) :: F_saxs
      real*8, intent(out) :: Esaxs
      real*8 :: r(3), r2, qr, qphys, cscale, start, finish
      real*8 :: E_den, E_num, F_pre, F_grain(3)
      real*8, dimension(0:max_q_point-1) :: I1, I0
      integer :: grain_i, grain_j, q

      ! Allocate and fill conformation arrays
      if (use_qbug) call cpu_time(start)
      call initialise_SAXS_conformation(pos) 
      F_saxs = 0.0d0
      E_den = 0.0d0
      E_num = 0.0d0

      ! Compute intensity
      logI(:) = SAXS_intensity(pos)
      I1(:) = 10**logI(:)
      I0(:) = 10**target_curve(:) 

      ! Compute E_saxs and F_saxs (q**2 weight)
!      cscale = 1
      cscale = I0(0)/I1(0)    !lal
      q_loopE: do q = 3, max_q_point - 1 !lal, changed 1 to 3, date 14/04/2021
        qphys = q * delta_q
        E_num = E_num + ((cscale * I1(q) - I0(q)) * qphys)**2
        F_pre = - 2 * cscale*(cscale * I1(q) - I0(q)) * qphys**2 / (max_q_point * I0(0)**2)
         grain_i_loop: do grain_i = 1, num_atoms
            grain_j_loop: do grain_j = grain_i+1, num_tot
               r = pos_TOT(grain_i*3-2:grain_i*3) - pos_TOT(grain_j*3-2:grain_j*3)
               r2 = DistanceMatrix_TOT(grain_j,grain_i)**2
               qr = qphys * DistanceMatrix_TOT(grain_j,grain_i)
               F_grain = ( r / r2 ) * F_CG_TOT(q, grain_i) * F_CG_TOT(q, grain_j) * ( cos(qr) - sin(qr) / qr )
               F_saxs(grain_i*3-2:grain_i*3) = F_saxs(grain_i*3-2:grain_i*3) + 2 * F_pre * F_grain
               if (grain_j .le. num_atoms) then
                   F_saxs(grain_j*3-2:grain_j*3) = F_saxs(grain_j*3-2:grain_j*3) - 2 * F_pre * F_grain
               endif
            enddo grain_j_loop
         enddo grain_i_loop
      enddo q_loopE
      Esaxs = E_num / (max_q_point * I0(0)**2) !lal, changed I1 to I0, date 08/04/2021
      
      print *, Esaxs, F_pre
      
      ! Print parameters
      if (use_qbug) then
        call cpu_time(finish)
        PRINT *, " "
        PRINT *, " "
        PRINT '("ENERGY AND FORCE TIME: ", f6.3)', FINISH-START
        PRINT *, "  Number of atoms: ", num_TOT
        PRINT *, "  Maximum value for scattering vector: ", saxs_vect_max 
        PRINT *, "  Energy: ", Esaxs
        PRINT *, " "
        PRINT *, " "
        close(222)
      endif

      ! Deallocate stuff
      deallocate(F_CG_TOT, AtomNames_TOT, DistanceMatrix_TOT, pos_TOT)
      if (add_hydration_shell) deallocate(pos_SOL)

      return
      
    end subroutine SAXS_energy_and_force

    ! lm759> SAXS energy wrt target and save curve
    function SAXS_energy(pos, write_to_unit) result(energy)

      implicit none
      real*8, dimension(1:3*num_atoms) :: pos
      real*8, dimension(0:max_q_point-1) :: logI1, I1, I0
      real*8 :: energy
      integer, optional :: write_to_unit

      ! Allocate and fill conformation arrays
      call initialise_SAXS_conformation(pos)

      ! Compute intensity
      logI1(:) = SAXS_intensity(pos)
      I1(:) = 10**logI1(:) 
      I0(:) = 10**target_curve(:)

      ! Compute approprate energy 
      select case (SAXS_norm_type)
        case (1)
          energy = norm1_curves(I1, I0)
        case (2)
          energy = norm2_curves(I1, I0)
        case (3)
          energy = norm_R_curves(I1, I0)
        case default
          STOP 'Undefined norm type. Aborting.'
          energy = 0.0d0
      end select

      ! Save SAXS curve
      if (present(write_to_unit)) call write_SAXS_curve_to_unit(logI1, write_to_unit)

      ! Deallocate stuff
      deallocate(F_CG_TOT, AtomNames_TOT, DistanceMatrix_TOT, pos_TOT)
      if (add_hydration_shell) deallocate(pos_SOL)

    end function SAXS_energy

    !> @brief Hydration shell computation
    subroutine hydrate_pos_outer_shell(pos)

      use defs, only: use_qbug
      implicit none
      real*8, intent(in), dimension(1:3*num_atoms) :: pos
      real*8, dimension(1:num_atoms) :: CutoffArrayInternal_pow2, CutoffArrayExternal_pow2
      integer :: lattice_x, lattice_y, lattice_z, atom_i, contact_num
      integer :: pdb_res_num, pdb_atom_num
      real*8 :: r_ij, x, y, z, cut_off_internal_2, cut_off_external_2, start, finish
      real*8, dimension(1:6) :: MinMaxCoordsArray
      integer, dimension(1:6) :: latticeArray

      ! Open file to save hydration shell
      if (use_qbug) then 
        call cpu_time(start)
        open(unit=222,file=trim(saxs_input_file) // '_hyd.pdb', form='formatted')
        pdb_res_num = 1
        pdb_atom_num = 1
      endif

      ! Initialise and fill stuff
      num_SOL = 0
      call fill_cutoff_array(atom_names, CutoffArrayInternal_pow2, CutoffArrayExternal_pow2)
      call box_dimensions(pos, MinMaxCoordsArray, sqrt(maxval(CutoffArrayExternal_pow2, num_atoms))) 
      latticeArray = int(MinMaxCoordsArray/dx)

      ! Compute coordinates for solvent grains
      x_loop: do lattice_x=latticeArray(1), latticeArray(2)
        x = lattice_x*dx
        y_loop: do lattice_y=latticeArray(3), latticeArray(4)
          y = lattice_y*dx
          z_loop: do lattice_z=latticeArray(5), latticeArray(6)
            z = lattice_z*dx
            contact_num = 0
            atom_loop: do atom_i=1,num_atoms
              cut_off_internal_2 = CutoffArrayInternal_pow2(atom_i)
              cut_off_external_2 = CutoffArrayExternal_pow2(atom_i)
              r_ij = (pos(3*atom_i-2)-x)**2
              if ( r_ij <= cut_off_external_2 ) then
                r_ij = r_ij + (pos(3*atom_i-1)-y)**2
                if ( r_ij <= cut_off_external_2 ) then
                  r_ij= r_ij + (pos(3*atom_i)-z)**2
                  clash: if(r_ij <= cut_off_internal_2 ) then
                    contact_num = 0
                    exit
                  endif clash
                  contact: if (r_ij <= cut_off_external_2 ) then
                    contact_num = contact_num + 1
                  endif contact
                endif
              endif
            enddo atom_loop
            add: if (contact_num /= 0 ) then
              num_SOL = num_SOL + 1
              pos_SOL(3*num_SOL -2) = x
              pos_SOL(3*num_SOL -1) = y
              pos_SOL(3*num_SOL   ) = z
              ! Save coordinate to hydration shell file
              if (use_qbug) then
                write (222,'(a6,i5,a1,a4,a1,a3,a2,i4,a4,3f8.3)') &
                   'ATOM  ',pdb_atom_num,' ', 'HOH', ' ', 'SOL', ' A', pdb_res_num,'    ', x,y,z
                pdb_atom_num = merge(pdb_atom_num + 1,1,pdb_atom_num .lt. 99999)
                pdb_res_num = merge(pdb_res_num + 1,1,pdb_res_num .lt. 9999)
              endif
            endif add
          enddo z_loop
        enddo y_loop
      enddo x_loop

      ! Print parameters
      if (use_qbug) then
        call cpu_time(finish)
        PRINT *, " "
        PRINT *, " "
        PRINT '("SHELL TIME: ", f6.3)', FINISH-START
        PRINT *, "  Water radius: ", dx / 2.0
        PRINT *, "  Water contrast: ", SAXS_w_shell
        PRINT *, "  Number of shell atoms: ", num_SOL
        PRINT *, " "
        PRINT *, " "
        close(222)
      endif

    end subroutine hydrate_pos_outer_shell

    ! lm759> Compares normalization of target and computed intensity
    subroutine check_SAXS_consistency(pos)

      use hash_int
      implicit none
      real*8, dimension(1:3*num_atoms), intent(in) :: pos
      real*8, dimension(:), allocatable :: F_CG_TOT_0
      real*8 :: I1, I0 
      integer :: grain_i, grain_j, grain_hi, q

      ! Allocate stuff: either hydrated or not, but only for q=0
      if_add_hydration_shell: if (add_hydration_shell) then
        allocate(pos_SOL(1:3*MAXSOL))
        pos_SOL = 0.0
        call hydrate_pos_outer_shell(pos)
        num_TOT = num_atoms + num_SOL
      else
        num_TOT = num_atoms
      endif if_add_hydration_shell
      allocate(AtomNames_TOT(num_TOT))
      allocate(F_CG_TOT_0(num_TOT))
      AtomNames_TOT(1:num_atoms) = atom_names(1:num_atoms)
      if (add_hydration_shell) then
        AtomNames_TOT(num_atoms+1:num_TOT) = 'HOH  '
      endif

      ! Fill form factors at q=0 only
      FF_at_0: do grain_i=1,num_TOT
        call hash_get(adjustl(trim(AtomNames_TOT(grain_i))), grain_hi)
        if (grain_hi /= 0 ) then
          F_CG_TOT_0(grain_i) = F_grains(0, grain_hi)
        else
          print *, "Couldn't find structure factor parameters for Grain : "
          print '(i5 a)', grain_i, atom_names(grain_i)
          F_CG_TOT_0(grain_i) = 0.0
          STOP 5
        endif
      enddo FF_at_0

      ! Compute intensity at q=0
      I1 = 0.0d0
      do grain_i = 1, num_TOT
        I1 = I1 + F_CG_TOT_0(grain_i)**2
        do grain_j = grain_i+1, num_TOT
          I1 = I1 + 2 * F_CG_TOT_0(grain_i) * F_CG_TOT_0(grain_j)
        enddo
      enddo
      I0 = 10**target_curve(0)
      if (use_mean_correction) I1 = 10**(log10(I1) + mean_correction_curve(0))

      if (abs(I1 - I0) >= saxs_match_threshold * I0) then
        print '(A ES15.7 A ES15.7 A ES10.1 A)', &
           &'check_SAXS_consistency(): The target SAXS profile (',I0,') &
           &does not match the calculated SAXS profile (',I1,') in q=0, &
           &with threshold ',saxs_match_threshold,', &
           &which is symptomatic of something being very wrong. ABORTING'
        STOP 5
      endif
 
      ! Deallocate stuff
      deallocate(AtomNames_TOT, F_CG_TOT_0)
      if (add_hydration_shell) deallocate(pos_SOL)

    end subroutine check_SAXS_consistency

    ! lm759> SAXS energy with q**2 weigthing factor
    function norm2_curves(curve1, curve2) result(norm)
      implicit none
      integer :: q
      real*8, dimension(0 : max_q_point - 1) :: curve1, curve2
      real*8 :: norm
      norm = 0.0d0
      do q = 1, max_q_point-1
        norm = norm + ((curve1(q) - curve2(q)) * (q * delta_q))**2
      enddo
      norm = norm / (max_q_point * (curve1(0))**2)
      return
    end function norm2_curves

    ! lm759> Linear SAXS energy
    function norm1_curves(curve1, curve2) result(norm)
      implicit none
      integer :: q
      real*8, dimension(0:max_q_point-1) :: curve1, curve2
      real*8 :: norm
      norm = 0.0d0
      do q = 0, max_q_point - 1
        norm = norm + dabs(curve1(q) - curve2(q))
      enddo
      norm = norm / max_q_point
      return
    end function norm1_curves

    ! lm759> R-factor SAXS energy
    function norm_R_curves(curve1, curve2) result(norm)
      implicit none
      integer :: q
      real*8, dimension(0:max_q_point-1) :: curve1, curve2
      real*8 :: norm, norm_num, norm_den, sc, sc_num, sc_den
      norm_num = 0.0d0
      norm_den = 0.0d0
      sc_num = 0.0d0
      sc_den = 0.0d0
      !do q = 1, max_q_point - 1
      !  sc_num = sc_num + ( q * delta_q )**2 * 10**curve1(q) * 10**curve2(q) 
      !  sc_den = sc_den + ( q * delta_q )**2 * ( 10**curve1(q) )**2 
      !end do
      !sc = sc_num / sc_den
      sc = 1 ! If scale factor needs to be used, uncomment previous lines, comment this
      do q = 1, max_q_point - 1
        norm_num = norm_num + (q * delta_q * (sc * 10**curve1(q) - 10**curve2(q)))**2
        norm_den = norm_den + (q * delta_q)**2 * sc * 10**curve1(q) * 10**curve2(q) 
      enddo
      norm = norm_num / norm_den 
      return
    end function norm_R_curves

    ! lm759> Save intensity and energy data 
    subroutine write_SAXS_curve_to_unit(curve1, unit_number)
      implicit none
      real*8, dimension(0:max_q_point-1), intent(in) :: curve1
      integer, intent(in) :: unit_number
      integer :: q_point
      do q_point = 0, max_q_point-1
        if (isNaN(curve1(q_point))) then
          print *, curve1
          print *, 'q_point: ', q_point
          STOP 'NaN Error'
        endif
       if (linear_intensity) then
         write (unit_number,'(f8.3 ES15.7)') q_point*delta_q, 10**curve1(q_point)
       else
         write (unit_number,'(f8.3 ES15.7)') q_point*delta_q, curve1(q_point)
       endif
      enddo
      write (unit_number,*) ' '
    end subroutine write_SAXS_curve_to_unit

    !> @brief Initialize grains with form factors
    subroutine init_hash_CG()

      use defs, only: use_qbug
      use hash_int
      implicit none
      integer :: grain_i, q
      real*8 :: bufferValue, bufferValue2, bufferValue3, bufferValue4, temp
      character*180 :: datFile, parameter_line
      character*5 :: bufferKey
      integer :: NF=18
      integer :: grain_hash_num, stat

      ! Get grain number from parameter_file
      grain_number = 0
      open(unit=NF, file=parameter_file, status="old")
      get_grain_number: do
        read(NF,'(a)', iostat=stat) parameter_line
        if (stat /= 0) exit
        if (parameter_line(1:1) /= comment_pattern) grain_number = grain_number + 1
      enddo get_grain_number
      close(NF)

     ! Allocate stuff
      if (.not. allocated(Grains)) allocate(Grains(1:grain_number))
      if (.not. allocated(F_Grains)) allocate(F_Grains(0:num_points-1,1:grain_number))
      if (.not. allocated(Grains_cutoff)) allocate(Grains_cutoff(1:grain_number))
     
      ! Make the hash with the grain names
      call hash_init()
      open(unit=NF, file=parameter_file, status="old")
      grain_i = 1
      build_grain_hash: do while (grain_i <= grain_number)
        read(nf,'(a)', iostat=stat) parameter_line
        if (parameter_line(1:1) /= comment_pattern) then
          read(parameter_line,'(a5)') bufferKey
          call hash_set(bufferKey, grain_i)
          Grains(grain_i)=bufferKey
          grain_i = grain_i + 1
        end if
      enddo build_grain_hash
      close(NF)

      ! Fill the array Grains_cutoff
      open(unit=NF, file=parameter_file, status="old")
      grain_i = 1
      build_grain_arrays: do while (grain_i <= grain_number)
        read (nf, '(a)') parameter_line
        if (parameter_line(1:1) /= comment_pattern) then
          read(parameter_line,'(a5 f4.1 f4.1 f6.1 f6.1 x a1)') &
            bufferKey, bufferValue, bufferValue2, bufferValue3, bufferValue4
          call hash_get(bufferKey, grain_hash_num)
          Grains_cutoff(grain_hash_num) = bufferValue2
          grain_i = grain_i + 1
        endif
      enddo build_grain_arrays
      close(NF)

      ! Build F_grains from the dat/ff* files
      read_dat_files: do grain_i=1, grain_number
        if ( in_solution_curve .and. .not. explicit_sol_contribution ) then
          datFile= 'dat/ff' // trim(Grains(grain_i)) // '.cor.awk.dat' ! Solvent-corrected FF
        else
          datFile= 'dat/ff' // trim(Grains(grain_i)) // '.awk.dat' ! In vacuo FF
        endif
        open(unit=NF, file=datFile, status='old')
        copy_q: do q=0, num_points - 1
          read(NF,'(F8.3 F10.8)') temp, F_grains(q, grain_i)
        enddo copy_q
        close(NF)
        if ( trim(Grains(grain_i)) .eq. 'HOH' ) then
          F_grains(:, grain_i) = SAXS_w_shell * F_grains(:, grain_i) 
          if (use_qbug) PRINT *, "Check solvent excess electron density: ", F_grains(0,grain_i)
        endif
      enddo read_dat_files

    end subroutine init_hash_CG

    !> @brief Fill coarse-grain form factors arrays
    subroutine fill_structureFactor_array(num_grain, GrainName, F_q_CG)
      use hash_int
      implicit none
      integer, intent(in) :: num_grain
      character*5, intent(in), dimension(1:num_grain) :: GrainName
      real*8, intent(out),dimension(0:num_points-1, 1:num_grain ) :: F_q_CG
      integer :: grain_hash_num
      integer :: grain_i
      match_cutoff: do grain_i=1,num_grain
        call hash_get(adjustl(trim(GrainName(grain_i))), grain_hash_num)
        if (grain_hash_num /= 0 ) then
          F_q_CG(0:num_points-1, grain_i) = F_grains(0:num_points-1, grain_hash_num)
        else
          print *, "Couldn't find structure factor parameters for Grain : "
          print '(i5 a)', grain_i, GrainName(grain_i)
          F_q_CG(0:num_points-1, grain_i) = 0.0
          STOP 5
        endif
      enddo match_cutoff
    end subroutine fill_structureFactor_array

    ! Fill coarse-grain cutoff array for hydration shell
    subroutine fill_cutoff_array(AtomNames, CutoffArrayInternal_pow2, CutoffArrayExternal_pow2)
      use hash_int
      implicit none
      character*5, intent(in), dimension(1:num_atoms) :: AtomNames
      real*8, dimension(1:num_atoms), intent(out) :: CutoffArrayInternal_pow2
      real*8, dimension(1:num_atoms), intent(out) :: CutoffArrayExternal_pow2
      integer :: atom_i, grain_hash_num
      character*5 :: atom_name
      match_cutoff: do atom_i=1,num_atoms
          atom_name = trim(AtomNames(atom_i))
          call hash_get(atom_name, grain_hash_num)
          if (grain_hash_num .eq. 0) then
            print *, "Couldn't find parameters for atom/grain : ",atom_name
            STOP 6
          endif
          CutoffArrayInternal_pow2(atom_i) = ( 0.5   * dx                    + Grains_cutoff(grain_hash_num) )**2
          CutoffArrayExternal_pow2(atom_i) = ( 0.99  * (n_shells + 0.5) * dx + Grains_cutoff(grain_hash_num) )**2
      enddo match_cutoff
    end subroutine fill_cutoff_array

    ! Set the limits for the water box
    subroutine box_dimensions(XCoords, MinMaxCoordsArray, padding_length)
      use defs, only: use_qbug
      implicit none
      integer :: i
      real*8,intent(in) :: padding_length ! In nm
      real*8,dimension(1:6),intent(out) ::  MinMaxCoordsArray
      real*8  XCoords(3*num_atoms)
      MinMaxCoordsArray=(/ Minval(XCoords(::3)), Maxval(XCoords(::3)), &
      & Minval(XCoords(2::3)), Maxval(XCoords(2::3)), &
      & Minval(XCoords(3::3)), Maxval(XCoords(3::3)) /)
      do i=1,6
        MinMaxCoordsArray(i) = merge(MinMaxCoordsArray(i)-padding_length, &
        & MinMaxCoordsArray(i)+padding_length, modulo(i,2) == 1)
      enddo
      return
    end subroutine box_dimensions

    ! Compute pair-wise distances
    subroutine fill_half_distance_matrix(num_grain, pos, DistanceMatrix)
      implicit none
      integer, intent(in) :: num_grain
      integer :: i,j
      real*8, dimension(1:num_grain,1:num_grain), intent(out) :: DistanceMatrix
      real*8, dimension(1:3*num_grain) ::  pos
      real*8 :: r_ij
      DistanceMatrix=0.d0
      i_loop: do i=1,num_grain
        j_loop: do j=i+1,num_grain
          r_ij = sqrt((pos(3*i-2)-pos(3*j-2))**2+(pos(3*i-1)-pos(3*j-1))**2+(pos(3*i)-pos(3*j))**2)
          DistanceMatrix(i,j)=r_ij
          DistanceMatrix(j,i)=r_ij
        enddo j_loop
      enddo i_loop
    end subroutine fill_half_distance_matrix

end module SAXS_scoring

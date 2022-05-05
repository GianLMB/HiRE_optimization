! This set of routines serves for the initialisation 
!
! Copyright Normand Mousseau January 2006

module md_defs
  use defs
  use geometric_corrections

  implicit none
  save

  real(8) :: timeunit, dt, dt1, dt2, timestep
  real(8), dimension(:), allocatable, target :: vel
  real(8), dimension(:), allocatable :: invmass
  real(8), dimension(:), allocatable :: dt1_mass, dt2_mass
  real(8) :: friction_coef
  real(8) :: langevin_scalev
  real(8), dimension(:), allocatable :: temperatures
  real(8), dimension(:), pointer :: vx, vy, vz
  real(8) :: initial_simulation_time
  integer(8) :: restart_simulation_time

  logical :: confining_sphere
  real(8) :: radius_confining_sphere,radius_confining_sphere2, k_confine

  integer(8) :: n_production
  integer :: n_correct_com, n_correct_rotation
  integer :: n_stats, n_stats_average, berendsen_time
  !Br1
  integer :: n_titration
  integer :: n_save_configs, n_equilibration, n_rescale_v_equil
  integer :: n_steps_thermalization, n_save_restarts
  real(8) :: simulation_time
  character(len=20) :: thermostat
  character(len=20) :: LOGFILE

  real(8) :: avstat, avpot, avkin, avetot
  real(8) :: avtemp, avpot2, avkin2, avetot2, avtemp2
  
  logical :: save_unrotated_structures
  
  ! Gaussian number
  integer :: gaussian_flag
  real(8) :: gaussian_number2

  ! for rattle
  logical :: rattle
  integer :: nrattle
  integer :: nbondh, nbonda, nbondt, nbonds

  ! For contrained motion of hydrogen
  logical :: control_hydrogen

  double precision, allocatable, dimension(:) :: beq, beq2, ibeq2, rij2
  double precision, allocatable, dimension(:,:) :: rij1
  double precision, allocatable, dimension(:) :: redu1_mass, redu2_mass, reduA, reduB, reduAA, reduBB
  integer, allocatable, dimension(:) :: ia, ib
  double precision :: epspos, epsvel
 
  logical :: thermo
  logical :: lj_test
  logical :: force_calc_test
  logical :: chk_ene

  ! IMD parameters
  logical :: is_imd
  integer :: imd_wait
  integer :: imd_port
  integer :: imd_debug
  real(8) :: imd_forcescale

  ! saxs serial parameters
  logical :: compute_saxs_serial
  integer :: saxs_serial_step, saxs_serial_step_min 
  ! Boolean deciding if we should compute a SAXS forces next step
  logical :: calc_SAXS_force = .false.
  ! Parameters for SAXS modulation profile
  logical :: modulate_SAXS_serial = .false. 
  logical :: calc_SAXS_modul = .false.
  double precision :: SAXS_wave, SAXS_invsig, SAXS_modstep, SAXS_onoff


  contains

  ! Routine that interfaces with OPEP
  subroutine calcforce(scale,xa,fa,etot)
    real(8), intent(in) :: scale
    real(8), intent(out) :: etot
    real(8), dimension(vecsize), intent(inout) :: xa
    real(8), dimension(vecsize), intent(out) :: fa

    real(8) :: E_prot, E_RNA, E_cross

    ! Call the force field
    if (prot_simulation) then
      call calcforce_protein(scale,xa(1:N_prot*3),fa(1:N_prot*3),E_prot)
    endif
    if (RNA_simulation) then
      call calcforce_RNA(scale,xa(N_prot*3+1:Natoms*3),fa(N_prot*3+1:Natoms*3),E_RNA,.false., 0.0d0)
    endif
    if (prot_simulation .and. RNA_simulation) then
      call calcforce_RNA_protein(scale, xa, fa, E_cross)
    endif
    etot = E_prot+E_RNA+E_cross
 
  end subroutine calcforce
  

  !> @brief  Routine that interfaces with OPEP
  !> @brief  this one is modified to allow for logging of the energy terms
  subroutine calcforce_log(scale,xa,fa,etot,logener, simuPercent)
    real(8), intent(in) :: scale
    real(8), intent(out) :: etot
    real(8), dimension(vecsize), intent(inout) :: xa
    real(8), dimension(vecsize), intent(out) :: fa

    double precision simuPercent
    logical, intent(in) :: logener
    real(8) :: E_prot, E_RNA, E_cross

    ! Call the force field
    if (prot_simulation) then
      call calcforce_protein(scale,xa(1:N_prot*3),fa(1:N_prot*3),E_prot)
    endif
    if (RNA_simulation) then
      call calcforce_RNA(scale,xa(N_prot*3+1:Natoms*3),fa(N_prot*3+1:Natoms*3),E_RNA,logener,simuPercent)
    endif
    if (prot_simulation .and. RNA_simulation) then
      call calcforce_RNA_protein(scale, xa, fa, E_cross)
    endif
    etot = E_prot+E_RNA+E_cross
 
  end subroutine calcforce_log

end module md_defs

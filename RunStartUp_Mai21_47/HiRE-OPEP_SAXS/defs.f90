module defs
  ! This module defines all variables used accross the program ART01
  !
  !  Copyright N. Mousseau, May 2001
  implicit none

#ifdef MPI
  include 'mpif.h'
#endif
  save
  ! ! Numero de version
  real(8), parameter :: version = 2.0 
  ! Lists of parameters
  real(8) :: target_temperature
  integer :: NFRAG 
  integer :: NATOMS 
  integer :: VECSIZE
  integer :: VECSIZE1                     ! Length of the force and position vectors of one fragment 

  double precision :: degree_freedom               ! degree of freedom, for nonlinear molecules, degree_freedom = n*3 - 6 - num_constraints

  integer(8) :: NUMBER_EVENTS                ! Total number of events in this run
  character(len=3)  :: SIMULATION_METHOD  ! Type of simulation ART or MD 
  character(len=10)  :: SIMULATION_TYPE    ! serial or replica (exchange) or tempering
  character(len=10) :: REPLICA_TYPE       ! T_Exchange or E_Scale
  character(len=10) :: MINIMIZATION_TYPE       ! MIN_damped or MIN_steep

  logical :: constrained_fragments, restrained_fragments
  real(8) :: k_spring

  character(len=5), dimension(:), allocatable :: atomic_type ! Atomic type
  real(8), dimension(:), allocatable, target :: force       ! Working forces on the atoms
  real(8), dimension(:), allocatable, target :: pos         ! Working positions of the atoms
  real(8), dimension(:), allocatable, target :: posref      ! Reference position
  real(8), dimension(:), allocatable, target :: mass        ! masses
  real(8) :: force_scaling_factor  ! Factor for rescaling the potential

  real(8), dimension(:), allocatable, target :: imd_forces ! Forces from IMD

  real(8), dimension(:), pointer :: x, y, z
  real(8), dimension(:), pointer :: xref, yref, zref
  real(8), dimension(:), pointer :: fx, fy, fz

  integer, dimension(:,:), allocatable :: list_fragments
  
  integer :: mincounter                     ! Counter for output files
  integer :: refcounter                     ! Id of reference file
  integer(8) :: niter
  integer :: nitt
  integer :: evalf_number  = 0              ! Number of force evalutions
  integer :: ndigits_filename

  real(8) :: total_energy, ref_energy       ! Energies
  character(len=20) :: restart

  logical :: usextc
  logical :: singlefile 
 
  integer :: T_id
  integer :: E_scale
  logical :: init_single_file
  logical :: PBC  ! periodic boundary condition
  logical :: C_M  ! center of mass for writing non-broken chains in pdb file
  real(8) :: BL   ! box length 

  logical :: prot_simulation, RNA_simulation
  integer :: N_prot, N_RNA

  integer :: flag_tit
  integer :: n_steps_tit

  ! Units for writing/reading
  integer, parameter :: FCONF        = 1         
  integer, parameter :: FCOUNTER     = 2        
  integer, parameter :: FLIST        = 3         
  integer, parameter :: FLOG         = 4         
  integer, parameter :: FREFCONFIG   = 11
  integer, parameter :: FENER        = 12
  integer, parameter :: FRESTART     = 13
  integer, parameter :: FREP         = 17
  integer, parameter :: FSAXS        = 19
  integer, parameter :: FTIT         = 20
  integer, parameter :: FTITpc       = 22
  integer, parameter :: HBFILE       = 25
  integer, parameter :: SAXSs        = 26
  integer, parameter :: SAXSc        = 27
  
  character(len=20) :: conf_saddle, conf_final, conf_initial
  character(len=20) :: debug_status
  logical :: use_qbug ! If true, calculates the forces the first time and stop
  logical :: use_tit
  logical :: use_nma  !! Normal mode Analysis
  logical :: use_ics  !! Internal Coordiante Space
  logical :: use_back !! Backbone variables
  logical :: use_base !! Base+sugar variables
  logical :: use_val  !! Valence angle variables 
  logical :: saxs_print
  character(len=20) :: CHAIN_FILE
  character(len=20) :: MASTER_LOGFILE
  character(len=20) :: EVENTSLIST
  character(len=20) :: REFCONFIG
  character(len=3)  :: FINAL
  character(len=3)  :: SADDLE
  character(len=11) :: COUNTER
  character(len=11) :: RESTARTFILE
  character(len=13) :: REPLICAFILE

  character(len=20) :: SAXSFILE

  character(len=4)  :: PDB_EXT   = '.pdb'
  
  integer :: N_REPLICA      ! number of temperature replica (T exchange)
  integer :: N_E_REPLICA    ! number of energy replica (hamiltonian exchange)
  integer :: N_TEMP         ! number of temperature tempering)

  integer(8) :: n_step_exchange
  integer :: n_rate_saxs
  real(8), dimension(:), allocatable :: T_replica
  real(8), dimension(:), allocatable :: T_tempering
  real(8), dimension(:), allocatable :: F_tempering
  real(8), dimension(:), allocatable :: E_tempering
  real(8), dimension(:), allocatable :: N_tempering

  integer :: n_step_E_exchange
  integer :: last_id
  real(8), dimension(:), allocatable :: E_replica

  type t_conformations
    character(len=20) :: path
    character(len=20) :: logfile
    integer ::  id, counter
    real(8), dimension(:), allocatable :: pos, posref, vel
    real(8), dimension(:), allocatable :: temperatures,scales
    real(8) :: energy, temperature, energyscale,free_energy,E_energy,N_energy, score
  end type t_conformations

!===New minimizer+Internal coordinates
  real(8), dimension(:), allocatable, target :: gra,gras,poss
  real(8), dimension(:), allocatable, target :: var,ovr,svars,vars,ovrs
  real(8), dimension(:), allocatable, target :: scl
  integer :: nvar
!==================
  
  integer :: ntasks, taskid, Error
#ifdef MPI
  integer :: status(MPI_STATUS_SIZE)
#endif
end module defs

module numerical_defs
    implicit none
    double precision :: pi, rad2deg, rad2,deg2rad
    data pi/3.141592653589793d+00/
    !parameter (rad2deg=180.0d0/3.141592653d0)
    parameter (rad2deg = 57.29577951308232088d0)
    parameter (deg2rad = 0.0174532925d0)
    parameter (rad2 = rad2deg*rad2deg)

    integer MAXPRE, MAXNAT, MAXTTY, MAXXC, MAXPNB
    integer MAXBO, MAXTH, MAXPHI, MAXPAI
    parameter (MAXPRE = 1500)               !! maximum number of residues
    parameter (MAXNAT = MAXPRE*6)           !! maximum number of atoms
    parameter (MAXTTY = 50000)              !! maximum number of residue name types
    parameter (MAXXC = 3*MAXNAT)            !! maximum number of cart coord
    parameter (MAXPNB = 3*MAXPRE*MAXPRE)    !! max number of SC-SC interactions
    parameter (MAXBO  = MAXNAT)             !! maximum number of bonds
    parameter (MAXTH = MAXNAT*3)            !! maximum number of bond angles
    parameter (MAXPHI = MAXNAT*4)           !! maximum number of torsional angles
    parameter (MAXPAI = MAXNAT*(MAXNAT+1)/2)!! max number of nonbonded-pairs

    save
end module numerical_defs

module score
    implicit none
    double precision, dimension(272) :: score_prot
    double precision, dimension(47)  :: score_RNA

    save
end module score

module PBC_defs
    implicit none
    ! common/pbcBL/box_length, inv_box_length
    real(8) box_length, inv_box_length

    ! common/PBC_R/periodicBC,CM
    logical periodicBC,CM

    save
end module PBC_defs

module PDBtext
    use numerical_defs
    implicit none
    ! common/textt/text2,text3,text4
    ! common/nnumres/numres,Id_atom
    ! character*7 text2(MAXNAT)
    ! character*5 text3(MAXNAT)
    ! character*7 text4(MAXNAT)
    ! integer     numres(MAXNAT),Id_atom(MAXNAT)
    character(7) :: text2(MAXNAT)
    character(5) :: text3(MAXNAT)
    character(7) :: text4(MAXNAT)
    integer :: numres(MAXNAT),Id_atom(MAXNAT)
    
    save
end module PDBtext

module fragments
    use numerical_defs
    implicit none
    !common/frags/nfrag,lenfrag(MAXPRE),ichain(MAXNAT)
    !integer nfrag,lenfrag,ichain 
    integer :: nfrag, nfrag_prot, nfrag_rna,lenfrag(MAXPRE),ichain(MAXNAT)

    save
end module fragments


!================
! PROTEIN MODULES
!================
 
module system_defs_prot
    use numerical_defs
    implicit none
    ! misc1
    integer NRES,NBONH,NBONA,NTHETH,NTHETA,NPHIH,natom3,&
                NPHIA,NNB,NTYPES,MBONA,MTHETA,MPHIA

    !COMMON/MISC2/AMASS(MAXNAT),IAC(MAXNAT),NNO(MAXTTY)
    double precision amass(maxnat)
    integer iac(maxnat), nno(maxtty)


    !COMMON/MISC3/IPRES(MAXPRE),IGRAPH(MAXNAT),LABRES(MAXNAT)
    integer IPRES(MAXPRE)
    character(4) :: IGRAPH(MAXNAT),LABRES(MAXNAT)


    !  COMMON/NBPARA/CUT,SCNB,SCEE,IDIEL,DIELC
    ! --- SET SOME PARAMETERS
    double precision :: &
      CUT    = 100.0,&    !! NO CUTOFF
      SCNB   = 80,&       !! divide the 1-4 VDW interactions by 8.0
      SCEE   = 80.0,&     !! divide the 1-4 ELEC interactions by 2.0
      IDIEL  = 0.0,&      !! 2.0 dielectric constant epsilon = 2r
      DIELC  = 1.0        !! ...............................

    save
end module system_defs_prot

module param_defs_prot
    use numerical_defs
    implicit none

    !  COMMON/PARM1/RK(MAXBO),REQ(MAXBO),TK(MAXTH),TEQ(MAXTH),
    ! $             PK(MAXPHI),PN(MAXPHI),
    ! $             PHASE(MAXPHI),CN1(MAXTTY),CN2(MAXTTY),SOLTY(60),
    ! $             GAMC(MAXPHI),GAMS(MAXPHI),IPN(MAXPHI),FMN(MAXPHI)
    double precision :: &
      RK(MAXBO),REQ(MAXBO),TK(MAXTH),TEQ(MAXTH),&
      PK(MAXPHI),PN(MAXPHI),&
      PHASE(MAXPHI),CN1(MAXTTY),CN2(MAXTTY),SOLTY(60),&
      GAMC(MAXPHI),GAMS(MAXPHI),FMN(MAXPHI)
    integer :: IPN(MAXPHI)

    !  COMMON/ENER1/IB(MAXBO),JB(MAXBO),ICB(MAXBO),IBH(MAXBO),JBH(MAXBO),
    ! $             ICBH(MAXBO)
    integer :: &
      IB(MAXBO),JB(MAXBO),ICB(MAXBO),IBH(MAXBO),JBH(MAXBO),&
      ICBH(MAXBO)

    !  COMMON/ENER2/IT(MAXTH),JT(MAXTH),KT(MAXTH),ICT(MAXTH),ITH(MAXTH),
    ! $             JTH(MAXTH),KTH(MAXTH),ICTH(MAXTH)
    integer :: &
      IT(MAXTH),JT(MAXTH),KT(MAXTH),ICT(MAXTH),ITH(MAXTH),&
      JTH(MAXTH),KTH(MAXTH),ICTH(MAXTH)

    !  COMMON/ENER3/IP(MAXPHI),JP(MAXPHI),KP(MAXPHI),LP(MAXPHI),
    ! 1             ICP(MAXPHI)
    integer :: &
      IP(MAXPHI),JP(MAXPHI),KP(MAXPHI),LP(MAXPHI),ICP(MAXPHI)

    !  COMMON/ENER4/
    ! 1 IPH(MAXPHI),JPH(MAXPHI),KPH(MAXPHI),LPH(MAXPHI),ICPH(MAXPHI)
    integer :: &
      IPH(MAXPHI),JPH(MAXPHI),KPH(MAXPHI),LPH(MAXPHI),ICPH(MAXPHI)

    !  COMMON/NONBON/NUMEX(MAXNAT),NATEX(MAXTTY)
    integer :: NUMEX(MAXNAT),NATEX(MAXTTY)

    !  COMMON/PRMLIM/NUMBND,NUMANG,NPTRA,NPHB,NIMPRP
    integer :: NUMBND,NUMANG,NPTRA,NPHB


    !  COMMON/NPAIR/NPAIR2,IPAIR(MAXPAI),JPAIR(MAXPAI)
    integer :: NPAIR2,IPAIR(MAXPAI),JPAIR(MAXPAI)

    save
end module param_defs_prot

module cutoffs_prot
    implicit none
    double precision :: &
      rcut2_caca_scsc_out, rcut2_caca_scsc_in,&
      rcut2_hb_mcmc_out, rcut2_hb_mcmc_in,&
      rcut2_4b_out, rcut2_4b_in,&
      rcut2_lj_out, rcut2_lj_in

    save
end module cutoffs_prot

module charge_prot
    use numerical_defs
    implicit none
    !common/charge/cg
    !CG(MAXNAT)
    double precision :: CG(MAXNAT)

    save
end module charge_prot



!============
! RNA MODULES
!============
 
module system_defs_RNA
    use numerical_defs
    implicit none
    ! misc1
    integer NRES,NBONH,NBONA,NTHETH,NTHETA,NPHIH,natom3,&
                NPHIA,NNB,NTYPES,MBONA,MTHETA,MPHIA

    !COMMON/MISC2/AMASS(MAXNAT),IAC(MAXNAT),NNO(MAXTTY)
    double precision amass(maxnat)
    integer iac(maxnat), nno(maxtty)


    !COMMON/MISC3/IPRES(MAXPRE),IGRAPH(MAXNAT),LABRES(MAXNAT)
    integer IPRES(MAXPRE)
    character(4) :: IGRAPH(MAXNAT),LABRES(MAXNAT)


    !  COMMON/NBPARA/CUT,SCNB,SCEE,IDIEL,DIELC
    ! --- SET SOME PARAMETERS
    double precision :: &
      CUT    = 100.0,&    !! NO CUTOFF
      SCNB   = 80,&       !! divide the 1-4 VDW interactions by 8.0
      SCEE   = 80.0,&     !! divide the 1-4 ELEC interactions by 2.0
      IDIEL  = 0.0,&      !! 2.0 dielectric constant epsilon = 2r
      DIELC  = 1.0        !! ...............................

    save
end module system_defs_RNA

module param_defs_RNA
    use numerical_defs
    implicit none

    !  COMMON/PARM1/RK(MAXBO),REQ(MAXBO),TK(MAXTH),TEQ(MAXTH),
    ! $             PK(MAXPHI),PN(MAXPHI),
    ! $             PHASE(MAXPHI),CN1(MAXTTY),CN2(MAXTTY),SOLTY(60),
    ! $             GAMC(MAXPHI),GAMS(MAXPHI),IPN(MAXPHI),FMN(MAXPHI)
    double precision :: &
      RK(MAXBO),REQ(MAXBO),TK(MAXTH),TEQ(MAXTH),&
      PK(MAXPHI),PN(MAXPHI),&
      PHASE(MAXPHI),CN1(MAXTTY),CN2(MAXTTY),SOLTY(60),&
      GAMC(MAXPHI),GAMS(MAXPHI),FMN(MAXPHI)
    integer :: IPN(MAXPHI)

    !  COMMON/ENER1/IB(MAXBO),JB(MAXBO),ICB(MAXBO),IBH(MAXBO),JBH(MAXBO),
    ! $             ICBH(MAXBO)
    integer :: &
      IB(MAXBO),JB(MAXBO),ICB(MAXBO),IBH(MAXBO),JBH(MAXBO),&
      ICBH(MAXBO)

    !  COMMON/ENER2/IT(MAXTH),JT(MAXTH),KT(MAXTH),ICT(MAXTH),ITH(MAXTH),
    ! $             JTH(MAXTH),KTH(MAXTH),ICTH(MAXTH)
    integer :: &
      IT(MAXTH),JT(MAXTH),KT(MAXTH),ICT(MAXTH),ITH(MAXTH),&
      JTH(MAXTH),KTH(MAXTH),ICTH(MAXTH)

    !  COMMON/ENER3/IP(MAXPHI),JP(MAXPHI),KP(MAXPHI),LP(MAXPHI),
    ! 1             ICP(MAXPHI)
    integer :: &
      IP(MAXPHI),JP(MAXPHI),KP(MAXPHI),LP(MAXPHI),ICP(MAXPHI)

    !  COMMON/ENER4/
    ! 1 IPH(MAXPHI),JPH(MAXPHI),KPH(MAXPHI),LPH(MAXPHI),ICPH(MAXPHI)
    integer :: &
      IPH(MAXPHI),JPH(MAXPHI),KPH(MAXPHI),LPH(MAXPHI),ICPH(MAXPHI)

    !  COMMON/NONBON/NUMEX(MAXNAT),NATEX(MAXTTY)
    integer :: NUMEX(MAXNAT),NATEX(MAXTTY)

    !  COMMON/PRMLIM/NUMBND,NUMANG,NPTRA,NPHB,NIMPRP
    integer :: NUMBND,NUMANG,NPTRA,NPHB


    !  COMMON/NPAIR/NPAIR2,IPAIR(MAXPAI),JPAIR(MAXPAI)
    integer :: NPAIR2,IPAIR(MAXPAI),JPAIR(MAXPAI)

    save
end module param_defs_RNA

module cutoffs_RNA
    implicit none
    double precision :: &
      rcut2_caca_scsc_out, rcut2_caca_scsc_in,&
      rcut2_hb_mcmc_out, rcut2_hb_mcmc_in,&
      rcut2_4b_out, rcut2_4b_in,&
      rcut2_lj_out, rcut2_lj_in

    save
end module cutoffs_RNA

module charge_RNA
    use numerical_defs
    implicit none
    !common/charge/cg
    !CG(MAXNAT)
    double precision :: CG(MAXNAT)

    save
end module charge_RNA

module int_defs
  
  implicit none
  integer, dimension(:), allocatable :: iseg1    !  number of the first residue of the segment
  integer, dimension(:), allocatable :: iseg2    !  number of the last residue of the segment
  integer, dimension(:), allocatable :: ibseg    !  number of the chain of the segment
  integer, dimension(:), allocatable :: lc       !  last atom number of each body
  integer, dimension(:), allocatable :: ipiv     !  pivoting P for each body
  integer, dimension(:), allocatable :: kch      !  last residue index of each body
  integer, dimension(:), allocatable :: nbasepp  !  number of atoms to move for each base
  integer, dimension(:), allocatable :: nbbvar   !  number of backbone variable for the 5' of each body
  logical,  dimension(:), allocatable ::  log3t  ! Logical for the 3' terminus
  logical,  dimension(:), allocatable ::  log5t  ! Logical for the 5' terminus
  integer :: nloop
  
end module int_defs

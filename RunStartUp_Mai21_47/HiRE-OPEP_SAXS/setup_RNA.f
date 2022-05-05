        module restraint_params
        integer :: Nrests,  Nposres
! Distance restraints
        integer, allocatable, dimension(:) :: resti, restj, trest
        double precision, allocatable, dimension(:) :: restk, restl,
     &   drestl

! Position restraints
        integer, allocatable, dimension(:) :: pri
        double precision, allocatable, dimension(:) :: prk
        double precision, allocatable, dimension(:,:) :: prx, prdx
        save
      end module restraint_params

      module RNAnb
        double precision, allocatable, dimension(:,:) :: nbcoef,nbct2
        double precision, allocatable, dimension(:) :: chatm         ! Br1
        double precision, allocatable, dimension(:) :: chrg
        integer, allocatable, dimension(:,:) :: nbscore
        character basenum(4)
        double precision alpa(4,4), alpb(4,4)
        integer bcoef(4,4)
        save
      end module RNAnb

      module rnabase
          integer, allocatable, dimension(:) :: blist,btype,bprot,bocc
          real, allocatable, dimension(:) :: bpch
          save
      end module rnabase

      module RNAHBparams
          double precision, dimension(6,4,4) :: dREF, alpam, alpbm, 
     &    calpam, salpam, calpbm, salpbm
          double precision, dimension(6,4,4,4,4) :: s !Br2 make BP dependent on protonation
          integer :: Nparam(4,4)
          double precision :: wc, wcCanonic, noWc, prot, tit, noWCq, sha
          double precision, dimension(4,3) :: planarityDistEq
          save
      end module RNAHBparams

!Br2Fr add new module for titration
      module RNAbaseTit
          integer :: fatortit(1000)
      end module RNAbaseTit 

      module RNAtorstype

             integer, allocatable, dimension(:) :: nphitype, nthettype

      end module RNAtorstype
      
      
      SUBROUTINE INITIALISE_RNA(P_B_C,BoxL,CofM,NATOM_TOT,
     &                              XFORT,AMASSES,
     &                              ATOMIC_TYPE)
C
C   1.READ THE TOPOLOGY FILE (INTERNAL COORDINATES,
C     PARAMETERS OF THE OPEP FUNCTION) AND THE STARTING
C     PDB FILE: (UNIT 14).
C
C   2.READ THE PARAMETERS FOR THE SIDE-CHAIN SIDE-CHAIN INTERACTIONS
C     (UNIT 28).
C
C   3.COMPUTE ENERGY and FORCES (NO CUTOFF FOR THE NONBONDED-INTERACTIONS. 
C
C     THIS PROGRAM WAS WRITTEN BY P. DERREUMAUX while at IBPC, UPR 9080
C     CNRS, PARIS, FRANCE.  DECEMBER 1998
C
      use defs, only: N_prot, N_RNA, use_qbug
      use numerical_defs
      use system_defs_RNA
      use param_defs_RNA
      use fragments
      use restraint_params
      use RNAnb
      use RNAHBparams
      use fileio
      use rnabase
      use score
      use cutoffs_RNA
      use charge_RNA
      use PDBtext
      use PBC_defs 
      use RNAtorstype

      implicit none 

      integer, intent(in) :: NATOM_TOT
      double precision :: xfort(3*NATOM_TOT), amasses(NATOM_TOT)
      character*5  ATOMIC_TYPE(NATOM_TOT)

C---------------------------------------------------------------------
      logical P_B_C,CofM, QDET
      real(8) BoxL
C---------------------------------------------------------------------

      double precision X(MAXXC)


      real(8) score_temp

      character(50) dummy, dummych

      integer i, j, idum, ifer, ikl, idatom
      integer nf, jn, i3, j3, k3t, l3t, it1, it2, it0, it3

      logical fileexists

C________________________________________________________________________
      character*30 path_for_inits
      logical      single_conf_init
      periodicBC = P_B_C
      CM = CofM
      write(*,*) 'PBC and CM are : ',P_B_C, CofM
C________________________________________________________________________



C---  read TOPOLOGY file
      nf=21
C     This was the original configuration
      OPEN(UNIT=nf,FILE="parametres_RNA.top",status="old")
      CALL RDTOP_RNA(nf,CG)
      CLOSE(UNIT=nf)
 

C---  read INPUT PDB file
      CALL init_remd(path_for_inits,single_conf_init)
      if(.not. single_conf_init)  then 
         nf=41
         OPEN(UNIT=nf,FILE=path_for_inits,status="old")
      else  
         nf=14
         OPEN(UNIT=nf,FILE="conf_initiale_RNA.pdb",status="old")
      endif
      CALL RDPDB(nf,N_RNA,X,text2(N_prot+1:N_RNA),
     &  Id_atom(N_prot+1:N_RNA),text3(N_prot+1:N_RNA),
     &  text4(N_prot+1:N_RNA),numres(N_prot+1:N_RNA))
      close(nf)

!! Br2 : read ions from PDB and add them to RNA list
!! Br2 : define masses for ions (not from top file)

      !print *, 'pippo reading pdb'

c---  NEW JANV05
C---- Set-up the vector associating each atom with a fragment of the protein
C----    length_fragment: table containing the length of each fragment
C----    ichain : the fragment to which atom idatom belongs to

C---  read ichain file
      nf=55
      OPEN(UNIT=nf,FILE="ichain_RNA.dat",status="old")
      read(55,*)  nfrag_rna
      NFRAG = NFRAG+nfrag_rna
      do i=nfrag_prot+1, nfrag
         read(55,*) idum, lenfrag(i)
      enddo
      idatom = N_prot+1
      do i=nfrag_prot+1, NFRAG
         do j=1, lenfrag(i)
            ichain(idatom) = i
            idatom = idatom +1
         end do
      end do
!      read(55,*) N_ION                   !Br2 : add number of ions as last line of ichain
      close(55)

!! Br2 : adjust ichain to account for ions ? possibly not needed, check Ebond ecc.

       !print *, 'pippo reading ichain'

C---  read the restraints
      inquire(file="restraints.dat", exist=fileexists)
      if(fileexists .eqv. .true.) then
c      if(testfile("restraints.dat") .eqv. .true.) then
        open(unit=55, file="restraints.dat", status="old")
        read(55,*) Nrests, Nposres
        allocate(resti(Nrests), restj(Nrests), trest(Nrests))
        allocate(restk(Nrests), restl(Nrests), drestl(Nrests))
        allocate(pri(Nposres), prk(Nposres))
        allocate(prx(3,Nposres), prdx(3,Nposres))
!       read a blank line for comments at the top of the list
        read(55, *) 
        read(55, *) (resti(i), restj(i), restk(i), restl(i),
     $   drestl(i), trest(i), i = 1,Nrests)
!       read a blank line for comments at the top of the posres list
        if(Nposres .gt. 0) then
          read(55, *) 
          read(55, *) (pri(i), prk(i), prx(1,i), prx(2,i), prx(3,i),
     $       prdx(1,i), prdx(2,i), prdx(3,i), i = 1,Nposres)
        endif
      else
        Nrests = 0
        Nposres = 0
      endif

c        read(55,'(2i4, f15.10)') (resti(I), restj(i), restlens(i),
c     $      i = 1,Nrests)
      if( Nrests .gt. 0) then
        write (*,*) "Distance restraints"
        write (*,"(2a5,3a10,a5)") "ri","rj","k","length","dl","type"
        do i = 1, Nrests
          write(*,"(2i5, 3f10.6, i2)") resti(i), restj(i), restk(i), 
     &      restl(i), drestl(i), trest(i)
        enddo
      endif
      if(Nposres .gt. 0) then
        write (*,*) "Position restraints"
        write (*,"(2a5, 6a10)") "pi", "pk", "px", "py", "pz",
     &     "pdx", "pdy", "pdz"
        do i = 1, Nposres
            write(*,"(i5, 7f10.6)") pri(i), prk(i),
     &          prx(1,i), prx(2,i), prx(3,i),
     &          prdx(1,i), prdx(2,i), prdx(3,i)
        enddo
      endif


C---  read list of bases, for RNA
      if (.not. allocated(blist)) then
        allocate( blist(NRES))
        allocate( btype(NRES))
        allocate( bprot(NRES))
        allocate( bocc(NRES))
        allocate( bpch(NRES))
      endif
      open(unit=38,file="bblist.dat",status="old")   
      read(38,*) (blist(i), btype(i), bprot(i), i =1, NRES) 
      close(38)
      !print *, 'pippo reading baselist'

!---  read weight file for RNA
      OPEN(UNIT=55,FILE="scale_RNA.dat",status="old")
      do i=1,47
         read(nf,'(i4,f15.10,A)') ikl, score_temp, dummy
         score_RNA(i)=dble(score_temp)
      enddo
      close(55)

      ! NTYPES is the number of particle types, found in the .top file
      ! Currently, it's taken to be :
      ! 1: C5*  2: O5*  3: P  4: CA  5: CY
      ! 6,7: G1,2  8,9: A1,2  10: U1  11: C1
      ! 12: D 13: MG 14: NA  15:CL
      if (.not. allocated(nbcoef)) then
        allocate( nbcoef(NTYPES, NTYPES))
        allocate( nbct2(NTYPES, NTYPES))
        allocate( nbscore(NTYPES, NTYPES))
        allocate( chrg(NTYPES))
        allocate( chatm(N_RNA))
      endif
      nbcoef = 1.0
      nbct2 = 3.6 ! 3.65 ! 3.6 3.2
      nbscore = 4
      chrg = 0
      chatm = 0                                !Br1

!      ! C4 are bigger
      nbct2(4,:) = 4.0
      nbct2(:,4) = 4.0
!      ! CY are bigger
!      nbct2(5,:) = 4.0
!      nbct2(:,5) = 4.0
!      ! bases beads size
      nbct2(6:11,6:11) = 3.2 !3.25  !3.2
      ! P-P interactions
      nbcoef(3,3) = 1.0
      nbct2(3,3) = 4.0   !4.0
      nbscore(3,3) = 5
!     dummy residues are larger
      if (NTYPES .eq. 11) then
        nbct2(12,:) = 8
        nbct2(:,12) = 8
      endif

      ! Mg++ interactions
      nbscore(3,13) = 5
      nbscore(13,3) = 5
      nbscore(13,13) = 5
      nbscore(13,14) = 5
      nbscore(13,15) = 5
      chrg(3) = -1
      chrg(13) = 2

      ! Na+ interactions
      nbscore(3,14) = 5
      nbscore(14,3) = 5
      nbscore(14,13) = 5
      nbscore(14,15) = 4
      chrg(14) = 1

      ! Cl- interactions
      nbscore(3,15) = 5
      nbscore(15,3) = 5
      nbscore(15,13) = 5
      nbscore(15,14) = 4
      chrg(15) = -1 

      !print *, 'pippo declare particles'

!      do i = 1,N_RNA                                  !Br1
!         chatm(i) = chrg(IAC(i)   )                   !Br1
!      enddo

      !print *, N_RNA
      open(unit=39,file="chargeatm_RNA.dat",status="old")   ! Br2
      !print *, 'open charges'
      read(39,*) (dummych, chatm(i), i =1, N_RNA)
      !print *, 'pippo reading charge 1x1'
      close(39)

c$$$      do i = 1,N_RNA                                   !Br2
c$$$         print *, chatm(i), chrg(IAC(i)   )            !Br2
c$$$      enddo
     
       !print *, 'pippo reading charge'

!        do i = 1,N_RNA                                   !Br2
!           print *, chatm(i), chrg(IAC(i)   )            !Br2
!        enddo

c$$$      do i = 1, NTYPES
c$$$        do j = 1, NTYPES
c$$$          nbcoef(i,j) = nbcoef(i,j)*score_RNA(nbscore(i,j))
c$$$        enddo
c$$$      enddo

      alpam = 2
      alpbm = 2
      bcoef = 35

      basenum = (/ 'G', 'A', 'C', 'U' /) 

      !                G   A   C  U
      bcoef(1,:) = (/ 34, 20, 19, 21 /)  ! ???? what's this ??? ever used?
      bcoef(2,:) = (/ 20, 32, 22, 18 /)
      bcoef(3,:) = (/ 19, 22, 33, 27 /)
      bcoef(4,:) = (/ 21, 18, 27, 35 /)

      call fillHBparams()
      alpam = alpam+pi
      alpbm = alpbm+pi
      calpam = cos(alpam)
      salpam = sin(alpam)
      calpbm = cos(alpbm)
      salpbm = sin(alpbm)

      !print *, 'pippo loaded'
      
      QDET = .FALSE. .or. use_qbug

      allocate(nthettype(NTHETA))
      do jn = 1,NTHETA
         i3 = it(jn)/3 +1
         j3 = jt(jn)/3 +1
         k3t = kt(jn)/3 +1
         it0 = iac(i3)
         it1 = iac(j3)
         it2 = iac(k3t)

         if (it0.eq.4.and.it1.eq.5.and.(it2.eq.6.or.it2.eq.8
     &    .or.it2.eq.10.or.it2.eq.11))then
            nthettype(jn)=0

         else if (it0.eq.5.and.(it1.eq.6.and.it2.eq.7.or.it1.
     &    eq.8.and.it2.eq.9))then
            nthettype(jn)=1
         
         else if (it0.eq.3.and.it1.eq.2.and.it2.eq.1)then
            nthettype(jn)=2
         
         else if (it0.eq.2.and.it1.eq.1.and.it2.eq.4)then
            nthettype(jn)=3
         
         else if (it0.eq.1.and.it1.eq.4.and.it2.eq.3)then
            nthettype(jn)=4
         
         else if (it0.eq.4.and.it1.eq.3.and.it2.eq.2)then
            nthettype(jn)=5
         
         else if (it0.eq.1.and.it1.eq.4.and.it2.eq.5)then
            nthettype(jn)=6
         
         else if (it0.eq.5.and.it1.eq.4.and.it2.eq.3)then
            nthettype(jn)=7
         
         endif
      enddo



      allocate(nphitype(NPHIA))
      do jn = 1,NPHIA                ! torsion info backbone or base
         i3 = ip(jn)/3+1
         j3 = jp(jn)/3+1
         k3t = kp(jn)/3+1
         l3t = lp(jn)/3+1
         it0=iac(i3)
         it1=iac(j3)
         it2=iac(k3t)
         it3=iac(l3t)

         if(it0.eq.4.and.it1.eq.5.and.(it2.eq.6.or.it2.eq.8)
     &    .and.(it3.eq.7.or.it3.eq.9))then
            nphitype(jn)=0

         else if(it0.eq.4.and.(it1.eq.6.or.it1.eq.8).and.(it2.eq.7
     &    .or.it2.eq.9).and.it3.eq.5)then
            nphitype(jn)=1

         else if(it0.eq.1.and.it1.eq.4.and.it2.eq.5.and.
     &    (it3.eq.6.or.it3.eq.8.or.it3.eq.10.or.it3.eq.11))then
            nphitype(jn)=2

         else if(it0.eq.3.and.it1.eq.4.and.it2.eq.5.and.
     &    (it3.eq.6.or.it3.eq.8.or.it3.eq.10.or.it3.eq.11))then
            nphitype(jn)=3

         else if(it0.eq.1.and.it1.eq.4.and.it2.eq.3.and.it3.eq.2)then
            nphitype(jn)=4

         else if(it0.eq.5.and.it1.eq.4.and.it2.eq.3.and.it3.eq.2)then
            nphitype(jn)=5

         else if(it0.eq.2.and.it1.eq.1.and.it2.eq.4.and.it3.eq.3)then
            nphitype(jn)=6

         else if(it0.eq.2.and.it1.eq.1.and.it2.eq.4.and.it3.eq.5)then
            nphitype(jn)=7

         else if(it0.eq.3.and.it1.eq.2.and.it2.eq.1.and.it3.eq.4)then
            nphitype(jn)=8

         else if(it0.eq.4.and.it1.eq.3.and.it2.eq.2.and.it3.eq.1)then
            nphitype(jn)=9

         endif
          if(QDET) then
           write(7,'(8i4)') jn,i3,j3,k3t,l3t,it1,it2,nphitype(jn) 
          endif
      enddo   
!        write(*,"(2i5, 3f10.6, i2)") resti(i), restj(i), restk(i), 
C---------------------------------------------------------------

      do i=1, 3*N_RNA
        xfort(N_prot*3+i) = x(i)
      enddo

      do i=1, N_RNA
        amasses(N_prot+i) = amass(i)
        atomic_type(N_prot+i) = text3(N_prot+i)
      end do 

      open(57, file="cutoffs_RNA.dat")
      read(57, *)  rcut2_caca_scsc_out, rcut2_caca_scsc_in,
     $             rcut2_hb_mcmc_out, rcut2_hb_mcmc_in,
     $             rcut2_4b_out, rcut2_4b_in,
     $             rcut2_lj_out, rcut2_lj_in
      close(57)

c     We copy the box length and define the inverse box-length
      box_length = BoxL
      inv_box_length = 1.0d0 / box_length

      return
      end
      
C     --------------------------------------------------------------

      subroutine fillHBparams()
          common/MOLECULETYPE/ molecule_type
          character(len=20) molecule_type

          if(molecule_type=="DNA") then
            call fillDNAHBparams()
          else if(molecule_type=="RNA") then
            call fillRNAHBparams()
          endif

      end subroutine fillHBparams

C     --------------------------------------------------------------

      subroutine fillRNAHBparams()
          use RNAHBparams
          use score

          implicit none
          double precision z

          wc = score_RNA(18)
          wcCanonic = score_RNA(19)
          noWc = score_RNA(20)
          tit = score_RNA(21)
          noWCq = score_RNA(22)
          sha = score_RNA(26)
          z=0.0d0

          
          
!         A-A
          dREF(1:3,2,2) =  (/ 5.63, 6.84, 5.92/)       !WWt, HHt, HWt                    
          alpam(1:3,2,2) = (/ 2.40, 1.05, 2.72/)
          alpbm(1:3,2,2) = (/ 2.40, 1.08, 1.29/)
          s(1:3,2,2,1,1) = (/2*wc, 2*noWc, 2*noWc/)        !Br2  1,1 -> q1=0, q2=0
          s(1:3,2,2,1,2) = (/    z,2*noWcq, z/)            !Br2  1,2 -> q1=0, q2=1
          s(1:3,2,2,2,1) = (/    z,2*noWcq, 2*noWcq/)      !Br2  2,1 -> q1=1, q2=0
          s(1:3,2,2,2,2) = (/    z,2*noWcq, z/)            !Br2  2,2 -> q1=1, q2=1
          s(1:3,2,2,4,1) = (/ z, 2*noWc, 2*noWc/)          !4->DMS/CforMCT
	  s(1:3,2,2,4,2) = (/ z, 2*noWc, 2*noWc/)          !4->DMS/CMCT
          s(1:3,2,2,1,4) = (/ z, 2*noWc, 2*noWc/)	   !4->DMS/CMCT
	  s(1:3,2,2,2,4) = (/ z, 2*noWc, 2*noWc/)	   !4->DMS/CMCT
          s(1:3,2,2,3,1) = (/ 2*wc*sha, 2*noWc, 2*noWc/)   !3->SHAPE
	  s(1:3,2,2,3,2) = (/ 2*wc*sha, 2*noWc, 2*noWc/)   !3->SHAPE
          s(1:3,2,2,1,3) = (/ 2*wc*sha, 2*noWc, 2*noWc/)   !3->SHAPE
	  s(1:3,2,2,2,3) = (/ 2*wc*sha, 2*noWc, 2*noWc/)   !3->SHAPE
	  Nparam(2,2) = 3                !Br3  CHECK angle sign when order is inverted when one is negative

!         A-C
          dREF(1:3,2,3) = (/7.26, 5.78, 4.78/)         !WSc, HWt, +Wc  
          alpam(1:3,2,3) = (/2.36, 1.36, 2.64/)
          alpbm(1:3,2,3) = (/0.75, 2.30, 1.78/)   
          s(1:3,2,3,1,1) = (/ 2*noWc,  2*noWc, z /) !Br2  1,1 -> q1=0, q2=0
          s(1:3,2,3,1,2) = (/ 2*noWcq, z, z/)                 !Br2  1,2 -> q1=0, q2=1
          s(1:3,2,3,2,1) = (/ z, 2*noWcq, 2*tit/)            !Br2  2,1 -> q1=1, q2=0
          s(1:3,2,3,2,2) = (/ z, z, z/)                      !Br2  2,2 -> q1=1, q2=1
          s(1:3,2,3,4,1) = (/ 2*noWc,  2*noWc, z /)          !4->DMS/CMCT
	  s(1:3,2,3,4,2) = (/ 2*noWc,  2*noWc, z /)          !4->DMS/CMCT
          s(1:3,2,3,1,4) = (/ 2*noWc,  2*noWc, z /)          !4->DMS/CMCT
	  s(1:3,2,3,2,4) = (/ 2*noWc,  2*noWc, z /)          !4->DMS/CMCT		
          s(1:3,2,3,3,1) = (/ 2*noWc,  2*noWc, z /) 	     !3->SHAPE
	  s(1:3,2,3,3,2) = (/ 2*noWc,  2*noWc, z /) 	     !3->SHAPE
          s(1:3,2,3,1,3) = (/ 2*noWc,  2*noWc, z /) 	     !3->SHAPE
	  s(1:3,2,3,2,3) = (/ 2*noWc,  2*noWc, z /) 	     !3->SHAPE
          Nparam(2,3) = 3
          dREF(:,3,2) = dREF(:,2,3)
          alpam(:,3,2) = alpbm(:,2,3)
          alpbm(:,3,2) = alpam(:,2,3)
          s(:,3,2,1,1) = s(:,2,3,1,1)
          s(:,3,2,1,2) = s(:,2,3,2,1)
          s(:,3,2,2,1) = s(:,2,3,1,2)
          s(:,3,2,2,2) = s(:,2,3,2,2)
          s(:,3,2,4,1) = s(:,2,3,4,1)
	  s(:,3,2,4,2) = s(:,2,3,4,2)
          s(:,3,2,1,4) = s(:,2,3,1,4)
	  s(:,3,2,2,4) = s(:,2,3,2,4)
          s(:,3,2,3,1) = s(:,2,3,3,1)
	  s(:,3,2,3,2) = s(:,2,3,3,2)
          s(:,3,2,1,3) = s(:,2,3,1,3)
	  s(:,3,2,2,3) = s(:,2,3,2,3)
          Nparam(3,2) = Nparam(2,3)

!         A-G
          dREF(1:4,2,1) = (/ 4.88, 6.63, 6.17, 6.05/)         !  WW_c, HSt, sst, +Hc  
          alpam(1:4,2,1) = (/ 3.04, 1.19, -2.01, 2.62/)
          alpbm(1:4,2,1) = (/ 2.58, -1.69, -1.57, 0.89/)
          s(1:4,2,1,1,1) = (/ 2*wc, 2*noWc, 2*noWc, z/)        !Br2  1,1 -> q1=0, q2=0
          s(1:4,2,1,1,2) = (/ z, 2*noWcq, 2*noWcq, z/)            !Br2  1,2 -> q1=0, q2=1
          s(1:4,2,1,2,1) = (/ z, 2*noWcq, 2*noWcq, 2*tit/)      !Br2  2,1 -> q1=1, q2=0
          s(1:4,2,1,2,2) = (/ z, 2*noWcq, 2*noWcq, 2*tit/)      !Br2  2,2 -> q1=1, q2=1
          s(1:4,2,1,4,1) = (/ z, 2*noWc, 2*noWc, z/) 		!4->DMS/CMCT
	  s(1:4,2,1,4,2) = (/ z, 2*noWc, 2*noWc, z/) 		!4->DMS/CMCT
          s(1:4,2,1,1,4) = (/ z, 2*noWc, 2*noWc, z/)		!4->DMS/CMCT
	  s(1:4,2,1,2,4) = (/ z, 2*noWc, 2*noWc, z/)		!4->DMS/CMCT
          s(1:4,2,1,3,1) = (/ 2*wc*sha, 2*noWc, 2*noWc, z/)	!3->SHAPE
	  s(1:4,2,1,3,2) = (/ 2*wc*sha, 2*noWc, 2*noWc, z/)	!3->SHAPE
	  s(1:4,2,1,1,3) = (/ 2*wc*sha, 2*noWc, 2*noWc, z/)	!3->SHAPE
	  s(1:4,2,1,2,3) = (/ 2*wc*sha, 2*noWc, 2*noWc, z/)	!3->SHAPE
          Nparam(2,1) = 4


          dREF(:,1,2) = dREF(:,2,1)
          alpam(:,1,2) = alpbm(:,2,1)
          alpbm(:,1,2) = alpam(:,2,1)
          s(:,1,2,1,1) = s(:,2,1,1,1)
          s(:,1,2,1,2) = s(:,2,1,2,1)
          s(:,1,2,2,1) = s(:,2,1,1,2)
          s(:,1,2,2,2) = s(:,2,1,2,2)
          s(:,1,2,4,1) = s(:,2,1,4,1)
	  s(:,1,2,4,2) = s(:,2,1,4,2)
          s(:,1,2,1,4) = s(:,2,1,1,4)
	  s(:,1,2,2,4) = s(:,2,1,2,4)
	  s(:,1,2,3,1) = s(:,2,1,3,1)
	  s(:,1,2,3,1) = s(:,2,1,3,1)
	  s(:,1,2,1,3) = s(:,2,1,1,3)
	  s(:,1,2,2,3) = s(:,2,1,2,3)
          Nparam(1,2) = Nparam(2,1)


!         A-U
          dREF(1:3,2,4) = (/ 4.92, 6.00, 6.00/)        ! WWc, HWt, HWc  !lal0420 modified 6
          alpam(1:3,2,4) = (/ 2.84, 0.91, 0.93/)
          alpbm(1:3,2,4) = (/ 2.36, 1.82, 2.46/)
          s(1:3,2,4,1,1) = (/ 2.2*wcCanonic, 2*noWc, 2*noWc/ )     !Br2  1,1 -> q1=0, q2=0      !OKKIO : cambiato a mano 14.4 -> 16!!!
          s(1:3,2,4,1,2) = (/ z, z, z/)                     !Br2  1,2 -> q1=0, q2=1
          s(1:3,2,4,2,1) = (/ z, 2*noWcq, 2*noWcq/)          !Br2  2,1 -> q1=1, q2=0
          s(1:3,2,4,2,2) = (/ z, z, z/)               !Br2  2,2 -> q1=1, q2=1
          s(1:3,2,4,4,1) = (/ z, 2*noWc, 2*noWc/ )			!4->DMS/CMCT
	  s(1:3,2,4,4,2) = (/ z, 2*noWc, 2*noWc/ )			!4->DMS/CMCT
          s(1:3,2,4,1,4) = (/ z, 2*noWc, 2*noWc/ )			!4->DMS/CMCT
	  s(1:3,2,4,2,4) = (/ z, 2*noWc, 2*noWc/ )			!4->DMS/CMCT
          s(1:3,2,4,3,1) = (/ 2.2*wcCanonic*sha, 2*noWc, 2*noWc/ )	!3->SHAPE
	  s(1:3,2,4,3,2) = (/ 2.2*wcCanonic*sha, 2*noWc, 2*noWc/ )	!3->SHAPE
	  s(1:3,2,4,1,3) = (/ 2.2*wcCanonic*sha, 2*noWc, 2*noWc/ )	!3->SHAPE
	  s(1:3,2,4,2,3) = (/ 2.2*wcCanonic*sha, 2*noWc, 2*noWc/ )	!3->SHAPE
          Nparam(2,4) = 3
          
          dREF(:,4,2) = dREF(:,2,4)
          alpam(:,4,2) = alpbm(:,2,4)
          alpbm(:,4,2) = alpam(:,2,4)
          s(:,4,2,1,1) = s(:,2,4,1,1)
          s(:,4,2,1,2) = s(:,2,4,2,1)
          s(:,4,2,2,1) = s(:,2,4,1,2)
          s(:,4,2,2,2) = s(:,2,4,2,2)
          s(:,4,2,4,1) = s(:,2,4,4,1)
	  s(:,4,2,4,2) = s(:,2,4,4,2)
          s(:,4,2,1,4) = s(:,2,4,1,4)
	  s(:,4,2,2,4) = s(:,2,4,2,4)
	  s(:,4,2,3,1) = s(:,2,4,3,1)
	  s(:,4,2,3,2) = s(:,2,4,3,2)
	  s(:,4,2,1,3) = s(:,2,4,1,3)
	  s(:,4,2,2,3) = s(:,2,4,2,3)
          Nparam(4,2) = Nparam(2,4)

!          C-C
!           dREF(1:1,3,3) =  (/ 4.91/)  !Br3 NO significant CC pairing
!           alpam(1:1,3,3) = (/ 2.22/)
!           alpbm(1:1,3,3) = (/ 2.24/)
!           s(1:1,3,3,1,1) = (/ z /) !Br2  1,1 -> q1=0, q2=0
!           s(1:1,3,3,1,2) = (/ z /) !Br2  1,2 -> q1=0, q2=1
!           s(1:1,3,3,2,1) = (/ z /) !Br2  2,1 -> q1=1, q2=0
!           s(1:1,3,3,2,2) = (/ z /) !Br2  2,2 -> q1=1, q2=1
!           Nparam(3,3) = 1

!         C-G
          dREF(1:3,3,1) = (/ 4.75, 5.28, 5.68/)              ! WWc, WWt, +Wc
          alpam(1:3,3,1) = (/ 2.17, 1.64, 1.92/)
          alpbm(1:3,3,1) = (/ 2.71, -3.07, 2.34/)
          s(1:3,3,1,1,1) = (/ 2.6*wcCanonic, 2*wc, z/)       !Br2  1,1 -> q1=0, q2=0
          s(1:3,3,1,1,2) = (/ z, z, z/)                      !Br2  1,2 -> q1=0, q2=1
          s(1:3,3,1,2,1) = (/ z, z, 1*tit/)                  !Br2  2,1 -> q1=1, q2=0
          s(1:3,3,1,2,2) = (/ z, z, z/)                      !Br2  2,2 -> q1=1, q2=1
          s(1:3,3,1,4,1) = (/ z, z, z/) 		       !4->DMS/CMCT
	  s(1:3,3,1,4,2) = (/ z, z, z/) 		       !4->DMS/CMCT
          s(1:3,3,1,1,4) = (/ z, z, z/)			       !4->DMS/CMCT
	  s(1:3,3,1,2,4) = (/ z, z, z/)			       !4->DMS/CMCT
	  s(1:3,3,1,3,1) = (/ 2.6*wcCanonic*sha, 2*wc*sha, z/) !3->SHAPE
	  s(1:3,3,1,3,2) = (/ 2.6*wcCanonic*sha, 2*wc*sha, z/) !3->SHAPE
	  s(1:3,3,1,1,3) = (/ 2.6*wcCanonic*sha, 2*wc*sha, z/) !3->SHAPE
	  s(1:3,3,1,2,3) = (/ 2.6*wcCanonic*sha, 2*wc*sha, z/) !3->SHAPE
          Nparam(3,1) = 3
          
          dREF(:,1,3) = dREF(:,3,1)
          alpam(:,1,3) = alpbm(:,3,1)
          alpbm(:,1,3) = alpam(:,3,1)
          s(:,1,3,1,1) = s(:,3,1,1,1)
          s(:,1,3,1,2) = s(:,3,1,2,1)
          s(:,1,3,2,1) = s(:,3,1,1,2)
          s(:,1,3,2,2) = s(:,3,1,2,2)
          s(:,1,3,4,1) = s(:,3,1,4,1)
	  s(:,1,3,4,2) = s(:,3,1,4,2)
          s(:,1,3,1,4) = s(:,3,1,1,4)
	  s(:,1,3,2,4) = s(:,3,1,2,4)
	  s(:,1,3,3,1) = s(:,3,1,3,1)
	  s(:,1,3,3,2) = s(:,3,1,3,2)
	  s(:,1,3,1,3) = s(:,3,1,1,3)
	  s(:,1,3,2,3) = s(:,3,1,2,3)
          Nparam(1,3) = Nparam(3,1)

!         C-U
          dREF(1:1,3,4) = (/ 4.81 /)             ! WWc
          alpam(1:1,3,4) = (/ 2.02 /)
          alpbm(1:1,3,4) = (/-2.50 /)            !Br3 CHECK when particles are inverted
          s(1:1,3,4,1,1) = (/2*wc/ )               !Br2  1,1 -> q1=0, q2=0     
          s(1:1,3,4,1,2) = (/ z /)                 !Br2  1,2 -> q1=0, q2=1
          s(1:1,3,4,2,1) = (/ z /)                 !Br2  2,1 -> q1=1, q2=0
          s(1:1,3,4,2,2) = (/ z /)                 !Br2  2,2 -> q1=1, q2=1
          s(1:1,3,4,4,1) = (/ z /)                 !4->DMS/CMCT
	  s(1:1,3,4,4,2) = (/ z /)                 !4->DMS/CMCT
          s(1:1,3,4,1,4) = (/ z /) 		   !4->DMS/CMCT
	  s(1:1,3,4,2,4) = (/ z /) 		   !4->DMS/CMCT
	  s(1:1,3,4,3,1) = (/ 2*wc*sha /)          !3->SHAPE
	  s(1:1,3,4,3,2) = (/ 2*wc*sha /)          !3->SHAPE
	  s(1:1,3,4,1,3) = (/ 2*wc*sha /)          !3->SHAPE
	  s(1:1,3,4,2,3) = (/ 2*wc*sha /)          !3->SHAPE
          Nparam(3,4) = 1
          
          dREF(:,4,3) = dREF(:,3,4)
          alpam(:,4,3) = alpbm(:,3,4)
          alpbm(:,4,3) = alpam(:,3,4)
          s(:,4,3,1,1) = s(:,3,4,1,1)
          s(:,4,3,1,2) = s(:,3,4,2,1)
          s(:,4,3,2,1) = s(:,3,4,1,2)
          s(:,4,3,2,2) = s(:,3,4,2,2)
          s(:,4,3,4,1) = s(:,3,4,4,1)
	  s(:,4,3,4,2) = s(:,3,4,4,2)
          s(:,4,3,1,4) = s(:,3,4,1,4)
	  s(:,4,3,2,4) = s(:,3,4,2,4)
	  s(:,4,3,3,1) = s(:,3,4,3,1)
	  s(:,4,3,3,2) = s(:,3,4,3,2)
	  s(:,4,3,1,3) = s(:,3,4,1,3)
	  s(:,4,3,2,3) = s(:,3,4,2,3)
          Nparam(4,3) = Nparam(3,4)

!         G-G
          dREF(1:3,1,1) =  (/ 6.00, 6.00, 6.73 /)        !     HWc, HWt, SSt   6.22 -> 6.0, 6.25 -> 6.00
          alpam(1:3,1,1) = (/ 1.27, 3.02, -1.78 /)
          alpbm(1:3,1,1) = (/ 2.90, 0.82, -1.83 /)
          s(1:3,1,1,1,1) = (/ 2*noWc, 2*noWc, 2*noWc /)   !Br2  1,1 -> q1=0, q2=0     
          s(1:3,1,1,1,2) = (/ z, z, 2*noWcq /)             !Br2  1,2 -> q1=0, q2=1
          s(1:3,1,1,2,1) = (/  2*noWcq, 2*noWcq, 2*noWcq/)   !Br2  2,1 -> q1=1, q2=0
          s(1:3,1,1,2,2) = (/ z, z, 2*noWcq/)             !Br2  2,2 -> q1=1, q2=1
          s(1:3,1,1,4,1) = (/ 2*noWc, 2*noWc, 2*noWc /)   !4->DMS/CMCT
	  s(1:3,1,1,4,2) = (/ 2*noWc, 2*noWc, 2*noWc /)   !4->DMS/CMCT
          s(1:3,1,1,1,4) = (/ 2*noWc, 2*noWc, 2*noWc /)   !4->DMS/CMCT
	  s(1:3,1,1,2,4) = (/ 2*noWc, 2*noWc, 2*noWc /)   !4->DMS/CMCT
          s(1:3,1,1,3,1) = (/ 2*noWc, 2*noWc, 2*noWc /)   !3->SHAPE
	  s(1:3,1,1,3,2) = (/ 2*noWc, 2*noWc, 2*noWc /)   !3->SHAPE
          s(1:3,1,1,1,3) = (/ 2*noWc, 2*noWc, 2*noWc /)   !3->SHAPE
	  s(1:3,1,1,2,3) = (/ 2*noWc, 2*noWc, 2*noWc /)   !3->SHAPE
          Nparam(1,1) = 3

!         G-U
          dREF(1:1,1,4) = (/ 5.05 /)               ! WWc         
          alpam(1:1,1,4) = (/ 2.29 /)
          alpbm(1:1,1,4) = (/ 1.68 /)
          s(1:1,1,4,1,1) = (/ 2.1*wcCanonic /)            !Br2  1,1 -> q1=0, q2=0 !OKKIO : cambiato a mano 14.7 -> 16!!!
          s(1:1,1,4,1,2) = (/ z /)                 !Br2  1,2 -> q1=0, q2=1
          s(1:1,1,4,2,1) = (/ z /)                 !Br2  2,1 -> q1=1, q2=0
          s(1:1,1,4,2,2) = (/ z /)                 !Br2  2,2 -> q1=1, q2=1
          s(1:1,1,4,4,1) = (/ z /)  		   !4->DMS/CMCT
	  s(1:1,1,4,4,2) = (/ z /)  		   !4->DMS/CMCT
          s(1:1,1,4,1,4) = (/ z /)                 !4->DMS/CMCT
	  s(1:1,1,4,2,4) = (/ z /)                 !4->DMS/CMCT
	  s(1:1,1,4,3,1) = (/ 2.1*wc*sha /)        !3->SHAPE
	  s(1:1,1,4,3,2) = (/ 2.1*wc*sha /)        !3->SHAPE
	  s(1:1,1,4,1,3) = (/ 2.1*wc*sha /)        !3->SHAPE
	  s(1:1,1,4,2,3) = (/ 2.1*wc*sha /)        !3->SHAPE
          Nparam(1,4) = 1
          
          dREF(:,4,1) = dREF(:,1,4)
          alpam(:,4,1) = alpbm(:,1,4)
          alpbm(:,4,1) = alpam(:,1,4)
          s(:,4,1,1,1) = s(:,1,4,1,1)
          s(:,4,1,1,2) = s(:,1,4,2,1)
          s(:,4,1,2,1) = s(:,1,4,1,2)
          s(:,4,1,2,2) = s(:,1,4,2,2)
          s(:,4,1,4,1) = s(:,1,4,4,1)
	  s(:,4,1,4,2) = s(:,1,4,4,2)
          s(:,4,1,1,4) = s(:,1,4,1,4)
	  s(:,4,1,2,4) = s(:,1,4,2,4)
	  s(:,4,1,3,1) = s(:,1,4,3,1)
	  s(:,4,1,3,2) = s(:,1,4,3,2)
	  s(:,4,1,1,3) = s(:,1,4,1,3)
	  s(:,4,1,2,3) = s(:,1,4,2,3)
          Nparam(4,1) = Nparam(1,4)

!         U-U
          dREF(1:3,4,4) =  (/ 4.94, 4.84, 5.63 /)           ! WWc, WWt, wht
          alpam(1:3,4,4) = (/ 1.85, -1.88, 2.36 /)
          alpbm(1:3,4,4) = (/ 2.48, 1.71, -2.57 /)
          s(1:3,4,4,1,1) = (/ 2*wc, 2*wc, 2*noWc /)        !Br2  1,1 -> q1=0, q2=0
          s(1:3,4,4,1,2) = (/ z, z, 2*noWcq /)                !Br2  1,2 -> q1=0, q2=1
          s(1:3,4,4,2,1) = (/ z, z, z /)                    !Br2  2,1 -> q1=1, q2=0
          s(1:3,4,4,2,2) = (/ z, z, z /)                    !Br2  2,2 -> q1=1, q2=1
          s(1:3,4,4,4,1) = (/ z, z, 2*noWc /)		    !4->DMS/CMCT
	  s(1:3,4,4,4,2) = (/ z, z, 2*noWc /)		    !4->DMS/CMCT
          s(1:3,4,4,1,4) = (/ z, z, 2*noWc /)		    !4->DMS/CMCT
	  s(1:3,4,4,2,4) = (/ z, z, 2*noWc /)		    !4->DMS/CMCT
	  s(1:3,4,4,3,1) = (/ 2*wc*sha, 2*wc*sha, 2*noWc /) !3->SHAPE
	  s(1:3,4,4,3,2) = (/ 2*wc*sha, 2*wc*sha, 2*noWc /) !3->SHAPE
	  s(1:3,4,4,1,3) = (/ 2*wc*sha, 2*wc*sha, 2*noWc /) !3->SHAPE
	  s(1:3,4,4,2,3) = (/ 2*wc*sha, 2*wc*sha, 2*noWc /) !3->SHAPE
          Nparam(4,4) = 3

 
!     ------------------------------------------------------------------
!         G: CY G1 G2
          planarityDistEq(1,1:3) = (/ 0, 0, 0/)
!         A: CY A1 A2
          planarityDistEq(2,1:3) = (/ 0, 0, 0/)
!         C: CA CY C1
          planarityDistEq(3,1:3) = (/ 0, 0, 0/)
!         U: CA CY U1
          planarityDistEq(4,1:3) = (/ 0, 0, 0/)

      end subroutine fillRNAHBparams

C     --------------------------------------------------------------

      subroutine fillDNAHBparams()
          use RNAHBparams
          use score

          implicit none
          double precision z
          
          wc = score_RNA(18)
          wcCanonic = score_RNA(19)
          noWc = score_RNA(20)


                    
          
!         A-A
          dREF(1:3,2,2) =  (/ 5.63, 6.84, 5.92/)       !WWt, HHt, HWt                    
          alpam(1:3,2,2) = (/ 2.40, 1.05, 2.72/)
          alpbm(1:3,2,2) = (/ 2.40, 1.08, 1.29/)
          s(1:3,2,2,1,1) = (/2*wc, 2*noWc, 2*noWc/)     !Br2  1,1 -> q1=0, q2=0
          s(1:3,2,2,1,2) = (/    z,2*noWcq, z/)            !Br2  1,2 -> q1=0, q2=1
          s(1:3,2,2,2,1) = (/    z,2*noWcq, 2*noWcq/)       !Br2  2,1 -> q1=1, q2=0
          s(1:3,2,2,2,2) = (/    z,2*noWcq, z/)            !Br2  2,2 -> q1=1, q2=1
          Nparam(2,2) = 3                !Br3  CHECK angle sign when order is inverted when one is negative

!         A-C
          dREF(1:3,2,3) = (/7.26, 5.78, 4.78/)         !wsc, HWt, +Wc  
          alpam(1:3,2,3) = (/2.36, 1.36, 2.64/)
          alpbm(1:3,2,3) = (/0.75, 2.30, 1.78/)   
          s(1:3,2,3,1,1) = (/ 2*noWc,  2*noWc, z /) !Br2  1,1 -> q1=0, q2=0
          s(1:3,2,3,1,2) = (/ 2*noWcq, z, z/)                 !Br2  1,2 -> q1=0, q2=1
          s(1:3,2,3,2,1) = (/ z, 2*noWcq, 2.2*tit/)            !Br2  2,1 -> q1=1, q2=0
          s(1:3,2,3,2,2) = (/ z, z, z/)                      !Br2  2,2 -> q1=1, q2=1

          Nparam(2,3) = 3
          dREF(:,3,2) = dREF(:,2,3)
          alpam(:,3,2) = alpbm(:,2,3)
          alpbm(:,3,2) = alpam(:,2,3)
          s(:,3,2,1,1) = s(:,2,3,1,1)
          s(:,3,2,1,2) = s(:,2,3,2,1)
          s(:,3,2,2,1) = s(:,2,3,1,2)
          s(:,3,2,2,2) = s(:,2,3,2,2)
          Nparam(3,2) = Nparam(2,3)

!         A-G
!          dREF(1:4,2,1) = (/ 4.88, 6.63, 6.17, 6.05/)         !  WW_c, HSt, sst, +Hc  
!          alpam(1:4,2,1) = (/ 3.04, 1.19, -2.01, 2.62/)
!          alpbm(1:4,2,1) = (/ 2.58, -1.69, -1.57, 0.89/)
!          s(1:4,2,1,1,1) = (/ 2*wc, 2*noWc, 2*noWc, z/)        !Br2  1,1 -> q1=0, q2=0
!          s(1:4,2,1,1,2) = (/ z, 2*noWcq, 2*noWcq, z/)            !Br2  1,2 -> q1=0, q2=1
!          s(1:4,2,1,2,1) = (/ z, 2*noWcq, 2*noWcq, 2*tit/)      !Br2  2,1 -> q1=1, q2=0
!          s(1:4,2,1,2,2) = (/ z, 2*noWcq, 2*noWcq, 2*tit/)      !Br2  2,2 -> q1=1, q2=1

!          Nparam(2,1) = 4
!          dREF(:,1,2) = dREF(:,2,1)
!          alpam(:,1,2) = alpbm(:,2,1)
!          alpbm(:,1,2) = alpam(:,2,1)
!          s(:,1,2,1,1) = s(:,2,1,1,1)
!          s(:,1,2,1,2) = s(:,2,1,2,1)
!          s(:,1,2,2,1) = s(:,2,1,1,2)
!          s(:,1,2,2,2) = s(:,2,1,2,2)
          
!          Nparam(1,2) = Nparam(2,1)


!         A-G
          dREF(1:3,2,1) = (/ 4.88, 6.63, 6.05/)         !  WW_c, HSt, +Hc
          alpam(1:3,2,1) = (/ 3.04, 1.19,  2.62/)
          alpbm(1:3,2,1) = (/ 2.58, -1.69,  0.89/)
          s(1:3,2,1,1,1) = (/ 2*wc, 2*noWc, z/)        !Br2  1,1 -> q1=0, q2=0
          s(1:3,2,1,1,2) = (/ z, 2*noWcq, z/)            !Br2  1,2 -> q1=0, q2=1
          s(1:3,2,1,2,1) = (/ z, 2*noWcq, 1.6*tit/)     !Br2  2,1 -> q1=1, q2=0
          s(1:3,2,1,2,2) = (/ z, 2*noWcq, 1.6*tit/)      !Br2  2,2 -> q1=1, q2=1

          Nparam(2,1) = 3
          dREF(:,1,2) = dREF(:,2,1)
          alpam(:,1,2) = alpbm(:,2,1)
          alpbm(:,1,2) = alpam(:,2,1)
          s(:,1,2,1,1) = s(:,2,1,1,1)
          s(:,1,2,1,2) = s(:,2,1,2,1)
          s(:,1,2,2,1) = s(:,2,1,1,2)
          s(:,1,2,2,2) = s(:,2,1,2,2)

          Nparam(1,2) = Nparam(2,1)

!         A-U
          dREF(1:3,2,4) = (/ 4.92, 5.78, 5.89/)        ! WWc, HWt, HWc
          alpam(1:3,2,4) = (/ 2.84, 0.91, 0.93/)
          alpbm(1:3,2,4) = (/ 2.36, 1.82, 2.46/)
          s(1:3,2,4,1,1) = (/ 2.2*wcCanonic, 2*noWc, 2*noWc/ )     !Br2  1,1 -> q1=0, q2=0      !OKKIO : cambiato a mano 14.4 -> 16!!!
          s(1:3,2,4,1,2) = (/ z, z, z/)                     !Br2  1,2 -> q1=0, q2=1
          s(1:3,2,4,2,1) = (/ z, 2*noWcq, 2*noWcq/)          !Br2  2,1 -> q1=1, q2=0
          s(1:3,2,4,2,2) = (/ z, z, z/)                     !Br2  2,2 -> q1=1, q2=1
          Nparam(2,4) = 3
          
          dREF(:,4,2) = dREF(:,2,4)
          alpam(:,4,2) = alpbm(:,2,4)
          alpbm(:,4,2) = alpam(:,2,4)
          s(:,4,2,1,1) = s(:,2,4,1,1)
          s(:,4,2,1,2) = s(:,2,4,2,1)
          s(:,4,2,2,1) = s(:,2,4,1,2)
          s(:,4,2,2,2) = s(:,2,4,2,2)
          Nparam(4,2) = Nparam(2,4)

! !         C-C
!           dREF(1:1,3,3) =  (/ 4.91/)  !Br3 NO significant CC pairing
!           alpam(1:1,3,3) = (/ 2.22/)
!           alpbm(1:1,3,3) = (/ 2.24/)
!           s(1:1,3,3,1,1) = (/ z /) !Br2  1,1 -> q1=0, q2=0
!           s(1:1,3,3,1,2) = (/ z /) !Br2  1,2 -> q1=0, q2=1
!           s(1:1,3,3,2,1) = (/ z /) !Br2  2,1 -> q1=1, q2=0
!           s(1:1,3,3,2,2) = (/ z /) !Br2  2,2 -> q1=1, q2=1
!           Nparam(3,3) = 1

!         C-G
          dREF(1:3,3,1) = (/ 4.75, 5.28, 5.68/)              ! WWc, WWt, +Wc
          alpam(1:3,3,1) = (/ 2.17, 1.64, 1.92/)
          alpbm(1:3,3,1) = (/ 2.71, -3.07, 2.34/)
          s(1:3,3,1,1,1) = (/ 2.6*wcCanonic, 2*wc, z/)       !Br2  1,1 -> q1=0, q2=0
          s(1:3,3,1,1,2) = (/ z, z, z/)                      !Br2  1,2 -> q1=0, q2=1
          s(1:3,3,1,2,1) = (/ z, z, 0.5*tit/)                  !Br2  2,1 -> q1=1, q2=0
          s(1:3,3,1,2,2) = (/ z, z, z/)                      !Br2  2,2 -> q1=1, q2=1
          Nparam(3,1) = 3
          
          dREF(:,1,3) = dREF(:,3,1)
          alpam(:,1,3) = alpbm(:,3,1)
          alpbm(:,1,3) = alpam(:,3,1)
          s(:,1,3,1,1) = s(:,3,1,1,1)
          s(:,1,3,1,2) = s(:,3,1,2,1)
          s(:,1,3,2,1) = s(:,3,1,1,2)
          s(:,1,3,2,2) = s(:,3,1,2,2)
          Nparam(1,3) = Nparam(3,1)

!         C-U
          dREF(1:1,3,4) = (/ 4.81 /)             ! WWc
          alpam(1:1,3,4) = (/ 2.02 /)
          alpbm(1:1,3,4) = (/-2.50 /)            !Br3 CHECK when particles are inverted
          s(1:1,3,4,1,1) = (/2*wc/ )               !Br2  1,1 -> q1=0, q2=0     
          s(1:1,3,4,1,2) = (/ z /)                 !Br2  1,2 -> q1=0, q2=1
          s(1:1,3,4,2,1) = (/ z /)                 !Br2  2,1 -> q1=1, q2=0
          s(1:1,3,4,2,2) = (/ z /)                 !Br2  2,2 -> q1=1, q2=1
          Nparam(3,4) = 1
          
          dREF(:,4,3) = dREF(:,3,4)
          alpam(:,4,3) = alpbm(:,3,4)
          alpbm(:,4,3) = alpam(:,3,4)
          s(:,4,3,1,1) = s(:,3,4,1,1)
          s(:,4,3,1,2) = s(:,3,4,2,1)
          s(:,4,3,2,1) = s(:,3,4,1,2)
          s(:,4,3,2,2) = s(:,3,4,2,2)
          Nparam(4,3) = Nparam(3,4)

!         G-G
          dREF(1:3,1,1) =  (/ 6.22, 6.25, 6.73 /)        !     HWc, HWt, SSt
          alpam(1:3,1,1) = (/ 1.27, 3.02, -1.78 /)
          alpbm(1:3,1,1) = (/ 2.90, 0.82, -1.83 /)
          s(1:3,1,1,1,1) = (/ 2*noWc, 2*noWc, 2*noWc /)   !Br2  1,1 -> q1=0, q2=0     
          s(1:3,1,1,1,2) = (/ z, z, 2*noWcq /)             !Br2  1,2 -> q1=0, q2=1
          s(1:3,1,1,2,1) = (/  2*noWcq, 2*noWcq, 2*noWcq/)   !Br2  2,1 -> q1=1, q2=0
          s(1:3,1,1,2,2) = (/ z, z, 2*noWcq/)             !Br2  2,2 -> q1=1, q2=1
          Nparam(1,1) = 3

!         G-U
          dREF(1:1,1,4) = (/ 5.05 /)               ! WWc         
          alpam(1:1,1,4) = (/ 2.29 /)
          alpbm(1:1,1,4) = (/ 1.68 /)
          s(1:1,1,4,1,1) = (/ 2.1*wc /)            !Br2  1,1 -> q1=0, q2=0 !OKKIO : cambiato a mano 14.7 -> 16!!!
          s(1:1,1,4,1,2) = (/ z /)                 !Br2  1,2 -> q1=0, q2=1
          s(1:1,1,4,2,1) = (/ z /)                 !Br2  2,1 -> q1=1, q2=0
          s(1:1,1,4,2,2) = (/ z /)                 !Br2  2,2 -> q1=1, q2=1
          Nparam(1,4) = 1
          
          dREF(:,4,1) = dREF(:,1,4)
          alpam(:,4,1) = alpbm(:,1,4)
          alpbm(:,4,1) = alpam(:,1,4)
          s(:,4,1,1,1) = s(:,1,4,1,1)
          s(:,4,1,1,2) = s(:,1,4,2,1)
          s(:,4,1,2,1) = s(:,1,4,1,2)
          s(:,4,1,2,2) = s(:,1,4,2,2)
          Nparam(4,1) = Nparam(1,4)

!         U-U
          dREF(1:3,4,4) =  (/ 4.94, 4.84, 5.63 /)           ! WWc, WWt, wht
          alpam(1:3,4,4) = (/ 1.85, -1.88, 2.36 /)
          alpbm(1:3,4,4) = (/ 2.48, 1.71, -2.57 /)
          s(1:3,4,4,1,1) = (/ 2*wc, 2*wc, 2*noWc /)        !Br2  1,1 -> q1=0, q2=0
          s(1:3,4,4,1,2) = (/ z, z, 2*noWcq /)                !Br2  1,2 -> q1=0, q2=1
          s(1:3,4,4,2,1) = (/ z, z, z /)                    !Br2  2,1 -> q1=1, q2=0
          s(1:3,4,4,2,2) = (/ z, z, z /)                    !Br2  2,2 -> q1=1, q2=1
          Nparam(4,4) = 3

 
!     ------------------------------------------------------------------
!         G: CY G1 G2
!          planarityDistEq(1,1:3) = (/ 0, 0, 0/)
!         A: CY A1 A2
!          planarityDistEq(2,1:3) = (/ 0, 0, 0/)
!         C: CA CY C1
!          planarityDistEq(3,1:3) = (/ 0, 0, 0/)
!!         U: CA CY U1
!          planarityDistEq(4,1:3) = (/ 0, 0, 0/)
          planarityDistEq(1,1:3) = (/ 5.2, 3.7, 2.7/)
!         A
          planarityDistEq(2,1:3) = (/ 4.5, 4.0, 2.8/)
!         C
          planarityDistEq(3,1:3) = (/ 1.9, 0.79, 0.38/)
!         U
          planarityDistEq(4,1:3) = (/ 2.7, 1.4, 0.30/)


      end subroutine fillDNAHBparams

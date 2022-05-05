! Subroutine mdmin  !! Identical to Genevieve version Dec04
!
! Minimizes the energy. It first uses a steepest descent algorithm and after 
! several hundred steps then uses a conjugate gradient method.

! This minimization is done with only a minimal knowledge of the physics
! of the problem so that it is portable
!
! The minimization uses a simple steepest descent with variable step size.
!
! This file contains 3 routines
!

MODULE minimization
! This module defines a number of parameters used during the minimization
  implicit none
  save
  
  integer :: MAXSTEP_MIN
  real(8) :: MINF 
  real(8) :: STEP_DMD
  real(8), parameter :: FRIC = 0.5  !! 
  real(8), parameter :: reps = 1.0e-12 
  real(8)  :: MINF2 
 
 contains

    subroutine calcforce(scale,xa,fa,etot, time)
    use defs
    real(8) :: scale, time
    real(8), intent(out) :: etot
    real(8), dimension(vecsize), intent(inout) :: xa
    real(8), dimension(vecsize), intent(out) :: fa

    real(8) :: E_prot, E_RNA, E_cross
    E_prot  = 0
    E_RNA   = 0
    E_cross = 0

    ! Call the force field
    if (prot_simulation) then
      call calcforce_protein(scale,xa(1:N_prot*3),fa(1:N_prot*3),E_prot)
    endif
    if (RNA_simulation) then
      call calcforce_RNA(scale,xa(N_prot*3+1:Natoms*3),fa(N_prot*3+1:Natoms*3),E_RNA,.true., time)
    endif
    if (prot_simulation .and. RNA_simulation) then
      call calcforce_RNA_protein(scale, xa, fa, E_cross)
    endif
    etot = E_prot+E_RNA+E_cross
 
  end subroutine calcforce
  

END MODULE minimization

MODULE minimization1

  use defs
  implicit none

  real(8), parameter :: acc=1.d-6 !lal changed from -6 to -12
  real(8) :: eref

END MODULE minimization1

subroutine finilise_minimization()
  use defs
  use int_defs
  use constraints
  
  implicit none

  if (.not. MINIMIZATION_TYPE=='MIN_damped') then

    ! Deallocate variables used in min_converge1
    deallocate(gra)
    deallocate(scl)
    if (use_ics) then
        deallocate(gras)
        deallocate(vars)
        deallocate(ovrs)
        deallocate(poss)
    endif
  endif

  if (use_ics) then

    ! Deallocate variables used in bodies
    deallocate(lc)
    deallocate(kch)
    deallocate(nbbvar)
    deallocate(nbasepp)
    deallocate(ipiv)
    deallocate(iseg1)
    deallocate(iseg2)
    deallocate(ibseg)
    deallocate(log3t)
    deallocate(log5t)

    ! Deallocate variables used in varis
    deallocate(var)
    deallocate(ovr)
    deallocate(svars)

  endif

  if (constrained_fragments.or.restrained_fragments) then

    ! deallocate variables used in set_harmonic_constraints
    deallocate(constrained_atoms)
    deallocate(fc)
    deallocate(pc)

  endif 

end subroutine finilise_minimization

subroutine initialise_minimization()
  use defs
  use int_defs
  use constraints
  use system_defs_RNA
  use minimization
  use geometric_corrections

  implicit none
  integer :: ierror
  character(len=20) dummy


  if (use_ics) then 
    ! Allocate variables used in bodies
    if (.not. allocated(lc)) allocate(lc(nfrag))
    if (.not. allocated(kch)) allocate(kch(nfrag))
    if (.not. allocated(nbbvar)) allocate(nbbvar(nfrag))
    if (.not. allocated(nbasepp)) allocate(nbasepp(nres))
    if (.not. allocated(ipiv)) allocate(ipiv(nfrag))
    if (.not. allocated(iseg1)) allocate(iseg1(nfrag*2))
    if (.not. allocated(iseg2)) allocate(iseg2(nfrag*2))
    if (.not. allocated(ibseg)) allocate(ibseg(nfrag*2))
    if (.not. allocated(log3t)) allocate(log3t(nfrag*2))
    if (.not. allocated(log5t)) allocate(log5t(nfrag*2))
  endif

  call GET_ENVIRONMENT_VARIABLE('Max_Num_Iter_Minimization', dummy)
  if (dummy .eq. '') then
     MAXSTEP_MIN = 45000
  else  
     read(dummy,*) MAXSTEP_MIN
  endif
  
  call GET_ENVIRONMENT_VARIABLE('Force_Threshold_Minimization', dummy)
  if (dummy .eq. '') then
     MINF = 0.01
  else  
     read(dummy,*) MINF
  endif
  MINF2 = MINF*MINF  

  if(MINIMIZATION_TYPE=='MIN_damped')then
     call GET_ENVIRONMENT_VARIABLE('Time_Step_DampedMD', dummy)
     if (dummy .eq. '') then
        STEP_DMD = 0.1
     else  
        read(dummy,*) STEP_DMD
     endif
     
     if(taskid .eq. 0 ) then
        open(unit=FLOG,file=MASTER_LOGFILE,status='unknown',action='write',position='append',iostat=ierror) 
        write(flog,*) ' '
        write(flog,*) '****** Damped MD Minimization  parameters **************************************'
        write(flog,'(1X,A51,I8)')    ' Maximum number of iteration (minimization)      : ', MAXSTEP_MIN
        write(flog,'(1X,A51,F8.5)')    ' Force threshold for minimization                : ', MINF
        write(flog,'(1X,A51,F8.5)')    ' Time step for the damped MD minimization        : ', STEP_DMD
        write(flog,*) ' '
        close(flog)
     endif
  else
     if(taskid .eq. 0 ) then
        open(unit=FLOG,file=MASTER_LOGFILE,status='unknown',action='write',position='append',iostat=ierror) 
        write(flog,*) ' '
        write(flog,*) '****** Steepest Descent Minimization  parameters **************************************'
        write(flog,'(1X,A51,I8)')    ' Maximum number of iteration (minimization)      : ', MAXSTEP_MIN
        write(flog,'(1X,A51,F8.5)')    ' Force threshold for minimization                : ', MINF
        write(flog,*) ' '
        close(flog)
     endif
  endif
end subroutine initialise_minimization


subroutine min_converge(success, logfile, relax_name)
  use defs
  use md_defs, only:calc_SAXS_force, saxs_serial_step_min, saxs_serial_step, calc_SAXS_modul
  use minimization
  use geometric_corrections
  use writetofile
  
  implicit none
  
  character(len=20), intent(in) :: LOGFILE
  logical, intent(out) :: success
  integer :: iter, i, npart, ierror
  real(8) :: ftot, ftot2, VERSTEP, delr, rmsd
  real(8) :: dotvf, v2, frac, ref_rmsd, dif
  real(8) :: enerexit, scale
  real(8), dimension(VECSIZE) :: vel

  character(len=100),intent(in),optional :: relax_name
  character(len=100) :: fname
  if(present(relax_name)) then
    fname = trim(relax_name)
  else
    fname = 'relaxation.pdb'
  endif

!  type (t_conformations) :: conformation

!  scale = conformation%energyscale

!  write(*,*) "scale---mdmmin", scale
 
  scale = 1.0
 
  open(unit=FLOG,file=LOGFILE,status='unknown',action='write',position='append',iostat=ierror)

  enerexit=0
  evalf_number = 0   ! Set the number of force evaluations to zero
  
  
  calc_SAXS_force = .true.
  calc_SAXS_modul = .false.
  call calcforce(scale,pos,force,total_energy, 0.0d0)
  evalf_number = evalf_number + 1

  ftot2 = dot_product(force,force)
  
  vel = 0.0  ! Vectorial operation
  VERSTEP = STEP_DMD
  do iter=1, MAXSTEP_MIN   
     
     call minimize_rmsd(natoms,posref, pos, delr, rmsd, npart)
     
     ref_rmsd = rmsd
     ! v(t+ 1/2 * dt) = v(t)+(dt/2)*F(t)/m   
     do i=1, vecsize
        vel(i)= vel(i) + (0.5*VERSTEP/MASS(i)) * force(i)
     enddo
     
     ! x(t+dt) = x(t) + dt * v(t+ 1/2* dt)     
     v2 = 0.0
     do i=1, vecsize
        v2 = v2 + vel(i) * vel(i) 
     enddo
     frac = exp(-VERSTEP * FRIC * sqrt(v2))
!     frac = exp(-VERSTEP * FRIC * 10)
     
     if (frac .le. 0.0000001) then
        frac = 0.0001 
     endif
     
!     do i=1, vecsize
!        vel(i) = vel(i) * frac
!     enddo
     do i=1, vecsize
!        pos(i) = pos(i) + tanh(VERSTEP * vel(i)/(1.0 + sqrt(ftot2)))
        ! way more robust in case of clashes, and generally more efficient
        pos(i) = pos(i) + tanh(vel(i))
     enddo
     
     ! Calculating F(t+dt)
    calc_SAXS_force = .false.
    if (saxs_serial_step .gt. 0)then
        calc_SAXS_force = (mod(iter, saxs_serial_step_min) .eq. 0)
        calc_SAXS_modul = .false.
    endif 
     
     call calcforce(scale,pos,force,total_energy,0.0d0)
!     write(flog,*) 'total energy for F(t+dt)', total_energy
     evalf_number = evalf_number + 1
     ftot2 = 0.0
     do i=1, vecsize
        ftot2 = ftot2 + force(i) * force(i)
!        if (force(i) > 10) then
!      write(98,*) i, force(i)
!        endif
     enddo
     
     call minimize_rmsd(natoms, posref, pos, delr, rmsd, npart)
     
     dif = abs(rmsd - ref_rmsd)
!     if (dif < reps ) then
!     endif
     
     if (mod(iter, 50) == 0) then
       if (debug_status .eq. 'debug') write(*, "(' ','iter: ',i5,' energy: ',f12.4,' frac: ', f8.4,       &
             & ' ftot: ', f12.4, ' dotvf: ', f12.6, ' vel: ', f12.4,  '  delr: ',  &       
             & f7.4, ' rmsd: ',f7.3,' npart: ',i4,' verstep: ', f6.3)")  &
             & iter, total_energy, frac, sqrt(ftot2), dotvf, sqrt(v2),   &
             & delr, rmsd, npart, VERSTEP 
        write(flog, "(' ','iter: ',i5,' energy: ',f12.4,' frac: ', f8.4,       &
             & ' ftot: ', f12.4, ' dotvf: ', f12.6, ' vel: ', f12.4,  '  delr: ',  &       
             & f7.4, ' rmsd: ',f7.3,' npart: ',i4,' verstep: ', f6.3)")  &
             & iter, total_energy, frac, sqrt(ftot2), dotvf, sqrt(v2),   &
             & delr, rmsd, npart, VERSTEP 
     endif
     
     
     if ( (ftot2 < MINF2 .or. dif < reps) .and. (iter .gt. 5) ) then
        write(flog,*) sqrt(ftot2), sqrt(MINF2), dif, rmsd, ref_rmsd,reps, iter
        exit
     endif
     
     ! let velocities to be zero if (v dot F) is negative
     dotvf = 0.0
     do i=1, vecsize
        dotvf = dotvf + vel(i) * force(i)
     enddo
     if (dotvf < 0.0) then
        do i=1, vecsize
           vel(i) = 0.0
        enddo
        VERSTEP = 0.95 * VERSTEP
     endif
     
     if (VERSTEP < 0.01) VERSTEP = 1.8 * VERSTEP
     
     ! v(t+dt) = v(t+1/2*dt) + dt/2 * F(t+dt)/m 
     do i=1, vecsize  
        vel(i) = vel(i) + (0.5 * VERSTEP / MASS(i)) * force(i)
     enddo
     if (mod(iter,50) .eq. 0) then
       call write_to_file('relaxe',0,natoms,pos,posref,ndigits_filename,fname,iter,0.0d0, & 
         total_energy,.true.,singlefile)
     endif
  enddo
  
  ftot = sqrt(ftot2)
  if (ftot < MINF  .and. iter > 1) then
     success = .true.
     if (debug_status .eq. 'debug') then 
        write(*,*) 'Minimization successful   ftot : ', ftot
        write(*,*) 'Minimization Energy : ', total_energy 
     endif
     write(FLOG,*) 'Minimization successful   ftot : ', ftot
     write(FLOG,*) 'Minimization Energy : ', total_energy 
  else
     success = .false.
     if (debug_status .eq. 'debug') then 
        write(*,*) 'Minimization failed   ftot : ', ftot
        write(*,*) 'Minimization Energy : ', total_energy
     endif
     write(FLOG,*) 'Minimization failed   ftot : ', ftot
     write(FLOG,*) 'Minimization Energy : ', total_energy
  endif
  
  close(flog)
end subroutine min_converge

!==================================
!
!   Steepest Descent Minimization
!
!===================================

subroutine min_converge1(success, logfile, relax_name)

  use defs
  use md_defs, only:calc_SAXS_force, saxs_serial_step, calc_SAXS_modul
  use minimization
  use minimization1
  use geometric_corrections
  use writetofile
  implicit none
 
  
  character(len=20), intent(in) :: LOGFILE
  logical, intent(out) :: success
  integer :: iter, i, npart, ierror
  real(8) :: ftot, ftot2, VERSTEP, delr, rmsd
  real(8) :: dotvf, v2, frac, ref_rmsd, dif
  real(8) :: enerexit, scale
  real(8), dimension(VECSIZE) :: vel
  real(8), parameter :: tfac=3.0d0
  character(len=100),intent(inout),optional :: relax_name
  character(len=100) :: fname
  integer :: natoms3
  integer :: nd,nw,nxa,nga,nxb,ngb,nfun
  real(8) :: delmax
  real(8), dimension(:), allocatable  :: w1
  
  if(present(relax_name)) then
    fname = trim(relax_name)
  else
    fname = 'relaxation.pdb'
  endif

  natoms3 = NATOMS*3

  if(present(relax_name)) then
    fname = trim(relax_name)
  else
    fname = 'relaxation.pdb'
  endif

  open(unit=FLOG,file=LOGFILE,status='unknown',action='write',position='append',iostat=ierror)

  enerexit=0
  evalf_number = 0   ! Set the number of force evaluations to zero
  
  !allocate(gra(natoms3))
  calc_SAXS_force = .true.
  calc_SAXS_modul = .false.

  call calcforce(scale,pos,force,total_energy,0.0d0)
  
  if(.not.use_ics)then 
     
     if (.not.allocated(var))allocate(var(natoms3))
     if (.not.allocated(scl))allocate(scl(natoms3))
     if (.not.allocated(gra))allocate(gra(natoms3))
     
     gra(1:natoms3)=-force(1:natoms3)
     scl(1:natoms3)=0.01d0
     var(1:natoms3)=pos(1:natoms3)
     nvar=natoms3
  else
     call varis

     ! Allocate stuff
     PRINT*,"ALLOCATED GRAS", allocated(gras)
     if (.not. allocated(gra)) allocate(gra(nvar))
     if (.not. allocated(scl)) allocate(scl(nvar))
     if (.not. allocated(gras)) allocate(gras(nvar))
     if (.not. allocated(vars)) allocate(vars(nvar))
     if (.not. allocated(ovrs)) allocate(ovrs(nvar))
     if (.not. allocated(poss)) allocate(poss(natoms3))
     PRINT*,"ALLOCATED GRAS", allocated(gras)

     call assemb!gradients
     
     scl(1:nvar)=tfac*0.1d0
 
  endif

  !print*,'gra',gra(1),force(1)
  evalf_number = evalf_number + 1

  
  delmax=0.
  nd=1+(nvar*(nvar+1))/2
  nw=nd+nvar
  nxa=nw+nvar
  nga=nxa+nvar
  nxb=nga+nvar
  ngb=nxb+nvar
  !print*,natoms3,nd,nw,nxa,nga,nxb,ngb

  !allocate(w(ngb))
  !allocate(w1(ngb))
  !w1(1:ngb)=0.d0

  !print*,'hello before MINFOR',scl(1),use_ics

  !allocate(w1(ngb))

  call minfor(nvar,var,total_energy,gra,scl,acc,MAXSTEP_MIN,evalf_number,fname,ngb)
 ! call minfor(nvar,var,total_energy,gra,scl,acc,w1,w1(1:nvar),w1(1:nvar),w1(1:nvar),w1(1:nvar),&
  !     & w1(1:nvar),w1(1:nvar),MAXSTEP_MIN,evalf_number,fname)

  if(debug_status.eq.'mini') call gradt
!!!stop
  !print*,'heree'
  if(evalf_number.lt.MAXSTEP_MIN) then
     success=.true.
     write(FLOG,401)evalf_number,nvar
401  format(/2x,'--- Convergence at step ',i4,' for ',i5,' vars ---'/)
  else
     success=.false.
     write(FLOG,402)evalf_number,nvar
402  format(/2x,'--- Step limit (',i4,') for ',i5,' vars ---'/)
     !deallocate(scl,var,gra) ! Moved to finilise_minimization
  endif
  

  close(flog)


  !print*,'heree finish'

!   stop
end subroutine min_converge1


subroutine minfor(n,x,f,g,scale1,acc,maxfun,nfun,fname,nnnn)
!subroutine minfor(n,x,f,g,scale1,acc,h,d,w,xa,ga,xb,gb,maxfun,nfun,fname) 

  use minimization1, only: nvar

  implicit double precision (a-h,o-z)
  implicit integer (i-n)
  real(8), intent(inout),dimension(nvar) ::  g,scale1
  real(8), dimension(nvar) ::  ga,x,xa,xb,gb,w,d
  real(8), dimension(nnnn+1):: h
  character(len=100),intent(inout) :: fname
  !dimension h(*) !,w(*),gb(*),xb(*)
  !print *,'in minfor ',g(1)                                        

  nfun=0 
  itr=0 
  np=n+1

!     set the hessian to a diagonal matrix depending on scale(.)        
  c=0.0

  do i=1,n
     c=max(c,abs(g(i)*scale1(i)))
     !!!print*,i,g(i),scale1(i),g(i)*scale1(i),c
  enddo
!     do 30 i=1,n 
!30   c=max(c,abs(g(i)*scale1(i))) 
  if (c.le.0.) c=1.d0
  k=(n*np)/2
 ! print*,'k',k
 

  h(1:k)=0.d0

!  do 40 i=1,k 
!40   h(i)=0.d0
!     print*,k,h(k)
  k=1
     
  do 50 i=1,n
  
     h(k)=0.01*c/scale1(i)**2
    
50   k=k+np-i
  !   print*,'k new',k
   !  stop
     
!     set some variables for the first iteration                        
  dff=0. 
110 fa=f 
  isfv=1

  xa(1:n)=x(1:n)
  ga(1:n)=g(1:n)

!  do 120 i=1,n 
!     xa(i)=x(i) 
!      print*,'xa',xa(i)                                                
!120  ga(i)=g(i)


     
!     begin the iteration by giving the required printing               
130  itr=itr+1 
! calculate the search direction of the iteration                   

!  do 150 i=1,n 
!150  d(i)=-ga(i)
  d(1:n)=-ga(1:n)
  
  call mc11e (h,n,d,w,n) 
!     calculate a lower bound on the step-length                        
!     and the initial directional derivative                            
  c=0. 
  dga=0. 
  do 160 i=1,n 
     c=max(c,abs(d(i)/scale1(i)))
     !!!!!print*,i,c,d(i),d(i)/scale1(i)
160  dga=dga+ga(i)*d(i) 
!     test if the search direction is downhill                          
  if (dga.ge.0.) go to 240 
!     set the initial step-length of the line search                    
  stmin=0. 
  stepbd=0. 
  steplb=acc/c 
  fmin=fa 
  gmin=dga 
  step=1.d0
  !!!print*,'c',c
  if (dff.le.0.) step=min(step,1./c) 
  if (dff.gt.0.) step=min(step,(dff+dff)/(-dga)) 
170 c=stmin+step

  !print*,c,stmin,step
  !     test whether func has been called maxfun times             
  if (nfun.ge.maxfun) go to 250 
  nfun=nfun+1 
  !     calculate another function value and gradient                     
  do 180 i=1,n 
180  xb(i)=xa(i)+c*d(i) 
    
  
  call move(nfun,xb,fb,gb,fname)
     
  !     store this function value if it is the smallest so far            
  isfv=min(2,isfv) 
  if (fb.gt.f) go to 220 
  if (fb.lt.f) go to 200 
  gl1=0. 
  gl2=0. 
  do 190 i=1,n 
     gl1=gl1+(scale1(i)*g(i))**2 
190  gl2=gl2+(scale1(i)*gb(i))**2
     
     if (gl2.ge.gl1) go to 220 
200  isfv=3 
     f=fb 
  do 210 i=1,n 
     x(i)=xb(i) 
210  g(i)=gb(i)

    
!     calculate the directional derivative at the new point             
220  dgb=0.
     
  do 230 i=1,n 
230  dgb=dgb+gb(i)*d(i) 
!     branch if we have found a new lower bound on the step-length      
  if (fb-fa.le.0.1*c*dga) go to 280 
  !     finish the iteration if the current step is steplb                
  if (step.gt.steplb) go to 270 
240 if (isfv.ge.2) go to 110 
!     at this stage the whole calculation is complete                   
250 if(nfun.lt.maxfun) then 
     nfun=nfun+1 
     !print *,'in minfor ',g(1)                                        
     call move(nfun,x,f,g,fname)
     
  endif
  return
  
!     calculate a new step-length by cubic interpolation                
270 stepbd=step 
  c=gmin+dgb-3.*(fb-fmin)/step 
  c=gmin/(c+gmin-sqrt(c*c-gmin*dgb)) 
  step=step*max(0.1d0,c) 
  go to 170 
  !     set the new bounds on the step-length                             
280 stepbd=stepbd-step 
  stmin=c 
  fmin=fb 
  gmin=dgb 
  !     calculate a new step-length by extrapolation                      
  step=9.*stmin 
  if (stepbd.gt.0.) step=0.5*stepbd 
  c=dga+3.*dgb-4.*(fb-fa)/stmin 
  if (c.gt.0.) step=min(step,stmin*max(1.d0,-dgb/c)) 
  if (dgb.lt.0.7*dga) go to 170 
  !     test for convergence of the iterations                            
  isfv=4-isfv 
  if (stmin+step.le.steplb) go to 240 
  !     revise the second derivative matrix                               
  ir=-n 
  do 290 i=1,n 
     xa(i)=xb(i) 
     xb(i)=ga(i) 
     d(i)=gb(i)-ga(i) 
290  ga(i)=gb(i) 
  call mc11a(h,n,xb,1./dga,w,ir,1,0.d0) 
  ir=-ir 
  call mc11a (h,n,d,1./(stmin*(dgb-dga)),d,ir,0,0.d0) 
  !     branch if the rank of the new matrix is deficient                 
  if (ir.lt.n) go to 250 
  !     begin another iteration                                           
  dff=fa-fb 
  fa=fb 
  go to 130 

end subroutine  minfor


subroutine move(icyc,vrc,esum,grc,fname1)

  use defs
  use minimization
  use minimization1
  use geometric_corrections
  use writetofile
  use md_defs, only:calc_SAXS_force, saxs_serial_step, saxs_serial_step_min, calc_SAXS_modul
  implicit double precision (a-h,o-z)
  implicit integer (i-n)
  
  integer, intent(in)::icyc
  real(8) :: esum, mintime
  real(8), intent(in),dimension(nvar) ::  vrc
  real(8), intent(out),dimension(nvar) ::  grc
  real(8) :: scale
  !character(len=100),intent(in),optional :: relax_name1
 
  character(len=100),intent(in) :: fname1
  !if(present(relax_name1)) then
  !  fname = trim(relax_name1)
  !else
  !  fname = 'relaxation.pdb'
  !endif

  var(1:nvar)=vrc(1:nvar)
  !print*,'move',icyc,var(nvar),vrc(nvar)
 
  !print *,'In move',vrc(5)
  !print *,'In move',vrc(6)
  if(icyc.eq.1) eref=total_energy
  !-------------------------------------------------------------------step
  if(.not.use_ics) then
     pos(1:nvar)=var(1:nvar)
  else
     call microb
  endif
  
  !call energy
  scale=1.0
  calc_SAXS_force = .false.
  if (saxs_serial_step .gt. 0)then
    calc_SAXS_force  = (mod(icyc, saxs_serial_step_min) .eq. 0)
    calc_SAXS_modul = .false.
  endif 
  
  mintime = icyc/1.0
  call calcforce(scale,pos,force,total_energy,mintime)

  if(.not.use_ics) then
     gra(1:nvar)=-force(1:nvar)
     !!!print*,'gra a',gra(1)
  else
     call assemb
  endif

  grc(1:nvar)=gra(1:nvar)
  

  !      print *,"move"
  !print*,'energy',eref,total_energy
  delta=total_energy-eref
  eref=total_energy
  esum=total_energy
!!!!!-----------------------------------------------------------max gradient
  gx=0.d0
  do k=1,nvar
      g=abs(gra(k))
      if(g.gt.gx) then
         gx=g
     endif
  enddo
!  !print*,'max',gx
!!!!------------------------------------------------------------------cycle
  !if(.not.quiet)then
  if(debug_status.eq.'mini')then
  write(6,20) icyc,total_energy,delta,gx
20 format(2x,i5,') ',f12.2,' D= ',e8.2,' G= ',e8.2)
endif
  !call minimize_rmsd(natoms, posref, pos, delr, rmsd, npart)
  rmsd=2.d0
  dif = abs(rmsd - ref_rmsd)
  
  call write_to_file('relaxe',0,natoms,pos,posref,ndigits_filename,fname1,nfun,0.0d0, & 
       total_energy,.true.,singlefile)


  ftot2 = 0.0
  do i=1, nvar
     ftot2 = ftot2 + force(i) * force(i)
     !        if (force(i) > 10) then
     !      write(98,*) i, force(i)
     !        endif
  enddo
  
  write(flog, "(' ','iter: ',i5,' energy: ',f12.4, &
             & ' ftot: ', f12.4, '  delr: ',  &       
             & f7.4, ' rmsd: ',f7.3,' npart: ',i4)")  &
             & icyc, total_energy, sqrt(ftot2),   &
             & delr, rmsd, npart 
  
end subroutine move


subroutine mc11a(a,n,z,sig,w,ir,mk,eps)

  use minimization1, only: nvar

  implicit double precision (a-h,o-z)
  implicit integer (i-n)
  
  real(8), intent(inout),dimension(nvar) ::  z

  dimension w(*),a(*)
      
  if(n.gt.1) goto 1
  a(1)=a(1)+sig *z(1)**2
  ir=1
  if(a(1).gt.0.)return
  a(1)=0.
  ir=0
  return
1 continue
  np=n+1
  if(sig.gt.0.)goto 40
  if(sig.eq.0..or.ir.eq.0)return
  ti=1./sig
  ij=1
  if(mk.eq.0)goto 10
  do 7 i=1,n
     if(a(ij).ne.0.)ti=ti+w(i)**2/a(ij)
7    ij=ij+np-i
     goto20
10   continue
     do 11 i=1,n
11      w(i)=z(i)
        do 15 i=1,n
           ip=i+1
           v=w(i)
           if(a(ij).gt.0.)goto12
      w(i)=0.
      ij=ij+np-i
      goto 15
12    continue
      ti=ti+v**2/a(ij)
      if(i.eq.n)goto14
      do 13 j=ip,n
      ij=ij+1
13    w(j)=w(j)-v*a(ij)
14    ij=ij+1
15    continue
20    continue
      if(ir.le.0 )goto 21
      if(ti.gt.0.)goto 22
      if(mk-1)40,40,23
21    ti=0.
      ir=-ir-1
      goto23
22    ti=eps/sig
      if(eps.eq.0.)ir=ir-1
23    continue
      mm=1
      tim=ti
      do 30 i=1,n
      j=np-i
      ij=ij-i
      if(a(ij).ne.0.)tim=ti-w(j)**2/a(ij)
      w(j)=ti
30    ti=tim
      goto 41
40    continue
      mm=0
      tim=1./sig
41    continue
      ij=1
      do 66 i=1,n
      ip=i+1
      v=z(i)
      if(a(ij).gt.0.)goto 53
      if(ir.gt.0 .or.sig.lt.0..or.v.eq.0.)goto 52
      ir=1-ir
      a(ij)=v**2/tim
      if(i.eq.n)return
      do 51 j=ip,n
      ij=ij+1
51    a(ij)=z(j)/v
      return
52    continue
      ti=tim
      ij=ij+np-i
      goto66
53    continue
      al=v/a(ij)
      if(mm)54,54,55
54    ti=tim+v*al
      goto56
55    ti=w(i)
56    continue
      r=ti/tim
      a(ij)=a(ij)*r
      if(r.eq.0.)goto 70
      if(i.eq.n)goto 70
      b=al/ti
      if(r.gt.4.)goto 62
      do 61 j=ip,n
      ij=ij+1
      z(j)=z(j)-v*a(ij)
61    a(ij)=a(ij)+b*z(j)
      goto64
62    gm=tim/ti
      do 63 j=ip,n
      ij=ij+1
      y=a(ij)
      a(ij)=b*z(j)+y*gm
63    z(j)=z(j)-v*y
64    continue
      tim=ti
      ij=ij+1
66    continue
70    continue
      if(ir.lt.0)ir=-ir
      return
!!c-----------------------------------multiply vector z by inverse of factors in a
      entry mc11e(a,n,z,w,ir)
   
      if(ir.lt.n)return
      w(1)=z(1)
      if(n.gt.1)goto 400
      z(1)=z(1)/a(1)
      return
400   continue
      do 402 i=2,n
      ij=i
      i1=i-1
      v=z(i)
      do 401 j=1,i1
         
         v=v-a(ij)*z(j)
        ! print*,n,j,ij,z(j)
401      ij=ij+n-j
         !print*,'a',ij
      w(i)=v
402   z(i)=v
      z(n)=z(n)/a(ij)
      !print*,'hello',ij,a(ij)
      np=n+1
      do 411 nip=2,n
      i=np-nip
      ii=ij-nip    
      v=z(i)/a(ii)
      ip=i+1
      ij=ii
      do 410 j=ip,n
         ii=ii+1
        ! print*,'ii',a(ii)
410   v=v-a(ii)*z(j)
      
411   z(i)=v
      
      return
    end subroutine mc11a

!!!===================
!    
!   GRADIENT TEST
!    
!!!===================
    subroutine gradt
      use defs
      use minimization
      use minimization1, only: nvar
      use md_defs, only:calc_SAXS_force, calc_SAXS_modul
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
      
      real(8),dimension(nvar) :: fsav 
     ! natoms3 = NATOMS*3


      if(use_ics)then
         call assemb
         call putbac(0)!Store coordinates and gradients
      endif
      
      dstep=1.d-6
      calc_SAXS_force = .true.
      calc_SAXS_modul = .false.
      call calcforce(scale,pos,force,total_energy,0.0d0)

      write(6,5) total_energy
5     format(/2x,'Gradient test .... etot= ',f10.3/)
      esav=total_energy
      
      k=0
      
      if(use_ics)then
         do i=1,nvar
            var(i)=var(i)+dstep
            call microb
            calc_SAXS_force = .true.
            calc_SAXS_modul = .false.
            call calcforce(scale,pos,force,total_energy,0.0d0)
            !print*,'grad',i,gra(i)
            grad=(total_energy-esav)/dstep
            !print*,'tt',total_energy,esav
            write(6,10) i,j,gra(i),grad,gra(i)-grad
10          format(2x,i4,'/',i1,') ',3f15.3)
            if (abs(gra(i)-grad) > 1.d-3) then
               write(6,*) 'alarm: ',i,j,total_energy,esav,dstep
            endif
            
            call putbac(1)!Put back the stored variables 
         enddo
     
      else
         fsav(1:nvar)=-force(1:nvar)
         do i=1,natoms
            do j=1,3
               k=k+1
               !print*,'pos before',pos(k)
               pos(k)=pos(k)+dstep
               !call rebuild
               calc_SAXS_force = .true.
               calc_SAXS_modul = .false.
               !print*,'new pos gradt',k,pos(k)
               call calcforce(scale,pos,force,total_energy,0.0d0)
               
               grad=(total_energy-esav)/dstep
               pos(k)=pos(k)-dstep
               write(6,10) i,j,fsav(k),grad,fsav(k)-grad
!10             format(2x,i4,'/',i1,') ',3f15.3)
               if(abs(fsav(k)-grad) > 1.d-3) then
                  write(6,*) 'alarm: ',i,j,total_energy,esav,dstep
               endif
               
            enddo
         enddo
      endif
    end subroutine gradt

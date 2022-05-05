subroutine bodies

  use defs, only : pos,NATOMS
  use rnabase
  use PDBtext, only : numres
  use fragments
  use int_defs
  use system_defs_RNA
  
  implicit none
  
  integer :: inb,nbodies,i,j,k,istart,iend,ij,ii,m1,jj,kstart,kend,i1,kk
  real(8):: rmin,xco,yco,zco,dx2,dy2,dz2,r
  real(8),dimension(natoms):: xx,yy,zz
  real(8),dimension(natoms,3):: com

  !!nfrag: total number of fragment
  
  inb=nfrag
  
  ! Move?
  !allocate(lc(nfrag))
  !allocate(kch(nfrag))
  !allocate(nbbvar(nfrag))
  !allocate(nbasepp(NRES))
  !allocate(ipiv(nfrag))
  !allocate(iseg1(nfrag*2))
  !allocate(iseg2(nfrag*2))
  !allocate(ibseg(nfrag*2))
  !allocate(log3t(nfrag*2))
  !allocate(log5t(nfrag*2))
  !Finding n° of rigid bodies
  !cha=menc(nen(1))
  !x(:)=0.d0
  !y(:)=0.d0
  !z(:)=0.d0
  !kc(:)=0!n° of interacting points in each body
  !lc(:)=0!last atom number of each body
  !ipiv(:)=0 !pivoting P for each body
  !ibody(:)=1

 
  do i=1,nfrag
     if(i.eq.1)then
        lc(i)=lenfrag(i)
     else
        lc(i)=lc(i-1)+lenfrag(i)
     endif
  
     !print*,i,lc(i)
  end do

  do i=1,NRES
     if ((btype(i).eq.1).or.(btype(i).eq.2))then !! A or G
        nbasepp(i)=2
     else
        nbasepp(i)=1
     endif
  enddo
  
  do i=1,nfrag
     do j=1,NRES
        if(blist(j).eq.lc(i))then
           !print*,'ii',j,blist(j),btype(j)
           kch(i)=j
           !print*,'i',kch(i)
        endif
     enddo
  enddo

  do i=1,nfrag
     if(i.eq.1)then
        ij=blist(1)-nbasepp(1)-1
     else
        ij=blist(kch(i-1)+1)-nbasepp(kch(i-1)+1)
     endif
     if(ij.eq.4)then
        nbbvar(i)=2
     else
        nbbvar(i)=3
     endif
     !print*,'i',nbbvar(i),i
  enddo

  xx(:)=0.d0
  yy(:)=0.d0
  zz(:)=0.d0
  do i=1,nfrag
     if(i.eq.1)then
        istart=1
     else
        istart=lc(i-1)+1
     endif
     iend=lc(i)
     do j=istart,iend
        xx(i)=xx(i)+pos((j-1)*3+1)
        yy(i)=yy(i)+pos((j-1)*3+2)
        zz(i)=zz(i)+pos((j-1)*3+3)
        !print*,i,xx(i),j,pos((j-1)*3+1),yy(i)
     enddo
  enddo
        
  nbodies=nfrag
  !print*,'leng',lenfrag(nbodies)

  !print*,"nbodies",nbodies
  !write(6,26) nbodies
!26 format(/2x,'Input is composed of ',i4,' rigid body(/ies) ')

  !print*,pos(2:3*natoms:3)
  !Finding pivoting residues for each body
  do k=1,nbodies
     com(k,1)=xx(k)/lenfrag(k)
     com(k,2)=yy(k)/lenfrag(k)
     com(k,3)=zz(k)/lenfrag(k)
     
     if(k.eq.1)then
        istart=1
     else
        istart=kch(k-1)+1
     endif

     rmin=1000.d0
     xco = com(k,1)
     yco = com(k,2)
     zco = com(k,3)
     !print*,'com',xco,istart,kch(k)
     do i=istart,kch(k)
        if(i.eq.istart)then
           if(i.eq.1)then
              ii=1
           else
              ii=lc(k-1)+1
           endif
        else
           ii=blist(i-1)+1
        endif
        !print*,'ii',ii
 
        dx2=(pos((ii-1)*3+1)-xco)**2
        dy2=(pos((ii-1)*3+2)-yco)**2
        dz2=(pos((ii-1)*3+3)-zco)**2
        r=sqrt(dx2+dy2+dz2)
     
        if(ii.eq.istart)then
           !print*,'FIRST',r
           rmin=r
        endif
   
        if(r.le.rmin)then
           rmin=r
           m1 = ii
           !print*,r,rmin,ii
        endif
     enddo
     !print*,k,m1
     ipiv(k) = m1
  enddo

  k=0
  log3t(1:nfrag*2)=.false.
  log5t(1:nfrag*2)=.false.
  do jj=1,nfrag
     !print*,numres(ipiv(jj)),jj,ipiv(jj)
     if(jj.eq.1)then
        i=1
        kk=1
     else
        i=lc(jj-1)+1
        kk=kch(jj-1)+1
     endif
     i1=numres(ipiv(jj))-numres(i)+kk
     !print*,'i1',i1,numres(ipiv(jj)),numres(i),i
     if(jj.eq.1)then
        kstart=1
     else
        kstart=kch(jj-1)+1
     endif
     kend=kch(jj)
     !print*,kstart,kend
     k=k+1
     iseg1(k)=kstart
     iseg2(k)=i1
     log5t(k)=.true.
     !print*,'i',iseg1(k),iseg2(k)
     ibseg(k)=jj
     k=k+1
     iseg1(k)=i1
     iseg2(k)=kend
     log3t(k)=.true.
     ibseg(k)=jj
     !print*,'i',iseg1(k),iseg2(k)
  enddo
 

end subroutine bodies

!====Atoms index======================
!
!   If not P, O5',C5'(+1),CA(+2)
!   If P, P,O5'(+1),C5'(+2),CA(+3)
!
!   DEFINITION OF INTERNAL VARIABLES
!
!=====================================

subroutine varis

  use defs, only : pos,NATOMS,use_back,var,nvar,use_base,ovr,svars,use_val
  use numerical_defs
  use system_defs_RNA
  use fragments
  use int_defs
  use rnabase
  implicit none
   
  integer :: k,jj,kio,j,ii,i1,i2,i3,i4,ia,jf,ji,kk,i,idx,itot,ncl
  real(8) :: torp,ang

  k=0
  nloop=nfrag*2
  itot=0
  
  !print*,'NRES',NRES,NFRAG
  if(use_back)then
     !print*,'itot',itot
     itot=itot+(NRES-2*(NFRAG))*4+3*NFRAG
     !print*,itot
     do i=1,NFRAG
        itot=itot+5-nbbvar(i)
     enddo
     if(use_val)then
        itot=itot+(NRES-2*(NFRAG))*4+3*NFRAG
        do i=1,NFRAG
           itot=itot+5-nbbvar(i)
        enddo
     endif
  endif
  !print*,itot
  if(use_base)then
     do i=1,NRES
        itot=itot+nbasepp(i)
     enddo
     if(use_val)then
        do i=1,NRES
           itot=itot+nbasepp(i)+1
        enddo
     endif
  endif
  !print*,'nres',NRES,itot
  nvar=itot
  !print*,'nvar',nvar

  ! Allocate stuff
  if (.not. allocated(var)) allocate(var(nvar))
  if (.not. allocated(ovr)) allocate(ovr(nvar))
  if (.not. allocated(svars)) allocate(svars(nvar))

  if(use_base)then
     do kk=1,NRES
        ncl=nbasepp(kk)
        do i=1,ncl
           !            if(i.eq.1)then
           k=k+1
           i1=blist(kk)-ncl+i-3
           i2=blist(kk)-ncl+i-2
           i3=blist(kk)-ncl+i-1
           i4=blist(kk)-ncl+i
           !print*,'i4',i1,i2,i3,i4
           var(k)=torp(i1,i2,i3,i4)*rad2deg
           ovr(k)=var(k)
           svars(k)=var(k)
           
           !print*,'var',k,var(k),i1,i2,i3,i4
           !          else
           !           endif
        enddo
     enddo
     if(use_val)then
        do kk=1,NRES
           k=k+1
           
           ncl=nbasepp(kk)
           i1=blist(kk)-ncl-2
           i2=blist(kk)-ncl-1
           i3=blist(kk)-ncl
           var(k)=ang(i1,i2,i3)*rad2deg
           ovr(k)=var(k)
           svars(k)=var(k)
           !print*,'var',k,var(k),i1,i2,i3,i4

           do i=1,ncl
              !            if(i.eq.1)then
              k=k+1
              i1=blist(kk)-ncl+i-2
              i2=blist(kk)-ncl+i-1
              i3=blist(kk)-ncl+i
              !print*,'i4',i1,i2,i3,i4
              var(k)=ang(i1,i2,i3)*rad2deg
              ovr(k)=var(k)
              svars(k)=var(k)
              
              !print*,'var',k,var(k),i1,i2,i3
              !          else
              !           endif
           enddo
        enddo
     endif
  endif
  
   if(use_back)then
      do jj=1,nloop
!=====Backbone torsions
         do ii=iseg1(jj),iseg2(jj)
            if(log5t(jj).and.(ii.eq.iseg2(jj)))cycle
            if(log5t(jj))then
               if(ii.eq.iseg1(jj))then
                  ji=nbbvar(ibseg(jj))
                  jf=4
                  !print*,'ji',ibseg(jj)
               else
                  ji=1
                  jf=4
               endif
            elseif(log3t(jj))then
               if(ii.eq.iseg2(jj))then
                  ji=1
                  jf=3
               else
                  ji=1
                  jf=4
               endif
            endif
            !print*,'log5t',log5t(jj),iseg1(jj),iseg2(jj),ii,ji,jf
            if(ii.eq.1) then
               idx=0
            else
               idx=blist(ii-1)
            endif
            !print*,ji,jf
            do j=ji,jf 
               k=k+1
               if(j.eq.1)then
                  i1=blist(ii-1)-nbasepp(ii-1)-1
                  i2=idx+1
                  i3=idx+2
                  i4=idx+3
                  !print*,j,i1,i2,i3,i4
               elseif(j.eq.2)then
                  i1=idx+1
                  i2=idx+2
                  i3=idx+3
                  i4=idx+4
                  !print*,j,i1,i2,i3,i4
               elseif(j.eq.3)then
                  i1=blist(ii)-nbasepp(ii)-3
                  i2=blist(ii)-nbasepp(ii)-2
                  i3=blist(ii)-nbasepp(ii)-1
                  if((ii.eq.iseg2(jj)).and.log3t(jj))then
                     i4=blist(ii)-nbasepp(ii)
                  else
                     i4=blist(ii)+1
                  endif
                  !print*,j,i1,i2,i3,i4
               elseif(j.eq.4)then
                  i1=blist(ii)-nbasepp(ii)-2
                  i2=blist(ii)-nbasepp(ii)-1
                  i3=blist(ii)+1
                  i4=blist(ii)+2
                  !print*,j,i1,i2,i3,i4
               endif
               var(k)=torp(i1,i2,i3,i4)*rad2deg
               !print*,'var',k,var(k),i1,i2,i3,i4,log3t(jj),ii
               ovr(k)=var(k)
               svars(k)=var(k)
            
            enddo
         enddo      
      enddo
      !====Valence angles 
      if(use_val)then
         do jj=1,nloop
            do ii=iseg1(jj),iseg2(jj)
          
               if(log5t(jj).and.(ii.eq.iseg2(jj)))cycle
               
               if(log5t(jj))then
                  if(ii.eq.iseg1(jj))then
                     ji=nbbvar(ibseg(jj))
                     jf=4
                     !print*,'ji',ibseg(jj)
                  else
                     ji=1
                     jf=4
                  endif
               elseif(log3t(jj))then
                  if(ii.eq.iseg2(jj))then
                     ji=1
                     jf=3
                  else
                     ji=1
                     jf=4
                  endif
               endif
               !print*,'ji',k,ii,ji,jf
               if(ii.eq.1)then
                  idx=0
               else
                  idx=blist(ii-1)
               endif
               do j=ji,jf 
                  k=k+1
                  
                  if(j.eq.1)then
                     i1=blist(ii-1)-nbasepp(ii-1)-1
                  else
                     
                     if(ii.eq.iseg1(jj))then
                      
                        i1=idx+j+1-nbbvar(ibseg(jj))
                       
                     else
                        i1=blist(ii-1)+(j-1)
                     endif
                  endif
                  if(ii.eq.iseg1(jj))then
                     i2=idx+j+2-nbbvar(ibseg(jj))
                  else
                     i2=blist(ii-1)+j
                  endif
                   
                  if(j.eq.4)then
                     i3=blist(ii)+1
                  else
                     if(ii.eq.iseg1(jj))then
                        i3=idx+j+3-nbbvar(ibseg(jj))
                     else
                        i3=blist(ii-1)+1+j
                     endif
                  endif
                  var(k)=ang(i1,i2,i3)*rad2deg
                  !print*,'var',k,var(k),i1,i2,i3
                  ovr(k)=var(k)
                  svars(k)=var(k)
               enddo
            enddo
         enddo
      endif
   endif
   !print*,'k',k,nvar

   nvar=k

 end subroutine varis

 
 !======================================!
 !                                      !
 !   Gradients on internal variables    !
 !                                      !
 !======================================!

 subroutine assemb
   use defs, only : pos,NATOMS,use_back,var,nvar,use_base,gra,force,use_val
   use numerical_defs
   use system_defs_RNA
   use fragments
   use int_defs
   use rnabase
   
   implicit none

   integer::i1,i2,i3,i4,k,ncl,kk,i,kinc,j,jj,ii,ji,jf,iii,kkk
   integer :: jja,jjf,idx,kioi,jji,kmove,kki,jl,jii
   real*8 :: frx,fry,frz,u,ux,uy,uz,x0,y0,z0,dot,dx,dz,dy,dx1,dx2,dy1,dy2,dz1,dz2
   real(8),dimension(nvar):: gl
   real(8),dimension(NATOMS*3) :: for
   
   k=0
   for(1:NATOMS*3)=-force(1:NATOMS*3)
   
   if(use_base)then
      do kk=1,NRES
         ncl=nbasepp(kk)
         do i=1,ncl
            !            if(i.eq.1)then
            k=k+1
         
            i1=blist(kk)-ncl+i-2
            i2=blist(kk)-ncl+i-1
        
            frx=0.
            fry=0.
            frz=0.
            x0=pos((i2-1)*3+1)  !pivot
            y0=pos((i2-1)*3+2)
            z0=pos((i2-1)*3+3)
            kinc=0
            do jj=i,ncl-1+1

               kinc=kinc+1

               j=i2+kinc

               dx=pos((j-1)*3+1)-x0
               dy=pos((j-1)*3+2)-y0
               dz=pos((j-1)*3+3)-z0
               
               frx=frx+dy*for((j-1)*3+3)-dz*for((j-1)*3+2)
               fry=fry+dz*for((j-1)*3+1)-dx*for((j-1)*3+3)
               frz=frz+dx*for((j-1)*3+2)-dy*for((j-1)*3+1)
               !print*,'j',j,mena(j)
            enddo
            ux=x0-pos((i1-1)*3+1) !axis
            uy=y0-pos((i1-1)*3+2)
            uz=z0-pos((i1-1)*3+3)
            !print*,"ux",ux,uy,uz
            !print*,"x0",x0,y0,z0
            u=sqrt(ux**2+uy**2+uz**2)
            dot  = (frx*ux+fry*uy+frz*uz)/u
            gl(k)= dot/rad2deg
            !print*,'gl',k,gl(k)
         enddo
      enddo
      if(use_val)then
         do kk=1,NRES
            k=k+1
                 
            ncl=nbasepp(kk)
            i1=blist(kk)-ncl-2
            i2=blist(kk)-ncl-1
            i3=blist(kk)-ncl
            x0=pos((i2-1)*3+1)  !pivot
            y0=pos((i2-1)*3+2)
            z0=pos((i2-1)*3+3)
            kinc=0
            frx=0.d0
            fry=0.d0
            frz=0.d0
            !print*,ncl,ncl-1+1
            do jj=1,ncl
               
               kinc=kinc+1

               j=i2+kinc

               dx=pos((j-1)*3+1)-x0
               dy=pos((j-1)*3+2)-y0
               dz=pos((j-1)*3+3)-z0
               
               frx=frx+dy*for((j-1)*3+3)-dz*for((j-1)*3+2)
               fry=fry+dz*for((j-1)*3+1)-dx*for((j-1)*3+3)
               frz=frz+dx*for((j-1)*3+2)-dy*for((j-1)*3+1)
              ! print*,'j',j
            enddo
            dx1=x0-pos((i1-1)*3+1) !b i-1
            dy1=y0-pos((i1-1)*3+2)
            dz1=z0-pos((i1-1)*3+3)
            dx2=x0-pos((i3-1)*3+1) !-bi 
            dy2=y0-pos((i3-1)*3+2)
            dz2=z0-pos((i3-1)*3+3)
            ux=dy2*dz1-dz2*dy1
            uy=dz2*dx1-dx2*dz1
            uz=dx2*dy1-dy2*dx1
            u=sqrt(ux**2+uy**2+uz**2)
            dot=(frx*ux+fry*uy+frz*uz)/u
            gl(k)= dot/rad2deg
            !print*,'gl val',k,gl(k)
            !print*,'var',k,var(k),i1,i2,i3,i4
            do i=1,ncl
               k=k+1
               i1=blist(kk)-ncl+i-2
               i2=blist(kk)-ncl+i-1
               i3=blist(kk)-ncl+i
               x0=pos((i2-1)*3+1)  !pivot
               y0=pos((i2-1)*3+2)
               z0=pos((i2-1)*3+3)
               kinc=0
               frx=0.d0
               fry=0.d0
               frz=0.d0
               do jj=i,ncl-1+1
                  
                   kinc=kinc+1
                   
                   j=i2+kinc
                   
                   dx=pos((j-1)*3+1)-x0
                   dy=pos((j-1)*3+2)-y0
                   dz=pos((j-1)*3+3)-z0
                   
                   frx=frx+dy*for((j-1)*3+3)-dz*for((j-1)*3+2)
                   fry=fry+dz*for((j-1)*3+1)-dx*for((j-1)*3+3)
                   frz=frz+dx*for((j-1)*3+2)-dy*for((j-1)*3+1)
               !print*,'j',j,mena(j)
                enddo
                dx1=x0-pos((i1-1)*3+1) !b i-1
                dy1=y0-pos((i1-1)*3+2)
                dz1=z0-pos((i1-1)*3+3)
                dx2=x0-pos((i3-1)*3+1) !-bi 
                dy2=y0-pos((i3-1)*3+2)
                dz2=z0-pos((i3-1)*3+3)
                ux=dy2*dz1-dz2*dy1
                uy=dz2*dx1-dx2*dz1
                uz=dx2*dy1-dy2*dx1
                u=sqrt(ux**2+uy**2+uz**2)
                dot=(frx*ux+fry*uy+frz*uz)/u
                gl(k)= dot/rad2deg
               ! print*,'gl val2',k,gl(k)!,u,dx1
             enddo
          enddo
       endif
    endif
    
   if(use_back)then
      do jj=1,nloop
         if(log3t(jj))then
            do ii=iseg1(jj),iseg2(jj)
           
               if(ii.eq.iseg2(jj))then
                  ji=1
                  jf=3
               else
                  ji=1
                  jf=4
               endif
               if(ii.eq.1) then
                  idx=0
               else
                  idx=blist(ii-1)
               endif
               do iii=ji,jf 
                  k=k+1
                  if(iii.eq.1)then
                     i1=idx+1
                     i2=idx+2                
                     !print*,j,i1,i2,i3,i4
                  elseif(iii.eq.2)then
                     
                     i1=idx+2
                     i2=idx+3
                     
                     !print*,j,i1,i2,i3,i4
                  elseif(iii.eq.3)then
                     
                     i1=blist(ii)-nbasepp(ii)-2
                     i2=blist(ii)-nbasepp(ii)-1
                     
                  !print*,j,i1,i2,i3,i4
                  elseif(iii.eq.4)then
                     
                     i1=blist(ii)-nbasepp(ii)-1
                     i2=blist(ii)+1
                     !print*,j,i1,i2,i3,i4
                  endif
               
                  x0=pos((i2-1)*3+1)  !pivot
                  y0=pos((i2-1)*3+2)
                  z0=pos((i2-1)*3+3)
                  frx=0.d0
                  fry=0.d0
                  frz=0.d0
                  do kkk=ii,iseg2(jj)
                     if((kkk.eq.ii).and.(iii.eq.4)) cycle
                     if(kkk.eq.ii)then
                        if((ii.eq.iseg2(jj)).and.(iii.eq.3))then
                           jja=5
                        else
                           jja=iii+2
                        endif
                     else
                        if(((kkk-ii).eq.1).and.(iii.eq.4))then
                           jja=2
                        else
                           jja=1
                        endif
                     endif
                     jjf=5+nbasepp(kkk)
                     do jji=jja,jjf
                        kioi=blist(kkk-1)+jji
                        dx=pos((kioi-1)*3+1)-x0
                        dy=pos((kioi-1)*3+2)-y0
                        dz=pos((kioi-1)*3+3)-z0
                        frx=frx+dy*for((kioi-1)*3+3)-dz*for((kioi-1)*3+2)
                        fry=fry+dz*for((kioi-1)*3+1)-dx*for((kioi-1)*3+3)
                        frz=frz+dx*for((kioi-1)*3+2)-dy*for((kioi-1)*3+1)
                        !if(iii.eq.3)
                        !print*,'index',i1,i2,kioi,iii
                     enddo
                  enddo
                  ux=x0-pos((i1-1)*3+1) !axis
                  uy=y0-pos((i1-1)*3+2)
                  uz=z0-pos((i1-1)*3+3)
                  !print*,"ux",ux,uy,uz
                  !print*,"x0",x0,y0,z0
                  u=sqrt(ux**2+uy**2+uz**2)
                  dot  = (frx*ux+fry*uy+frz*uz)/u
                  gl(k)= dot/rad2deg
                  !print*,'k',k,gl(k)
               enddo
            enddo
!!!===5' TERMINUS
         else
            
            do ii=iseg1(jj),iseg2(jj)
           
               if(ii.eq.1) then
                  idx=0
               else
                  idx=blist(ii-1)
               endif
         
               if(ii.eq.iseg1(jj))then
                  kmove=1
                  ji=nbbvar(ibseg(jj))
                  jf=4
               elseif(ii.eq.iseg2(jj))then
                  cycle
               else
                  ji=1
                  jf=4
               endif
               do iii=ji,jf
                  k=k+1
                 
                  if(iii.eq.1)then
                     i1=idx+1
                     i2=idx+2
                     !print*,j,i1,i2,i3,i4
                  elseif(iii.eq.2)then
                     i1=idx+2
                     i2=idx+3
                     !print*,j,i1,i2,i3,i4
                  elseif(iii.eq.3)then
                     i1=blist(ii)-nbasepp(ii)-2
                     i2=blist(ii)-nbasepp(ii)-1
                     !print*,j,i1,i2,i3,i4
                  elseif(iii.eq.4)then
                     i1=blist(ii)-nbasepp(ii)-1
                     i2=blist(ii)+1
                     !print*,j,i1,i2,i3,i4
                  endif
                  x0=pos((i1-1)*3+1)  !pivot
                  y0=pos((i1-1)*3+2)
                  z0=pos((i1-1)*3+3)
                  frx=0.d0
                  fry=0.d0
                  frz=0.d0
                  do kki=1,kmove
                     kkk=iseg1(jj)+(kki-1)
                     if(kkk.eq.1)then
                        kioi=0
                     else
                        kioi=blist(kkk-1)
                     endif
                     !print*,ii,kki,kkk,kioi,kmove
                     if(kmove.eq.kki)then
                        if(iii.eq.1)then
                           cycle
                        elseif(iii.eq.2)then
                           jja=1
                           jjf=1
                        elseif(iii.eq.3)then
                           if(kki.eq.1)then
                              jja=1
                              jjf=4-nbbvar(ibseg(jj))
                              !print*,jja,jjf
                           else
                              jja=1
                              jjf=2
                           endif
                        else
                           jja=1
                           if(kki.eq.1)then
                              jjf=7-nbbvar(ibseg(jj))+nbasepp(kkk)
                              !print*,jja,jjf,nbbvar(ibseg(jj)),kkk
                           else
                              jjf=5+nbasepp(kkk)
                           endif
                        endif
                     else
                        jja=1
                        if(kki.eq.1)then
                           jjf=7-nbbvar(ibseg(jj))+nbasepp(kkk)
                           !print*,jja,jjf,nbbvar(ibseg(jj)),kkk
                        else
                           jjf=5+nbasepp(kkk)
                        endif
                     endif
                     do jl=jja,jjf
                        if((kki.eq.kmove).and.(iii.eq.4))then
                           if((kkk.eq.iseg1(jj)).and.(nbbvar(ibseg(jj)).eq.3))then
                              if(jl.eq.3) cycle
                           else
                              if(jl.eq.4)cycle
                           endif
                        endif
                        
                        jii=kioi+jl
                        !print*,'move',ii,i1,i2,jii
                        
                        dx=pos((jii-1)*3+1)-x0
                        dy=pos((jii-1)*3+2)-y0
                        dz=pos((jii-1)*3+3)-z0

                        frx=frx+dy*for((jii-1)*3+3)-dz*for((jii-1)*3+2)
                        fry=fry+dz*for((jii-1)*3+1)-dx*for((jii-1)*3+3)
                        frz=frz+dx*for((jii-1)*3+2)-dy*for((jii-1)*3+1)
                     enddo                                                                      
                  enddo
                  ux=x0-pos((i2-1)*3+1) !axis
                  uy=y0-pos((i2-1)*3+2)
                  uz=z0-pos((i2-1)*3+3)
                  !print*,"ux",ux,uy,uz
                  !print*,"x0",x0,y0,z0
                  u=sqrt(ux**2+uy**2+uz**2)
                  dot  = (frx*ux+fry*uy+frz*uz)/u
                  gl(k)= dot/rad2deg
                  !print*,'k',gl(k)
                  if(iii.eq.4)then
                     kmove=kmove+1 
                  endif
                  
               enddo
            enddo
         endif
         
      enddo
   endif
   if(use_back.and.use_val)then
      do jj=1,nloop
         if(log3t(jj))then
            do ii=iseg1(jj),iseg2(jj)
               if(ii.eq.iseg2(jj))then
                  ji=1
                  jf=3
               else
                  ji=1
                  jf=4
               endif
            
               !print*,'ji',k,ii,ji,jf
               if(ii.eq.1)then
                  idx=0
               else
                  idx=blist(ii-1)
               endif

               do iii=ji,jf 
                  k=k+1
                  
                  if(iii.eq.1)then
                     i1=blist(ii-1)-nbasepp(ii-1)-1
                  else
                     i1=blist(ii-1)+(iii-1)
                  endif
               
                  i2=blist(ii-1)+iii
                  
                  if(iii.eq.4)then
                     i3=blist(ii)+1
                  else
                     i3=blist(ii-1)+1+iii  
                  endif
                  
                  frx=0.
                  fry=0.
                  frz=0.
                  x0=pos((i2-1)*3+1)  !pivot
                  y0=pos((i2-1)*3+2)
                  z0=pos((i2-1)*3+3)
                  do kkk=ii,iseg2(jj)
                     if((kkk.eq.ii).and.(iii.eq.4)) cycle
                     if(kkk.eq.ii)then
                        jja=iii+1                       
                     else
                       
                        jja=1
                        
                     endif
                     jjf=5+nbasepp(kkk)
                     do jji=jja,jjf
                        kioi=blist(kkk-1)+jji
                        dx=pos((kioi-1)*3+1)-x0
                        dy=pos((kioi-1)*3+2)-y0
                        dz=pos((kioi-1)*3+3)-z0
                        frx=frx+dy*for((kioi-1)*3+3)-dz*for((kioi-1)*3+2)
                        fry=fry+dz*for((kioi-1)*3+1)-dx*for((kioi-1)*3+3)
                        frz=frz+dx*for((kioi-1)*3+2)-dy*for((kioi-1)*3+1)
                        !if(iii.eq.3)
                        !print*,'index',k,i1,i2,kioi,iii
                     enddo
                  enddo

                  !===Definition unit vector 
                  dx1=x0-pos((i1-1)*3+1) !b i-1
                  dy1=y0-pos((i1-1)*3+2)
                  dz1=z0-pos((i1-1)*3+3)
                  dx2=x0-pos((i3-1)*3+1) !-bi 
                  dy2=y0-pos((i3-1)*3+2)
                  dz2=z0-pos((i3-1)*3+3)
                  ux=dy2*dz1-dz2*dy1
                  uy=dz2*dx1-dx2*dz1
                  uz=dx2*dy1-dy2*dx1
                  u=sqrt(ux**2+uy**2+uz**2)
                  dot=(frx*ux+fry*uy+frz*uz)/u
                  gl(k)= dot/rad2deg
                  !print*,'gl val back',k,gl(k)!,u,dx1
               enddo
            enddo
         else
  
!====5 terminus
            do ii=iseg1(jj),iseg2(jj)
               if(ii.eq.iseg2(jj))cycle
               
                if(ii.eq.iseg1(jj))then
                  kmove=1
                  ji=nbbvar(ibseg(jj))
                  jf=4
               elseif(ii.eq.iseg2(jj))then
                  cycle
               else
                  ji=1
                  jf=4
               endif
          
         
               !print*,'ji',k,ii,ji,jf
               if(ii.eq.1)then
                  idx=0
               else
                  idx=blist(ii-1)
               endif
               do iii=ji,jf 
                  k=k+1
                  
                  if(iii.eq.1)then
                     i1=blist(ii-1)-nbasepp(ii-1)-1
                  else  
                     if(ii.eq.iseg1(jj))then
                        i1=idx+iii+1-nbbvar(ibseg(jj))
                     else
                        i1=blist(ii-1)+(iii-1)
                     endif
                  endif
                  if(ii.eq.iseg1(jj))then
                     i2=idx+iii+2-nbbvar(ibseg(jj))
                     !print*,'i2',i2
                  else
                     i2=blist(ii-1)+iii
                  endif
                  if(iii.eq.4)then
                     i3=blist(ii)+1
                   
                  else
                     if(ii.eq.iseg1(jj))then
                        i3=idx+iii+3-nbbvar(ibseg(jj))
                     else
                        i3=blist(ii-1)+1+iii
                     endif
                  endif
                               
                  frx=0.
                  fry=0.
                  frz=0.
                  x0=pos((i2-1)*3+1)  !pivot
                  y0=pos((i2-1)*3+2)
                  z0=pos((i2-1)*3+3)
                  do kki=1,kmove
                     kkk=iseg1(jj)+(kki-1)
                     if(kkk.eq.1)then
                        kioi=0
                     else
                        kioi=blist(kkk-1)
                     endif
                     if(kmove.eq.kki)then
                        if(iii.eq.1)then
                           cycle
                        elseif(iii.eq.2)then
                           jja=1
                           jjf=1
                        elseif(iii.eq.3)then
                           if(kki.eq.1)then
                              jja=1
                              jjf=4-nbbvar(ibseg(jj))
                           else
                              jja=1
                              jjf=2
                           endif
                        else
                           if(kki.eq.1)then
                              jja=1
                              jjf=5-nbbvar(ibseg(jj))
                           else
                              jja=1
                              jjf=3
                           endif
                        endif
                     else
                        jja=1
                        if(kki.eq.1)then
                           jjf=7-nbbvar(ibseg(jj))+nbasepp(kkk)
                        else
                           jjf=5+nbasepp(kkk)
                        endif
                     endif
                     do jl=jja,jjf
                        jii=kioi+jl
                        !print*,'k',k,i1,i2,i3,jii
                        dx=pos((jii-1)*3+1)-x0
                        dy=pos((jii-1)*3+2)-y0
                        dz=pos((jii-1)*3+3)-z0
                        frx=frx+dy*for((jii-1)*3+3)-dz*for((jii-1)*3+2)
                        fry=fry+dz*for((jii-1)*3+1)-dx*for((jii-1)*3+3)
                        frz=frz+dx*for((jii-1)*3+2)-dy*for((jii-1)*3+1)
                     enddo
                     !===Definition unit vector 
                     dx1=x0-pos((i1-1)*3+1) !b i-1
                     dy1=y0-pos((i1-1)*3+2)
                     dz1=z0-pos((i1-1)*3+3)
                     dx2=x0-pos((i3-1)*3+1) !-bi 
                     dy2=y0-pos((i3-1)*3+2)
                     dz2=z0-pos((i3-1)*3+3)
                     ux=dy2*dz1-dz2*dy1
                     uy=dz2*dx1-dx2*dz1
                     uz=dx2*dy1-dy2*dx1
                     u=sqrt(ux**2+uy**2+uz**2)
                     dot=-(frx*ux+fry*uy+frz*uz)/u
                     gl(k)= dot/rad2deg
                     !print*,'gl val back',k,gl(k)!,u,dx1
                     
                  enddo
                  if(iii.eq.4)then
                     kmove=kmove+1 
                  endif
               enddo
            enddo
         endif
      enddo
   endif
   
   gra(1:nvar)=gl(1:nvar)
   !do i=1,nvar
   !   print*,'gra',gra(i)
   !enddo
   
   
 end subroutine assemb

 
 !!!===Microb: to build the structure after each step of internal coordinates minimization 

 subroutine microb
   use defs, only : pos,NATOMS,use_back,var,nvar,use_base,ovr,use_val
   use numerical_defs
   use system_defs_RNA
   use fragments
   use int_defs
   use rnabase
   
   implicit none

   integer::i1,i2,i3,i4,k,ncl,kk,i,kinc,j,jj,ii,ji,jf,iii,kkk
   integer :: jja,jjf,idx,kioi,jji,sss,kmove,kki,jl,jii,jjj
   real(8) :: frx,fry,frz,u,ux,uy,uz,x0,y0,z0,dot,dx,dz,dy,del,c,s
   real(8) :: dx1,dx2,dy1,dy2,dz1,dz2
   real(8) :: xx,yy,zz,rx,ry,rz
   
   !print*,'here'
   k=0
   if(use_base)then
      do kk=1,NRES
         ncl=nbasepp(kk)
         do i=1,ncl
            !if(i.eq.1)then
            k=k+1

            del=var(k)-ovr(k)
            ovr(k) = var(k)
            !print*,'delta',del,k
            del=del*deg2rad
            c=cos(del)
            s=sin(del)
            
            i1=blist(kk)-ncl+i-2
            i2=blist(kk)-ncl+i-1

            x0=pos((i2-1)*3+1)  !pivot
            y0=pos((i2-1)*3+2)
            z0=pos((i2-1)*3+3)
            ux=x0-pos((i1-1)*3+1)
            uy=y0-pos((i1-1)*3+2)
            uz=z0-pos((i1-1)*3+3)
            u=sqrt(ux**2+uy**2+uz**2)
            ux = ux/u     !axis
            uy = uy/u
            uz = uz/u
            sss=0
            ! print*,'pivot',i2,mena(i2),ncl
            kinc=0
            do jj=i,ncl-1+1

               kinc=kinc+1

               j=i2+kinc

               xx=pos((j-1)*3+1)-x0
               yy=pos((j-1)*3+2)-y0
               zz=pos((j-1)*3+3)-z0
               rx=(c+ux**2*(1-c))*xx+(ux*uy*(1-c)-uz*s)*yy+(uy*s+ux*uz*(1-c))*zz
               ry=(uz*s+ux*uy*(1-c))*xx+(c+uy**2*(1-c))*yy+(-ux*s+uy*uz*(1-c))*zz
               rz=(-uy*s+ux*uz*(1-c))*xx+(ux*s+uy*uz*(1-c))*yy+(c+uz**2*(1-c))*zz
               pos((j-1)*3+1)=x0+rx
               pos((j-1)*3+2)=y0+ry
               pos((j-1)*3+3)=z0+rz
               !print*,'j',j,mena(j)
            enddo
         
         enddo
      enddo
      if(use_val)then
         do kk=1,NRES
            ncl=nbasepp(kk)
            k=k+1
            del=var(k)-ovr(k)
            ovr(k) = var(k)
            !print*,'delta',del,k
            del=del*deg2rad
            c=cos(del)
            s=sin(del)
            i1=blist(kk)-ncl-2
            i2=blist(kk)-ncl-1
            i3=blist(kk)-ncl
            x0=pos((i2-1)*3+1)  !pivot
            y0=pos((i2-1)*3+2)
            z0=pos((i2-1)*3+3)
            dx1=x0-pos((i1-1)*3+1) !b i-1
            dy1=y0-pos((i1-1)*3+2)
            dz1=z0-pos((i1-1)*3+3)
            dx2=x0-pos((i3-1)*3+1) !-bi 
            dy2=y0-pos((i3-1)*3+2)
            dz2=z0-pos((i3-1)*3+3)
            ux=dy2*dz1-dz2*dy1
            uy=dz2*dx1-dx2*dz1
            uz=dx2*dy1-dy2*dx1
            u=sqrt(ux**2+uy**2+uz**2)
            ux=ux/u
            uy=uy/u
            uz=uz/u
            kinc=0
        
            !print*,ncl,ncl-1+1
            do jj=1,ncl
               
               kinc=kinc+1

               j=i2+kinc
               
               xx=pos((j-1)*3+1)-x0
               yy=pos((j-1)*3+2)-y0
               zz=pos((j-1)*3+3)-z0
               rx=(c+ux**2*(1-c))*xx+(ux*uy*(1-c)-uz*s)*yy+(uy*s+ux*uz*(1-c))*zz
               ry=(uz*s+ux*uy*(1-c))*xx+(c+uy**2*(1-c))*yy+(-ux*s+uy*uz*(1-c))*zz
               rz=(-uy*s+ux*uz*(1-c))*xx+(ux*s+uy*uz*(1-c))*yy+(c+uz**2*(1-c))*zz
               pos((j-1)*3+1)=x0+rx
               pos((j-1)*3+2)=y0+ry
               pos((j-1)*3+3)=z0+rz
              
              ! print*,'j',j
            enddo
            do i=1,ncl
               k=k+1
               del=var(k)-ovr(k)
               ovr(k) = var(k)
               !print*,'delta',del,k
               del=del*deg2rad
               c=cos(del)
               s=sin(del)
               i1=blist(kk)-ncl+i-2
               i2=blist(kk)-ncl+i-1
               i3=blist(kk)-ncl+i
               x0=pos((i2-1)*3+1)  !pivot
               y0=pos((i2-1)*3+2)
               z0=pos((i2-1)*3+3)
               dx1=x0-pos((i1-1)*3+1) !b i-1
               dy1=y0-pos((i1-1)*3+2)
               dz1=z0-pos((i1-1)*3+3)
               dx2=x0-pos((i3-1)*3+1) !-bi 
               dy2=y0-pos((i3-1)*3+2)
               dz2=z0-pos((i3-1)*3+3)
               ux=dy2*dz1-dz2*dy1
               uy=dz2*dx1-dx2*dz1
               uz=dx2*dy1-dy2*dx1
               u=sqrt(ux**2+uy**2+uz**2)
               kinc=0
               ux=ux/u
               uy=uy/u
               uz=uz/u
               do jj=i,ncl-1+1
                  
                  kinc=kinc+1
                  
                  j=i2+kinc
                  xx=pos((j-1)*3+1)-x0
                  yy=pos((j-1)*3+2)-y0
                  zz=pos((j-1)*3+3)-z0
                  rx=(c+ux**2*(1-c))*xx+(ux*uy*(1-c)-uz*s)*yy+(uy*s+ux*uz*(1-c))*zz
                  ry=(uz*s+ux*uy*(1-c))*xx+(c+uy**2*(1-c))*yy+(-ux*s+uy*uz*(1-c))*zz
                  rz=(-uy*s+ux*uz*(1-c))*xx+(ux*s+uy*uz*(1-c))*yy+(c+uz**2*(1-c))*zz
                  pos((j-1)*3+1)=x0+rx
                  pos((j-1)*3+2)=y0+ry
                  pos((j-1)*3+3)=z0+rz
                  
                    !print*,'j',j,mena(j)
               enddo
            enddo
       
         enddo
      endif
   endif

   if(use_back)then
      do jjj=1,nloop
         if(log3t(jjj))then
            do ii=iseg1(jjj),iseg2(jjj)
               !print*,"3t",ii
               if(ii.eq.iseg2(jjj))then
                  ji=1
                  jf=3
               else
                  ji=1
                  jf=4
               endif
               if(ii.eq.1) then
                  idx=0
               else
                  idx=blist(ii-1)
               endif
               do iii=ji,jf
                  k=k+1
                  del=var(k)-ovr(k)
                  ovr(k) = var(k)
                  !print*,'delta',del,k
                  del=del*deg2rad
                  !print*,'del',del
                  c=cos(del)
                  s=sin(del)
                  
                  if(iii.eq.1)then
                     i1=idx+1
                     i2=idx+2
                     !print*,j,i1,i2,i3,i4
                  elseif(iii.eq.2)then

                     i1=idx+2
                     i2=idx+3

                     !print*,j,i1,i2,i3,i4
                  elseif(iii.eq.3)then
                     
                     i1=blist(ii)-nbasepp(ii)-2
                     i2=blist(ii)-nbasepp(ii)-1

                     !print*,j,i1,i2,i3,i4
                  elseif(iii.eq.4)then

                     i1=blist(ii)-nbasepp(ii)-1
                     i2=blist(ii)+1
                     !print*,j,i1,i2,i3,i4
                  endif
                  
                  x0=pos((i2-1)*3+1)  !pivot
                  y0=pos((i2-1)*3+2)
                  z0=pos((i2-1)*3+3)
                  ux=x0-pos((i1-1)*3+1)
                  uy=y0-pos((i1-1)*3+2)
                  uz=z0-pos((i1-1)*3+3)
                  u=sqrt(ux**2+uy**2+uz**2)
                  ux = ux/u     !axis
                  uy = uy/u
                  uz = uz/u
                  
                  do kkk=ii,iseg2(jjj)
                     if((kkk.eq.ii).and.(iii.eq.4)) cycle
                     if(kkk.eq.ii)then
                        if((ii.eq.iseg2(jjj)).and.(iii.eq.3))then
                           jja=5
                        else
                           jja=iii+2
                        endif
                     else
                        if(((kkk-ii).eq.1).and.(iii.eq.4))then
                           jja=2
                        else
                           jja=1
                        endif
                     endif
                     jjf=5+nbasepp(kkk)
                     do jji=jja,jjf
                        kioi=blist(kkk-1)+jji

                        xx=pos((kioi-1)*3+1)-x0
                        yy=pos((kioi-1)*3+2)-y0
                        zz=pos((kioi-1)*3+3)-z0
                        !print*,'pos',pos((kioi-1)*3+1)
                        
                        rx=(c+ux**2*(1-c))*xx+(ux*uy*(1-c)-uz*s)*yy+(uy*s+ux*uz*(1-c))*zz
                        ry=(uz*s+ux*uy*(1-c))*xx+(c+uy**2*(1-c))*yy+(-ux*s+uy*uz*(1-c))*zz
                        rz=(-uy*s+ux*uz*(1-c))*xx+(ux*s+uy*uz*(1-c))*yy+(c+uz**2*(1-c))*zz
                        pos((kioi-1)*3+1)=x0+rx
                        pos((kioi-1)*3+2)=y0+ry
                        pos((kioi-1)*3+3)=z0+rz
                        
                        !print*,'at',ii,i1,i2,kioi,'pos',pos((kioi-1)*3+1),c,s,ux,uy
                     enddo
                  enddo
               enddo
            enddo
         else
            
            do ii=iseg1(jjj),iseg2(jjj)
               if(ii.eq.1) then
                  idx=0
               else
                  idx=blist(ii-1)
               endif

               if(ii.eq.iseg1(jjj))then
                  kmove=1
                  ji=nbbvar(ibseg(jjj))
                  jf=4
               elseif(ii.eq.iseg2(jjj))then
                  cycle
               else
                  ji=1
                  jf=4
               endif
               !print*,'NEW K',k,ii,iseg1(jjj),iseg2(jjj)
               do iii=ji,jf
                  k=k+1
                  del=var(k)-ovr(k)
                  ovr(k) = var(k)
                  !print*,'delta 5',del,k,log5t(jjj)
                  del=del*deg2rad
                  c=cos(del)
                  s=sin(del)
                  if(iii.eq.1)then
                     i1=idx+1
                     i2=idx+2
                     !print*,j,i1,i2,i3,i4
                  elseif(iii.eq.2)then
                     i1=idx+2
                     i2=idx+3
                     !print*,j,i1,i2,i3,i4
                  elseif(iii.eq.3)then
                     i1=blist(ii)-nbasepp(ii)-2
                     i2=blist(ii)-nbasepp(ii)-1
                     !print*,j,i1,i2,i3,i4
                  elseif(iii.eq.4)then
                     i1=blist(ii)-nbasepp(ii)-1
                     i2=blist(ii)+1
                     !print*,j,i1,i2,i3,i4
                  endif
                  x0=pos((i1-1)*3+1)  !pivot
                  y0=pos((i1-1)*3+2)
                  z0=pos((i1-1)*3+3)
                  ux=x0-pos((i2-1)*3+1)
                  uy=y0-pos((i2-1)*3+2)
                  uz=z0-pos((i2-1)*3+3)
                  u=sqrt(ux**2+uy**2+uz**2)
                  ux = ux/u        !axis
                  uy = uy/u
                  uz = uz/u
                  sss=0
                 
                  do kki=1,kmove
                     kkk=iseg1(jjj)+(kki-1)
                     !print*,'kmove',kmove,'kkk',kkk,'iseg 1',iseg1(jjj),kki,iseg1(jjj)+(kki-1)
                     if(kkk.eq.1)then
                        kioi=0
                     else
                        kioi=blist(kkk-1)
                     endif
                     
                     if(kmove.eq.kki)then                                               
                        if(iii.eq.1)then                                                
                           cycle                                                        
                        elseif(iii.eq.2)then                                            
                           jja=1                                                        
                           jjf=1                                                        
                        elseif(iii.eq.3)then                                            
                           if(kki.eq.1)then                                             
                              jja=1                                                     
                              jjf=4-nbbvar(ibseg(jjj))                                   
                              !print*,jja,jjf                                           
                           else                                                         
                              jja=1                                                     
                              jjf=2                                                     
                           endif                                                        
                        else                                                            
                           jja=1                                                        
                           if(kki.eq.1)then                                             
                              jjf=7-nbbvar(ibseg(jjj))+nbasepp(kkk)                      
                              !print*,jja,jjf,nbbvar(ibseg(jj)),kkk                     
                           else                                                         
                              jjf=5+nbasepp(kkk)                                        
                           endif                                                        
                        endif                                                           
                     else                                                               
                        jja=1                                                           
                        if(kki.eq.1)then                                                
                           jjf=7-nbbvar(ibseg(jjj))+nbasepp(kkk)                         
                           !print*,jja,jjf,nbbvar(ibseg(jj)),kkk                        
                        else                                                            
                           jjf=5+nbasepp(kkk)                                           
                        endif                                                           
                     endif
                     !print*,'jjf',jjf,nbasepp(kkk),kkk
                     do jl=jja,jjf
                        !print*,'jl',jl,jja,jjf
                        if((kki.eq.kmove).and.(iii.eq.4))then                           
                           if((kkk.eq.iseg1(jjj)).and.(nbbvar(ibseg(jjj)).eq.3))then      
                              if(jl.eq.3) cycle                                         
                           else                                                         
                              if(jl.eq.4)cycle                                          
                           endif                                                        
                        endif                                                           
                                                                                        
                        jii=kioi+jl                                                     
                        !print*,'move',ii,i1,i2,jii
                        
                        xx=pos((jii-1)*3+1)-x0
                        yy=pos((jii-1)*3+2)-y0
                        zz=pos((jii-1)*3+3)-z0

                        rx=(c+ux**2*(1-c))*xx+(ux*uy*(1-c)-uz*s)*yy+(uy*s+ux*uz*(1-c))*zz
                        ry=(uz*s+ux*uy*(1-c))*xx+(c+uy**2*(1-c))*yy+(-ux*s+uy*uz*(1-c))*zz
                        rz=(-uy*s+ux*uz*(1-c))*xx+(ux*s+uy*uz*(1-c))*yy+(c+uz**2*(1-c))*zz
                        pos((jii-1)*3+1)=x0+rx
                        pos((jii-1)*3+2)=y0+ry
                        pos((jii-1)*3+3)=z0+rz
                         
                     enddo
                     
                  enddo
                  if(iii.eq.4)then
                     kmove=kmove+1
                  endif
               enddo
            enddo
         endif
      enddo
   endif
   if(use_back.and.use_val)then
      do jj=1,nloop
         if(log3t(jj))then
            do ii=iseg1(jj),iseg2(jj)
               if(ii.eq.iseg2(jj))then
                  ji=1
                  jf=3
               else
                  ji=1
                  jf=4
               endif
            
               !print*,'ji',k,ii,ji,jf
               if(ii.eq.1)then
                  idx=0
               else
                  idx=blist(ii-1)
               endif

        
               do iii=ji,jf
                  k=k+1
                  del=var(k)-ovr(k)
                  ovr(k) = var(k)
                  !print*,'delta',del,k
                  del=del*deg2rad
                  c=cos(del)
                  s=sin(del)

                  if(iii.eq.1)then
                     i1=blist(ii-1)-nbasepp(ii-1)-1
                  else
                     i1=blist(ii-1)+(iii-1)
                  endif
                  
                  i2=blist(ii-1)+iii
                  
                  if(iii.eq.4)then
                     i3=blist(ii)+1
                  else
                     i3=blist(ii-1)+1+iii  
                  endif
                  x0=pos((i2-1)*3+1)  !pivot
                  y0=pos((i2-1)*3+2)
                  z0=pos((i2-1)*3+3)
                  dx1=x0-pos((i1-1)*3+1) !b i-1
                  dy1=y0-pos((i1-1)*3+2)
                  dz1=z0-pos((i1-1)*3+3)
                  dx2=x0-pos((i3-1)*3+1) !-bi 
                  dy2=y0-pos((i3-1)*3+2)
                  dz2=z0-pos((i3-1)*3+3)
                  ux=dy2*dz1-dz2*dy1
                  uy=dz2*dx1-dx2*dz1
                  uz=dx2*dy1-dy2*dx1
                  u=sqrt(ux**2+uy**2+uz**2)
                  ux=ux/u
                  uy=uy/u
                  uz=uz/u
                  do kkk=ii,iseg2(jj)
                     if((kkk.eq.ii).and.(iii.eq.4)) cycle
                     if(kkk.eq.ii)then
                        jja=iii+1                       
                     else
                        
                        jja=1
                        
                     endif
                     jjf=5+nbasepp(kkk)
                     do jji=jja,jjf
                        
                        kioi=blist(kkk-1)+jji
                        xx=pos((kioi-1)*3+1)-x0
                        yy=pos((kioi-1)*3+2)-y0
                        zz=pos((kioi-1)*3+3)-z0
                        rx=(c+ux**2*(1-c))*xx+(ux*uy*(1-c)-uz*s)*yy+(uy*s+ux*uz*(1-c))*zz
                        ry=(uz*s+ux*uy*(1-c))*xx+(c+uy**2*(1-c))*yy+(-ux*s+uy*uz*(1-c))*zz
                        rz=(-uy*s+ux*uz*(1-c))*xx+(ux*s+uy*uz*(1-c))*yy+(c+uz**2*(1-c))*zz
                        pos((kioi-1)*3+1)=x0+rx
                        pos((kioi-1)*3+2)=y0+ry
                        pos((kioi-1)*3+3)=z0+rz
                        !print*,'check',k,i1,i2,i3,kioi
                        !if(iii.eq.3)
                        !print*,'index',k,i1,i2,kioi,iii
                     enddo
                  enddo
               enddo
            enddo
         else
            do ii=iseg1(jj),iseg2(jj)
               if(ii.eq.iseg2(jj))cycle
               
               if(ii.eq.iseg1(jj))then
                  kmove=1
                  ji=nbbvar(ibseg(jj))
                  jf=4
               elseif(ii.eq.iseg2(jj))then
                  cycle
               else
                  ji=1
                  jf=4
               endif
          
         
               !print*,'ji',k,ii,ji,jf
               if(ii.eq.1)then
                  idx=0
               else
                  idx=blist(ii-1)
               endif
               do iii=ji,jf 
                  k=k+1
                  del=var(k)-ovr(k)
                  ovr(k) = var(k)
                  !print*,'delta',del
                  del=del*deg2rad
                  !print*,'delta',del,k

                  c=cos(del)
                  s=sin(del)
                  
                  if(iii.eq.1)then
                     i1=blist(ii-1)-nbasepp(ii-1)-1
                  else  
                     if(ii.eq.iseg1(jj))then
                        i1=idx+iii+1-nbbvar(ibseg(jj))
                     else
                        i1=blist(ii-1)+(iii-1)
                     endif
                  endif
                  if(ii.eq.iseg1(jj))then
                     i2=idx+iii+2-nbbvar(ibseg(jj))
                     !print*,'i2',i2
                  else
                     i2=blist(ii-1)+iii
                  endif
                  if(iii.eq.4)then
                     i3=blist(ii)+1
                  else
                     if(ii.eq.iseg1(jj))then
                        i3=idx+iii+3-nbbvar(ibseg(jj))
                     else
                        i3=blist(ii-1)+1+iii
                     endif
                  endif
                               
               
                  x0=pos((i2-1)*3+1)  !pivot
                  y0=pos((i2-1)*3+2)
                  z0=pos((i2-1)*3+3)
                  dx1=x0-pos((i1-1)*3+1) !b i-1
                  dy1=y0-pos((i1-1)*3+2)
                  dz1=z0-pos((i1-1)*3+3)
                  dx2=x0-pos((i3-1)*3+1) !-bi 
                  dy2=y0-pos((i3-1)*3+2)
                  dz2=z0-pos((i3-1)*3+3)
                  ux=dy2*dz1-dz2*dy1
                  uy=dz2*dx1-dx2*dz1
                  uz=dx2*dy1-dy2*dx1
                  u=sqrt(ux**2+uy**2+uz**2)
                  ux=-ux/u
                  uy=-uy/u
                  uz=-uz/u
                  do kki=1,kmove
                     kkk=iseg1(jj)+(kki-1)
                     if(kkk.eq.1)then
                        kioi=0
                     else
                        kioi=blist(kkk-1)
                     endif
                     if(kmove.eq.kki)then
                        if(iii.eq.1)then
                           cycle
                        elseif(iii.eq.2)then
                           jja=1
                           jjf=1
                        elseif(iii.eq.3)then
                           if(kki.eq.1)then
                              jja=1
                              jjf=4-nbbvar(ibseg(jj))
                           else
                              jja=1
                              jjf=2
                           endif
                        else
                           if(kki.eq.1)then
                              jja=1
                              jjf=5-nbbvar(ibseg(jj))
                           else
                              jja=1
                              jjf=3
                           endif
                        endif
                     else
                        jja=1
                        if(kki.eq.1)then
                           jjf=7-nbbvar(ibseg(jj))+nbasepp(kkk)
                        else
                           jjf=5+nbasepp(kkk)
                        endif
                     endif
                     do jl=jja,jjf
                        jii=kioi+jl
                       ! print*,'k',k,i1,i2,i3,jii
                        xx=pos((jii-1)*3+1)-x0
                        yy=pos((jii-1)*3+2)-y0
                        zz=pos((jii-1)*3+3)-z0
                        rx=(c+ux**2*(1-c))*xx+(ux*uy*(1-c)-uz*s)*yy+(uy*s+ux*uz*(1-c))*zz
                        ry=(uz*s+ux*uy*(1-c))*xx+(c+uy**2*(1-c))*yy+(-ux*s+uy*uz*(1-c))*zz
                        rz=(-uy*s+ux*uz*(1-c))*xx+(ux*s+uy*uz*(1-c))*yy+(c+uz**2*(1-c))*zz
                        pos((jii-1)*3+1)=x0+rx
                        pos((jii-1)*3+2)=y0+ry
                        pos((jii-1)*3+3)=z0+rz
                     enddo
                  enddo
                  if(iii.eq.4)then
                     kmove=kmove+1
                  endif
               enddo
            enddo
         endif
      enddo
   endif
   
   
 end subroutine microb



 subroutine putbac(key)
   use defs, only : pos,NATOMS,nvar,ovr,svars,var,gra,poss,vars,gras,ovrs
   
   implicit none
   integer ::i,j,ntot,k
   integer,  intent(in) :: key

   ntot=NATOMS*3
  ! c--------------------------------------------------------------save conformation
   if(key.eq.0) then

      poss(1:ntot)=pos(1:ntot)
      vars(1:nvar)=var(1:nvar)
      gras(1:nvar)=gra(1:nvar)
      ovrs(1:nvar)=ovr(1:nvar)
      
   !   do i=1,inb
   !      if(i.eq.ibref)then
   !         do j=1,kuangtot
   !            uangs(i,j,1)=uang(i,j,1)
   !            uangs(i,j,2)=uang(i,j,2)
    !           uangs(i,j,3)=uang(i,j,3)
    !        enddo
    !     else
    !        uangs(i,1,1)=uang(i,1,1)
     !       uangs(i,1,2)=uang(i,1,2)
     !       uangs(i,1,3)=uang(i,1,3)
      !   endif
      !enddo
 !c-----------------------------------------------------------restore conformation
   else
      pos(1:ntot)=poss(1:ntot)
      var(1:nvar)=vars(1:nvar)
      gra(1:nvar)=gras(1:nvar)
      ovr(1:nvar)=ovrs(1:nvar)
      
   !   do i=1,inb
    !     if(i.eq.ibref)then
     !       do j=1,kuangtot
      !         uang(i,j,1)=uangs(i,j,1)
       !        uangs(i,j,2)=uangs(i,j,2)
        !       uangs(i,j,3)=uangs(i,j,3)
         !   enddo
        ! else
        !    uang(i,1,1)=uangs(i,1,1)
         !   uang(i,1,2)=uangs(i,1,2)
         !   uang(i,1,3)=uangs(i,1,3)
         !endif
      !enddo
   endif

  
 end subroutine putbac
 

 
 function torp(i1,i2,i3,i4) 
   use defs                                                        
   implicit none
   integer :: i1,i2,i3,i4    
   real*8 ::  h(3,3),k(3,3),dx1,dx2,dy1,dy2,dz1,dz2
   real*8 :: ux1,ux2,uy1,uy2,uz1,uz2,torp,ru2,ctor,ru1
  
  
   dx1=pos((i2-1)*3+1)-pos((i1-1)*3+1)
   dy1=pos((i2-1)*3+2)-pos((i1-1)*3+2) 
   dz1=pos((i2-1)*3+3)-pos((i1-1)*3+3)
   
   !print*,'pos1',i1,pos((i1-1)*3+1),pos((i1-1)*3+1),pos((i1-1)*3+3)
   !print*,'pos2',i2,pos((i2-1)*3+1),pos((i2-1)*3+2),pos((i2-1)*3+3)

   dx2=pos((i3-1)*3+1)-pos((i2-1)*3+1) 
   dy2=pos((i3-1)*3+2)-pos((i2-1)*3+2) 
   dz2=pos((i3-1)*3+3)-pos((i2-1)*3+3) 

   ux1=dy1*dz2-dz1*dy2 
   uy1=dz1*dx2-dx1*dz2 
   uz1=dx1*dy2-dy1*dx2 

   ru1=sqrt(ux1*ux1+uy1*uy1+uz1*uz1) 
   
   dx1=pos((i4-1)*3+1)-pos((i3-1)*3+1) 
   dy1=pos((i4-1)*3+2)-pos((i3-1)*3+2) 
   dz1=pos((i4-1)*3+3)-pos((i3-1)*3+3) 

   ux2=dz1*dy2-dy1*dz2 
   uy2=dx1*dz2-dz1*dx2 
   uz2=dy1*dx2-dx1*dy2 

   ru2=sqrt(ux2*ux2+uy2*uy2+uz2*uz2) 
   
   ctor=(ux1*ux2+uy1*uy2+uz1*uz2)/(ru1*ru2) 
   if(abs(ctor).gt.1.) ctor=sign(1.d0,ctor) 
   torp=acos(ctor) 
   if(ux1*(uy2*dz2-uz2*dy2)+uy1*(uz2*dx2-ux2*dz2)+uz1*(ux2*dy2-      &
        &  uy2*dx2).lt.0.) torp=-torp                                      
   
 END function torp


 function ang(i1,i2,i3)
   use defs                                                        
   implicit none
   integer :: i1,i2,i3
   real*8 :: dx,dy,dz,cang,rb,bx,by,bz,rd,ang
  
   dx=pos((i1-1)*3+1)-pos((i2-1)*3+1)
   dy=pos((i1-1)*3+2)-pos((i2-1)*3+2)
   dz=pos((i1-1)*3+3)-pos((i2-1)*3+3)
   rd=sqrt(dx*dx+dy*dy+dz*dz)
   bx=pos((i3-1)*3+1)-pos((i2-1)*3+1)
   by=pos((i3-1)*3+2)-pos((i2-1)*3+2)
   bz=pos((i3-1)*3+3)-pos((i2-1)*3+3)
   rb=sqrt(bx*bx+by*by+bz*bz)
   cang = (dx*bx+dy*by+dz*bz)/(rd*rb)
   if(abs(cang).gt.1.) cang=sign(1.d0,cang)
   ang=acos(cang)
  
    end function ang

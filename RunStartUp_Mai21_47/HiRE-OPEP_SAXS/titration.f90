!Br2 : TITRATION Coupling with Fernando's code
subroutine titration(nstep)
  use fileio
  use defs
  use system_defs_RNA
  use rnabase
  use writetofile
  use RNAnb
  use RANDOM
!Br2Fr add new module for titration
  use RNAbaseTit

  character(len=100) :: fname, header
  integer :: i, q
  real :: pch
  real(8) :: ran_number 
  real(8) :: ran3
  

  header = "HEADER"       !Br2
!  call writepdbq("on-the-fly_2.pdb",natoms,pos,header)       
  call writepdb_tit("on-the-fly.pdb",natoms,pos,header)                       !Br2
!  print *, 'on-the-fly pdb written'
!  write(FTIT,'(100i5)') (bprot(i),bocc(i)+5, i =1, 5) 
 
!  write(FTITpc,'(a,100i5)') '..occ', (bocc(i), i =1, NRES)  
!  write(FTIT,'(a,100i5)') 'bfchr', (bprot(i), i =1, NRES)  

  call execute_command_line ("./Titration.sh", wait=.true.)
!  print *, 'titration done'

  open(unit=39,file="charges_from_titration_RNA.dat",status="old")         !Br2
!  print *, 'read charges from titration', NRES, N_RNA
  do i = 1,NRES 
!Br2Fr changing the probability definition
!    write(*,*) fatortit(i)
     read(39,*) pch
!    write(*,*)' Before ',i,pch,chatm(blist(i)),bprot(i),fatortit(i)
!    only bases that can titrate (At, Ct, Ut and Gt) will be evaluated
     if(fatortit(i).gt.0)then
!      call random_number(random)
       ran_number  = ran3()
!      pch has the charge of the phosphates; bprot1 has not!!
!      bprot1 = int(pch+1)
!      aprot converts charges to be numbers between 0 (deprotonated) and 1 (protonated)
       aprot = pch + fatortit(i)
!      write(*,*)' aprot = ',aprot
!      write(*,*)'  random ',ran_number
       if(aprot.ge.ran_number)then
!        protonated for Fernando
!        write(*,*)' protonated '
         bprot1 = 2 - fatortit(i)
         bprot(i) = bprot1*(2-fatortit(i)) 
       else
!        write(*,*)' deprotonated '
         bprot1 =  fatortit(i) - 1
         bprot(i) = bprot1*(2-fatortit(i)) + fatortit(i) - 1
       endif
!      charges *not* considering the phospates!!!
       chatm(blist(i))=bprot1
     endif
!    write(*,*)' After ',i,pch,chatm(blist(i)),bprot(i),fatortit(i)
     bpch(i)=pch      ! charge for the nucleotide, including phosphate
!    bprot(i) = q*q   ! Br2 I can consider only protonated state instead of the charge by taking q*q
!Br2Fr
!    chatm(blist(i))=bprot1
!     print *, i, pch, pch+1, random, q, nint(pch+1)
!    if(i .eq. 25 .and. pch .eq. -1) stop
  enddo
! stop

  close(39)
! stop
!  write(FTIT,'(100i5)')nstep, (bprot(i), i =1, NRES)
  write(FTIT,'(100i5)') (bprot(i), i =1, NRES)  
  write(FTITpc,'(100i5)') (bocc(i), i =1, NRES)  
!  write(FTITpc,'(100f7.2)') (bpch(i), i =1, NRES) 
  write(HBFILE, '(A)') "TITCONF"
end subroutine titration

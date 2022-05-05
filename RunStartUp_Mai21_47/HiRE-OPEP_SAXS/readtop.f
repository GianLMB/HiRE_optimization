C-----------------------------------------------------------------------
c  copy/rewrite from rdtop(nf,chrg) in other-2005-last.f
c  as a helper subroutine
c
c  Xiao Dong, 10/2006
c
c
      SUBROUTINE readtop1(nbondh,nbonda,nbonds)

      ! local renaming of some variables with the [local_name] => [name] syntax,
      ! because of the clashes between names
      use system_defs_prot, only: nbonh_prot => nbonh,
     & nbona_prot => nbona
      use system_defs_RNA, only: nbonh_RNA => nbonh,
     & nbona_RNA => nbona
      use param_defs_prot, only: numbnd_prot => numbnd
      use param_defs_RNA, only: numbnd_RNA => numbnd
      implicit none

      integer, intent(out) :: nbondh, nbonda, nbonds



      nbondh = nbonh_prot + nbonh_RNA
      nbonda = nbona_prot + nbona_RNA
      nbonds = numbnd_prot + numbnd_RNA
      
      RETURN
      END

      subroutine readtop2(nbondt,bia,bib,beq,beq2)
!    $ redu1_mass, redu2_mass)
      
      !use numerical_defs
      use system_defs_prot, only: nbona
      use system_defs_prot, only: nbonh
      use param_defs_prot, only: req, ib, jb, icb
      use param_defs_prot, only: ibh, jbh, icbh
      use system_defs_rna, only: nbona_rna => nbona
      use param_defs_rna, only: req_rna => req,
     % ib_rna => ib, jb_rna => jb, icb_rna => icb
      
      implicit none



      integer :: nbondt
      integer, dimension(nbondt) :: bia, bib
      double precision, dimension(nbondt) :: beq,beq2

      integer :: i, j, i3, j3, ic, itot
      
      itot = 0
      do i = 1, nbonh
        itot = itot + 1
        i3 = ibh(i);  i3 = i3+1;  i3 = i3+2
        j3 = jbh(i);  j3 = j3+1;  j3 = j3+2
        if (mod(i3,3)/=0) stop "i3"
        if (mod(j3,3)/=0) stop "j3"
        i3=i3/3; bia(itot) = i3
        j3=j3/3; bib(itot) = j3
        
        ic = icbh(i)
        beq(itot) = req(ic)
        beq2(itot) = beq(itot)*beq(itot)
      end do


      do i = 1, nbona
        itot = itot + 1
        i3 = ib(i);  i3 = i3+1;  i3 = i3+2
        j3 = jb(i);  j3 = j3+1;  j3 = j3+2
        if (mod(i3,3)/=0) stop "i3"
        if (mod(j3,3)/=0) stop "j3"   
        i3=i3/3; bia(itot) = i3
        j3=j3/3; bib(itot) = j3

        ic = icb(i)
        beq(itot) = req(ic)
        beq2(itot) = beq(itot)*beq(itot)
      end do
      

      do i = 1, nbona_rna
        itot = itot + 1
        i3 = ib_rna(i);  i3 = i3+1;  i3 = i3+2
        j3 = jb_rna(i);  j3 = j3+1;  j3 = j3+2
        if (mod(i3,3)/=0) stop "i3"
        if (mod(j3,3)/=0) stop "j3"   
        i3=i3/3; bia(itot) = i3
        j3=j3/3; bib(itot) = j3

        ic = icb_rna(i)
        beq(itot) = req_rna(ic)
        beq2(itot) = beq(itot)*beq(itot)
      end do
      RETURN
      END


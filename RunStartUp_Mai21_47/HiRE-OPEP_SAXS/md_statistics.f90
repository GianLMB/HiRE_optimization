! This set of routines serves for the initialisation 
!
! Copyright Normand Mousseau January 2006


module md_statistics
  use md_defs
  use geometric_corrections

  implicit none
  save

  contains

  subroutine statistics(iter,status,constrained_energy_in, header_in)
    
    implicit none
    
    logical, intent(in), optional :: header_in
    logical :: header
    
    integer(8), intent(in) :: iter
    character(len=20), intent(in) ::  status
    real(8), intent(in), optional :: constrained_energy_in
    
    integer :: npart, i
    real(8) :: constrained_energy
    real(8) :: delr, rmsd, rmsd_alpha,kinetic, total_runtime,runtime, etot, current_temperature
    real(8), dimension(vecsize) :: posa
    character(len=100) :: fmt_accu
    
    if (present(header_in)) then
       header = header_in
    else
       header = .false.
    endif
    
    if (present(constrained_energy_in)) then
      constrained_energy = constrained_energy_in
    else
      constrained_energy = 0.0d0
    endif
    
    if (status .eq. 'accumulate') then 
       kinetic = 0.0d0
       do i=1, vecsize
          kinetic = kinetic + mass(i)*vel(i)*vel(i)
       end do
       kinetic = kinetic *0.5d0
       current_temperature = 2.0 * kinetic / degree_freedom ! (VECSIZE-6.0)
       
       etot = total_energy+kinetic
       
       posa(:) = pos(:)
       call minimize_rmsd(natoms,posref,posa,delr,rmsd,npart)
       call minimize_rmsd(natoms,posref,posa,delr,rmsd_alpha,npart,.true.,atomic_type)

       avstat = avstat +1
       avpot = avpot + total_energy
       avkin = avkin + kinetic
       avetot = avetot + etot
       avtemp = avtemp + current_temperature
       
       avpot2 = avpot2 + total_energy**2
       avkin2 = avkin2 + kinetic**2
       avetot2 = avetot2 + etot**2
       avtemp2 = avtemp2 + current_temperature**2
       
       runtime=iter*timestep*timeunit/1.0e6
       total_runtime = runtime + initial_simulation_time

       if (constrained_fragments) then
         fmt_accu = "(i14,1x,f11.4,1x,f11.4,7(1x,f12.5))"
         if (debug_status.eq.'debug') write(*,fmt_accu) iter, runtime,total_runtime, &
              constrained_energy, total_energy,kinetic, &
              etot,current_temperature,rmsd_alpha, rmsd
            write(FLOG,fmt_accu) iter, runtime,total_runtime,constrained_energy, &
              total_energy,kinetic,etot, &
              current_temperature,rmsd_alpha,rmsd
       else
        fmt_accu = "(i14,1x,f11.4,1x,f11.4,6(1x,f12.5))"
        if (debug_status.eq.'debug') write(*,fmt_accu) iter, runtime,total_runtime,total_energy,kinetic, &
             etot,current_temperature,rmsd_alpha, rmsd
        write(FLOG,fmt_accu) iter, runtime,total_runtime,total_energy,kinetic,etot, &
             current_temperature,rmsd_alpha,rmsd
       endif
       
    else if (status .eq. 'average') then
       
       if (avstat==0 ) avstat =1
       
       if (avstat .gt. 1 ) then     
          avpot2=sqrt((avpot2-avpot*avpot/avstat)/(avstat-1))
          avkin2=sqrt((avkin2-avkin*avkin/avstat)/(avstat-1))
          avetot2=sqrt((avetot2-avetot*avetot/avstat)/(avstat-1))
          avtemp2=sqrt((avtemp2-avtemp*avtemp/avstat)/(avstat-1))
       else
          avpot2=0.0
          avkin2=0.0
          avetot2=0.0
          avtemp2=0.0
       endif
       
       avpot = avpot/avstat
       avkin = avkin/avstat
       avetot = avetot/avstat
       avtemp = avtemp/avstat
       
       if (debug_status.eq. 'debug') then
          write(*,*) '   Statistics'
          write(*,*) '   **********'
          write(*,*) '   Number of steps  : ', avstat
          write(*,*) '   Potential energy : ', avpot,  ' +/- ', avpot2
          write(*,*) '   Kinetic energy   : ', avkin,  ' +/- ', avkin2
          write(*,*) '   Total energy     : ', avetot, ' +/- ', avetot2
          write(*,*) '   Temperature      : ', avtemp, ' +/- ', avtemp2
          write(*,*) ' ' 
       endif
       
       write(FLOG,*) '   Statistics'
       write(FLOG,*) '   **********'
       write(FLOG,*) '   Number of steps  : ', avstat
       write(FLOG,'(A23,F18.8,A5,F18.8)') '   Potential energy : ', avpot,  ' +/- ', avpot2
       write(FLOG,'(A23,F18.8,A5,F18.8)') '   Kinetic energy   : ', avkin,  ' +/- ', avkin2
       write(FLOG,'(A23,F18.8,A5,F18.8)') '   Total energy     : ', avetot, ' +/- ', avetot2
       write(FLOG,'(A23,F18.8,A5,F18.8)') '   Temperature      : ', avtemp, ' +/- ', avtemp2
       write(FLOG,*) ' ' 
       
    else if (status .eq. 'initialise') then 
       
       ! We initialise the accumulators 
       avstat  = 0
       avpot   = 0.0d0
       avkin   = 0.0d0
       avetot  = 0.0d0
       avtemp  = 0.0d0
       avpot2  = 0.0d0
       avkin2  = 0.0d0
       avetot2  = 0.0d0
       avtemp2 = 0.0d0
       
    endif
    
    if (header) then
      if (constrained_fragments) then
        if (debug_status .eq. 'debug') then
           write(*,'(10A)')  &
                '         Step    Time (ns)  Total time   Const. Ener   Potential     Kinetic',&
                '   Total Energy',&
                '   Temperature    RMSD-CA    RMSD-ALL  '
           write(*,'(10A)')  &
                '  ***********   **********  ***********  ***********  ***********  **********',&
                '  *************',&
                '  ***********   ********** ***********' 
        endif
        write(FLOG,'(10A)')  &
             '         Step    Time (ns)  Total time   Const. Ener  Potential     Kinetic',&
             '   Total Energy   Temperature    RMSD-CA    RMSD-ALL  '
        write(FLOG,'(10A)')  &
             '  ***********   **********  ***********  ***********  ***********  **********',&
             '  *************  ***********   ********** ***********' 
      else 
       if (debug_status .eq. 'debug') then
          write(*,'(10A)')  &
               '         Step    Time (ns)  Total time    Potential     Kinetic   Total Energy',&
               '   Temperature    RMSD-CA    RMSD-ALL  '
          write(*,'(10A)')  &
               '  ***********   **********  ***********  ***********  **********  *************',&
               '  ***********   ********** ***********' 
       endif
       write(FLOG,'(10A)')  &
            '         Step    Time (ns)  Total time    Potential     Kinetic   Total Energy',&
            '   Temperature    RMSD-CA    RMSD-ALL  '
       write(FLOG,'(10A)')  &
            '  ***********   **********  ***********  ***********  **********  *************',&
            '  ***********   ********** ***********' 
      endif
    endif
    
    return
  end subroutine statistics

end module md_statistics

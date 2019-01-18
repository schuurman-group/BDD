  module timingmod

    use constants
    
    implicit none

    save

    real(dp) :: t1,t2,dt

    contains

!#######################################################################

      subroutine get_time(it)

        use channels, only: ilog

        implicit none

        integer :: it

        if (it.eq.1) then
           call cpu_time(t1)
        else if (it.eq.2) then
           call cpu_time(t2)
           dt=t2-t1
        else
           write(ilog,'(/,2x,a,/)') 'Error in get_time'
           STOP
        endif

        return
        
      end subroutine get_time

!#######################################################################

      subroutine times(twall,tcpu)

        use constants

        implicit none

        real(dp) :: twall,tcpu
        integer  :: c,cr,cm

        ! cpu time
        call cpu_time(tcpu)

        ! wall time
        call system_clock(c,cr,cm)
        twall=real(c)/real(cr)
        
        return

      end subroutine times

!#######################################################################

  end module timingmod

module iomod

  save
  
  character(len=200) :: errmsg

contains

!#######################################################################

  subroutine freeunit(unit)

    implicit none

    integer         :: unit,i
    logical(kind=4) :: lopen

    ! N.B. Save the first 20 io units for standard files

    do i=20,1000
       inquire(unit=i,opened=lopen)
       if (.not.lopen) then
          unit=i
          exit
       endif
    enddo

    return

  end subroutine freeunit

!#######################################################################
!
! error_control: writes the passed string to the screen and, if open,
!                the log file, then terminates the program
!
!#######################################################################

  subroutine error_control

    use channels, only: ilog

    implicit none
      
    logical :: lopen
    
    ! Write error message to the screen
    write(6,'(/,2x,a,/)') trim(errmsg)
    
    ! If a log file is open, write the error message to the log file
    inquire(unit=ilog,opened=lopen)
    if (lopen) write(ilog,'(/,2x,a,/)') trim(errmsg)
    
    ! Terminate the program
    STOP
    
  end subroutine error_control
    
!#######################################################################

end module iomod

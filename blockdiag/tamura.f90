!######################################################################
! Projection diabatisation using the algorithm of Tamura. In the limit
! of equal numbers of adiabats and diabats, this method is identical
! to propagative block diagonalisation diabatisation. However, we
! may also use more adiabats than diabats, which might be advantageous
! in some situations
!######################################################################
module tamura

contains

!######################################################################

    subroutine tamura_diabatisation
    
    use constants
    use channels
    use bdglobal
    use adtmod
    
    implicit none

    
    
    return

  end subroutine tamura_diabatisation
    
!######################################################################
  
end module tamura

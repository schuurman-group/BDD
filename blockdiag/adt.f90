!######################################################################
! adt: routines for the calculation of the quasi-diabatic potential
!      matrix
!######################################################################
module adtmod

  implicit none
  
contains

!######################################################################

  subroutine get_adt

    use constants
    use channels
    use iomod
    use utils
    use bdglobal
    
    implicit none

    integer               :: i,j
    real(dp), allocatable :: invspsi(:,:),SST(:,:),sqrtSST(:,:)
    logical               :: pseudo
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(invspsi(nsta,nsta))
    invspsi=0.0d0

    allocate(SST(nsta,nsta))
    SST=0.0d0

    allocate(sqrtSST(nsta,nsta))
    sqrtSST=0.0d0
    
!----------------------------------------------------------------------
! Inversion of the electronic state overlap matrix S
!----------------------------------------------------------------------
    call invert_matrix(spsi,invspsi,nsta,pseudo)

    ! Print out a warning if the pseudoinverse of S is being used
    if (pseudo) write(ilog,'(/,a)') 'WARNING: Pseudoinverse used in &
         the construction of the ADT matrix!'
       
!----------------------------------------------------------------------
! Calculation of (S S^T)^1/2
!----------------------------------------------------------------------
    SST=matmul(spsi,transpose(spsi))
    call sqrt_matrix(SST,sqrtSST,nsta)

!----------------------------------------------------------------------
! Calculation of the ADT matrix: S^-1 (S S^T)^1/2
!----------------------------------------------------------------------
    adt=matmul(invspsi,sqrtSST)
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(invspsi)
    deallocate(SST)
    deallocate(sqrtSST)
    
    return
    
  end subroutine get_adt
  
!######################################################################
  
end module adtmod

!######################################################################
! adt: routines for the calculation of the quasi-diabatic potential
!      matrix
!######################################################################
module adt

  implicit none
  
contains

!######################################################################

  subroutine get_adt

    use constants
    use channels
    use iomod
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

  subroutine invert_matrix(mat,invmat,dim,pseudo)

    use constants
    use channels
    use iomod
    
    implicit none

    integer                      :: dim,info,i
    real(dp), dimension(dim,dim) :: mat,invmat,tmp,umat,vtmat,&
                                    sigmaplus
    real(dp), dimension(dim)     :: sigma
    real(dp), dimension(5*dim)   :: work
    real(dp), parameter          :: thrsh=1e-10_dp
    logical                      :: pseudo
    
!----------------------------------------------------------------------    
! SVD of the input matrix
!----------------------------------------------------------------------
    tmp=mat
    call dgesvd('A','A',dim,dim,tmp,dim,sigma,umat,dim,vtmat,dim,&
         work,5*dim,info)

    ! Exit if the SVD failed
    if (info.ne.0) then
       errmsg='SVD failed in subroutine invert_matrix'
       call error_control
    endif

!----------------------------------------------------------------------
! Moore-Penrose inverse
!----------------------------------------------------------------------
    sigmaplus=0.0d0
    pseudo=.false.
    do i=1,dim
       if (abs(sigma(i)).lt.thrsh) then
          sigmaplus(i,i)=0.0d0
          pseudo=.true.
       else
          sigmaplus(i,i)=1.0d0/sigma(i)
       endif
    enddo

    invmat=matmul(transpose(vtmat),matmul(sigmaplus,transpose(umat)))
    
    return
    
  end subroutine invert_matrix

!######################################################################

  subroutine sqrt_matrix(mat,sqrtmat,dim)

    use constants
    use channels
    use iomod
    
    implicit none

    integer                      :: dim,info,i
    real(dp), dimension(dim,dim) :: mat,sqrtmat,eigvec,dmat
    real(dp), dimension(dim)     :: lambda
    real(dp), dimension(3*dim)   :: work

!----------------------------------------------------------------------
! Diagonalisation of the input matrix
!----------------------------------------------------------------------
    eigvec=mat
    call dsyev('V','U',dim,eigvec,dim,lambda,work,3*dim,info)

    ! Exit if the diagonalisation failed
    if (info.ne.0) then
       errmsg='Diagonalisation failed in subroutine sqrt_matrix'
       call error_control
    endif

!----------------------------------------------------------------------
! Square root of the input matrix
!----------------------------------------------------------------------
    dmat=0.0d0
    do i=1,dim
       dmat(i,i)=sqrt(lambda(i))
    enddo

    sqrtmat=matmul(eigvec,matmul(dmat,transpose(eigvec)))

    return
    
  end subroutine sqrt_matrix
    
!######################################################################
  
end module adt

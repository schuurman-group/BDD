!**********************************************************************
! Routines to perform constant shifts and unitarty transformations of
! the diabatic potential matrix prior to fitting
!**********************************************************************
module transform

contains

!######################################################################
! shift: applies constant shifts to the elements of the diabatic
!        potential matrix
!######################################################################
  subroutine shift

    use constants
    use sysinfo
    use kdcglobal

    implicit none
    
    integer :: n,i,j
    
    ! Loop over geometries
    do n=1,ngeom

       ! Apply the shifts
       do j=1,nsta
          do i=1,nsta
             diabpot(i,j,n)=diabpot(i,j,n)+shift0(i,j)/eh2ev
          enddo
       enddo
       
    enddo
    
    return
    
  end subroutine shift
  
!######################################################################
! blockdiag: block diagonalisation of the diabatic potential matrix
!######################################################################
  subroutine blockdiag

    use constants
    use sysinfo
    use kdcglobal
    use symmetry
    use iomod
    
    implicit none

    integer              :: i,n,nI1,nI2
    integer, allocatable :: I1(:),I2(:)
    integer              :: isI1(nsta)
    real(dp)             :: wmat(nsta,nsta),T(nsta,nsta)
        
!----------------------------------------------------------------------
! Sets of state indices defining the transformation
!----------------------------------------------------------------------
    nI1=nbd(1)
    nI2=nbd(2)
    allocate(I1(nI1), I2(nI2))

    I1=ibd(1:nI1,1)
    I2=ibd(1:nI2,2)

!----------------------------------------------------------------------
! Sanity check 1: block diagonalisation transformations are only
! currently supported for C1 symmetry reference points This is due to
! the ability of the transformation to break state symmetries
!----------------------------------------------------------------------
    if (pntgrp /= 'c1') then
       errmsg='Error: block diagonalisation is only supported in C1'&
            //' symmetry'
    endif

!----------------------------------------------------------------------
! Sanity check 2: make sure that the two sets of state indices are
! disjoint    
!----------------------------------------------------------------------
    isI1=0
    do i=1,nI1
       isI1(I1(i))=1
    enddo

    do i=1,nI2
       if (isI1(I2(i)) == 1) then
          errmsg='Error in blockdiag: the two sets of state indices'&
               //' are not disjoint'
          call error_control
       endif
    enddo
    
!----------------------------------------------------------------------
! Perform the block diagonalisation at all poins
!----------------------------------------------------------------------
    ! Loop over geometries
    do n=1,ngeom

       ! Diabatic potential at this geometry
       wmat=diabpot(:,:,n)

       ! Compute the block diagonalising transformation T
       call block_diag_trans(nsta,wmat,nI1,nI2,I1,I2,T)

       ! Perform the transformation
       diabpot(:,:,n)=matmul(transpose(T),matmul(wmat,T))
       
    enddo

    return
    
  end subroutine blockdiag

!######################################################################
! block_diag_mat: given a matrix W and two sets of indices, I1 and I2,
!                 computes the corresponding block diagonalising
!                 transformation T subject to the constraint
!                 ||T-1|| = min
!######################################################################
  subroutine block_diag_trans(dim,A,nI1,nI2,I1,I2,T)

    use constants
    use utils
    use iomod
    
    implicit none

    ! Input matrix
    integer, intent(in)   :: dim
    real(dp), intent(in)  :: A(dim,dim)

    ! Indices defining the block structure of the transformation
    integer, intent(in)   :: nI1,nI2
    integer, intent(in)   :: I1(nI1),I2(nI2)

    ! Transformation matrix
    real(dp), intent(out) :: T(dim,dim)

    ! Working arrays
    real(dp)              :: work(3*dim)
    real(dp)              :: S(dim,dim),lambda(dim)
    real(dp)              :: S_BD(dim,dim)
    real(dp)              :: U(dim,dim),V(dim,dim)
    real(dp)              :: tmp(dim,dim)
    
    ! Everything else
    integer               :: i,j,ii,jj
    integer               :: error
    real(dp)              :: largest
    
!----------------------------------------------------------------------
! Compute the eigenvectors, S, of the input matrix
!----------------------------------------------------------------------
    S=A
    
    call dsyev('V','U',dim,S,dim,lambda,work,3*dim,error)

    if (error /= 0) then
       errmsg='Diagonalisation failed in subroutine block_diag_trans'
       call error_control
    endif
    
!----------------------------------------------------------------------
! Get the block diagonal part of the matrix of eigenvectors, S_BD
!----------------------------------------------------------------------
    S_BD=0.0d0

    ! First block
    do j=1,nI1
       jj=I1(j)
       do i=1,nI1
          ii=I1(i)
          S_BD(ii,jj)=S(ii,jj)
       enddo
    enddo

    ! Second block
    do j=1,nI2
       jj=I2(j)
       do i=1,nI2
          ii=I2(i)
          S_BD(ii,jj)=S(ii,jj)
       enddo
    enddo

!----------------------------------------------------------------------
! Compute U = S S_BD
!----------------------------------------------------------------------
    U=matmul(S,transpose(S_BD))
    
!----------------------------------------------------------------------
! Compute V=(S_BD S_BD^T)^-1/2
!----------------------------------------------------------------------
    tmp=matmul(S_BD,transpose(S_BD))

    call invsqrt_matrix(tmp,V,dim)

!----------------------------------------------------------------------
! Compute T = U V
!----------------------------------------------------------------------
    T=matmul(U,V)
    
    return
    
  end subroutine block_diag_trans
    
!######################################################################
  
end module transform

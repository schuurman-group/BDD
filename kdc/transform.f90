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
    !if (pntgrp /= 'c1') then
    !   errmsg='Error: block diagonalisation is only supported in C1'&
    !        //' symmetry'
    !   call error_control
    !endif

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
! Allocate the array of block diagonalising transformations
!----------------------------------------------------------------------
    allocate(Tmat(nsta,nsta,ngeom))
    Tmat=0.0d0
    
!----------------------------------------------------------------------
! Perform the block diagonalisation at all points
!----------------------------------------------------------------------
    ! Loop over geometries
    do n=1,ngeom

       ! Diabatic potential at this geometry
       wmat=diabpot(:,:,n)

       ! Compute the block diagonalising transformation T
       call block_diag_trans(nsta,wmat,nI1,nI2,I1,I2,T)

       ! Save the transformation
       Tmat(:,:,n)=T
       
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

!----------------------------------------------------------------------
! Rearrage T s.t. the i'th transformed state is in maximum coincidence
! with the i'th input state
!----------------------------------------------------------------------
    do i=1,dim
       largest=0.0d0
       do j=1,dim
          if (abs(T(i,j)) > largest) then
             largest = abs(T(i,j))
             ii=j
          endif
       enddo
       tmp(:,i)=T(:,ii)
    enddo

    T=tmp
    
    return
    
  end subroutine block_diag_trans
    
!######################################################################
! write_Tmat: writes the block diagonalising transformations to disk
!######################################################################
  subroutine write_Tmat

    use constants
    use sysinfo
    use kdcglobal
    use iomod
    
    implicit none

    integer               :: m,n,s1,s2,ndat
    integer               :: iscratch
    real(dp), allocatable :: q(:),Tij(:)
    character(len=60)     :: filename
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    ! Transformation matrix
    allocate(Tij(maxfiles1m))
    Tij=0.0d0

    ! Normal mode coordinates
    allocate(q(maxfiles1m))
    q=0.0d0
    
!----------------------------------------------------------------------
! Get the next free I/O unit
!----------------------------------------------------------------------
    call freeunit(iscratch)
    
!----------------------------------------------------------------------
! Transformation matrix elements along the 1-mode cuts
!----------------------------------------------------------------------
    ! Loop over modes
    do m=1,nmodes

       ! Cycle if there are no points for the current mode
       ndat=ngeom1m(m)
       if (ndat == 0) cycle

       ! Loop over elements of the transfomation matrix
       do s1=1,nsta
          do s2=1,nsta

             ! Fill in the coordinate transformation matrix element
             ! vectors
             do n=1,ndat
                q(n)=qvec(m,findx1m(m,n))
                Tij(n)=Tmat(s1,s2,findx1m(m,n))
             enddo

             ! Open the output file
             filename=''
             write(filename,'(a,i0,a,i0,a,i0,a)') &
                  'T_',s1,'_',s2,'_q',m,'.dat'
             open(iscratch,file=filename,form='formatted',&
                  status='unknown')

             ! Write the transformation matrix element to disk
             do n=1,ndat
                write(iscratch,'(F10.7,2x,ES11.4)') q(n),Tij(n)
             enddo
             
             ! Close the output file
             close(iscratch)
             
             
          enddo
       enddo
          
    enddo
       
    return
    
  end subroutine write_Tmat
    
!######################################################################
  
end module transform

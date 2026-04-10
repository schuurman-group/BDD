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

    integer              :: i,j,k,n
    integer              :: block_of(nsta)
    real(dp)             :: wmat(nsta,nsta),T(nsta,nsta)

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
! Sanity check 2: make sure that the sets of state indices are
! pairwise disjoint
!----------------------------------------------------------------------
    block_of=0
    do k=1,nblocks
       do j=1,nbd(k)
          i=ibd(j,k)
          if (block_of(i) /= 0) then
             errmsg='Error in blockdiag: the sets of state indices'&
                  //' are not disjoint'
             call error_control
          endif
          block_of(i)=k
       enddo
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
       select case(iblockdiag_alg)
       case(1)
          call block_diag_trans(nsta,wmat,nblocks,nbd,maxbd,ibd,T)
       case(2)
          call block_diag_trans_svd(nsta,wmat,nblocks,nbd,maxbd,ibd,T)
       end select

       ! Save the transformation
       Tmat(:,:,n)=T

       ! Perform the transformation
       diabpot(:,:,n)=matmul(transpose(T),matmul(wmat,T))

    enddo

    return

  end subroutine blockdiag

!######################################################################
! block_diag_trans: given a matrix A and sets of indices defining
!                   N blocks, computes the corresponding block
!                   diagonalising transformation T subject to the
!                   constraint ||T-1|| = min
!                   (Cederbaum, Schirmer, Meyer, J. Phys. A, 22,
!                    2427, 1989, Eq. 7)
!######################################################################
  subroutine block_diag_trans(dim,A,nblks,nbd,maxbd,ibd,T)

    use constants
    use utils
    use iomod

    implicit none

    ! Input matrix
    integer, intent(in)   :: dim
    real(dp), intent(in)  :: A(dim,dim)

    ! Indices defining the block structure of the transformation
    integer, intent(in)   :: nblks
    integer, intent(in)   :: nbd(nblks)
    integer, intent(in)   :: maxbd
    integer, intent(in)   :: ibd(maxbd,nblks)

    ! Transformation matrix
    real(dp), intent(out) :: T(dim,dim)

    ! Working arrays
    real(dp)              :: work(3*dim)
    real(dp)              :: S(dim,dim),lambda(dim)
    real(dp)              :: S_BD(dim,dim)
    real(dp)              :: U(dim,dim),V(dim,dim)
    real(dp)              :: tmp(dim,dim)

    ! Everything else
    integer               :: i,j,k,ii,jj
    integer               :: error

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

    do k=1,nblks
       do j=1,nbd(k)
          jj=ibd(j,k)
          do i=1,nbd(k)
             ii=ibd(i,k)
             S_BD(ii,jj)=S(ii,jj)
          enddo
       enddo
    enddo

!----------------------------------------------------------------------
! Compute U = S S_BD^T
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
! block_diag_trans_svd: SVD-based block diagonalising transformation.
!                       Computes T via the polar decomposition of
!                       each diagonal block S_nn of the eigenvector
!                       matrix, using the SVD to obtain the unitary
!                       factor directly without matrix inversion.
!                       (Cederbaum, Schirmer, Meyer, J. Phys. A, 22,
!                        2427, 1989, Eqs. 11-13)
!######################################################################
  subroutine block_diag_trans_svd(dim,A,nblks,nbd,maxbd,ibd,T)

    use constants
    use iomod

    implicit none

    ! Input matrix
    integer, intent(in)   :: dim
    real(dp), intent(in)  :: A(dim,dim)

    ! Indices defining the block structure of the transformation
    integer, intent(in)   :: nblks
    integer, intent(in)   :: nbd(nblks)
    integer, intent(in)   :: maxbd
    integer, intent(in)   :: ibd(maxbd,nblks)

    ! Transformation matrix
    real(dp), intent(out) :: T(dim,dim)

    ! Eigenvector matrix
    real(dp)              :: S(dim,dim),lambda(dim)
    real(dp)              :: work_eig(3*dim)

    ! SVD working arrays (sized for the largest block)
    integer               :: lwork_svd
    real(dp), allocatable :: Snn(:,:)
    real(dp), allocatable :: U_svd(:,:),VT_svd(:,:)
    real(dp), allocatable :: sigma(:)
    real(dp), allocatable :: work_svd(:)

    ! Block diagonal matrix F
    real(dp)              :: F(dim,dim)

    ! Everything else
    integer               :: i,j,k,ii,jj,dk
    integer               :: error

!----------------------------------------------------------------------
! Compute the eigenvectors, S, of the input matrix
!----------------------------------------------------------------------
    S=A

    call dsyev('V','U',dim,S,dim,lambda,work_eig,3*dim,error)

    if (error /= 0) then
       errmsg='Diagonalisation failed in subroutine '&
            //'block_diag_trans_svd'
       call error_control
    endif

!----------------------------------------------------------------------
! Allocate SVD work arrays for the largest block
!----------------------------------------------------------------------
    lwork_svd=5*maxbd
    allocate(Snn(maxbd,maxbd))
    allocate(U_svd(maxbd,maxbd))
    allocate(VT_svd(maxbd,maxbd))
    allocate(sigma(maxbd))
    allocate(work_svd(lwork_svd))

!----------------------------------------------------------------------
! Construct the block diagonal matrix F from the SVD of each
! diagonal block S_nn of the eigenvector matrix
!----------------------------------------------------------------------
    F=0.0d0

    do k=1,nblks

       dk=nbd(k)

       ! Extract the diagonal block S_nn
       do j=1,dk
          jj=ibd(j,k)
          do i=1,dk
             ii=ibd(i,k)
             Snn(i,j)=S(ii,jj)
          enddo
       enddo

       ! Compute the SVD: S_nn = U sigma V^T
       call dgesvd('A','A',dk,dk,Snn(1:dk,1:dk),dk,sigma(1:dk),&
            U_svd(1:dk,1:dk),dk,VT_svd(1:dk,1:dk),dk,&
            work_svd,lwork_svd,error)

       if (error /= 0) then
          errmsg='SVD failed in subroutine block_diag_trans_svd'
          call error_control
       endif

       ! Compute F_nn = V U^T and place into the full matrix F
       ! dgesvd returns U_svd and VT_svd (= V^T), so
       ! V(i,l) = VT_svd(l,i) and U(j,l) = U_svd(j,l)
       ! F_nn(i,j) = sum_l V(i,l) * U(j,l)
       !           = sum_l VT_svd(l,i) * U_svd(j,l)
       do j=1,dk
          jj=ibd(j,k)
          do i=1,dk
             ii=ibd(i,k)
             F(ii,jj)=dot_product(VT_svd(1:dk,i),U_svd(j,1:dk))
          enddo
       enddo

    enddo

!----------------------------------------------------------------------
! Compute T = S F
!----------------------------------------------------------------------
    T=matmul(S,F)

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(Snn)
    deallocate(U_svd)
    deallocate(VT_svd)
    deallocate(sigma)
    deallocate(work_svd)

    return

  end subroutine block_diag_trans_svd

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

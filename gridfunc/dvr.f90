module dvr

  implicit none

contains

!######################################################################

  subroutine dvr_grids

    use constants
    use iomod
    use sysinfo
    use gridglobal
    
    implicit none

    integer :: m

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    maxgrid=maxval(ndvr)
    allocate(grid(maxgrid,nmodes))
    grid=0.0d0

!----------------------------------------------------------------------    
! Compute the one-mode DVR grids
!----------------------------------------------------------------------
    ! Loop over modes
    do m=1,nmodes

       ! Cycle if a DVR basis was not specified for the current mode
       if (ndvr(m).eq.0) cycle

       ! Compute the DVR grid points for the current mode
       select case(idvr(m))

       case(1) ! HO DVR
          call ho_grid(m)
          
       end select
       
    enddo
    
    return
    
  end subroutine dvr_grids
    
!######################################################################

  subroutine ho_grid(m)

    use constants
    use iomod
    use sysinfo
    use gridglobal
    
    implicit none

    integer, intent(in)   :: m
    integer               :: j,k,j1,k1,ndim
    integer               :: error
    real(dp)              :: x0,omega,mass1
    real(dp), allocatable :: xmat(:,:)
    real(dp), allocatable :: work(:)
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    ! Number of grid points for the m'th mode
    ndim=ndvr(m)

    ! Matrix representation of the position operator
    allocate(xmat(ndim,ndim))
    xmat=0.0d0

    ! Diagonalisation work array
    allocate(work(3*ndim))
    work=0.0d0
    
!----------------------------------------------------------------------
! HO DVR parameters for the m'th mode
!----------------------------------------------------------------------
    ! Equilibrium position
    x0=dvrpar(m,1)
    
    ! Frequency
    omega=dvrpar(m,2)
    
    ! Mass
    mass1=dvrpar(m,3)

!----------------------------------------------------------------------
! Construct the matrix representation of the position operator
! This matrix is tridiagonal, but we loop over all elements because
! I am lazy and it doesn't matter
!----------------------------------------------------------------------
    xmat=0.0d0
    
    do j1=1,ndim
       j=j1-1
       do k1=1,ndim
          k=k1-1

          if (j.eq.k-1) then
             xmat(j1,k1)=(j+1)/(2.0d0*mass1*omega)
             xmat(j1,k1)=sqrt(xmat(j1,k1))
          else if (j.eq.k+1) then
             xmat(j1,k1)=j/(2.0d0*mass1*omega)
             xmat(j1,k1)=sqrt(xmat(j1,k1))
          else if (j.eq.k) then
             xmat(j1,k1)=x0
          endif
          
       enddo
    enddo

!----------------------------------------------------------------------
! Diagonalise the matrix representation of the position operator
! Again, this is a tridiagonal matrix, but its so small that we will
! not bother taking that into account
!----------------------------------------------------------------------
    ! Diagonalisation of xmat
    call dsyev('V','U',ndim,xmat,ndim,grid(1:ndim,m),work,3*ndim,error)
    
    ! Exit if the diagonalisation failed
    if (error.ne.0) then
       errmsg='Diagonalisation failed in ho_grid'
       call error_control
    endif
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(xmat)
    deallocate(work)
    
    return
    
  end subroutine ho_grid
    
!######################################################################
  
end module dvr

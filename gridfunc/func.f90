module func

  implicit none

contains

!######################################################################

  subroutine calc_func

    use constants
    use utils
    use sysinfo
    use gridglobal
    
    implicit none

    integer(8) :: ntotal,i
    integer    :: m,m1
    integer    :: indx(nfuncmode),iswap(nfuncmode)
    real(dp)   :: q(nfuncmode)
    
!----------------------------------------------------------------------
! Sort the direct product (sub) grid modes to be in ascending order
!----------------------------------------------------------------------
    call i4sortindxa1('A',nfuncmode,funcmode,indx)

    do m=1,nfuncmode
       iswap(m)=funcmode(indx(m))
    enddo

    funcmode=iswap

!----------------------------------------------------------------------
! Total number of points on the (sub) direct product grid
!----------------------------------------------------------------------
    ntotal=1
    do m1=1,nfuncmode
       m=funcmode(m1)
       ntotal=ntotal*ndvr(m)
    enddo

!----------------------------------------------------------------------
! Compute the function on the (sub) direct product grid
!----------------------------------------------------------------------
    ! Loop over grid points
    do i=1,ntotal

       ! Mode values at the current grid point
       call mode_values(i,q)

       ! Function value at the current grid point
       select case(ifunc)

          case(1) ! Projector onto an adiabatic state
             call adiabatic_projector(q)
             
       end select
       
    enddo
    
    return
    
  end subroutine calc_func

!######################################################################

  subroutine mode_values(i,q)

    use constants
    use gridglobal
    
    implicit none

    integer(8), intent(in) :: i
    integer(8)             :: j,jj
    integer                :: m1,m
    real(dp), intent(out)  :: q(nfuncmode)

!----------------------------------------------------------------------
! For the i'th point on the (sub) direct product grid, determine the
! mode values. Adapated from the decoder subroutine in the Quantics
! code.
!----------------------------------------------------------------------
    j=i-1
    do m1=1,nfuncmode
       m=funcmode(m1)
       jj=j/ndvr(m)
       q(m1)=j-jj*ndvr(m)+1
       j=jj
    enddo
    
    return
    
  end subroutine mode_values
  
!######################################################################

  subroutine adiabatic_projector(q)

    use constants
    use sysinfo
    use potfuncs
    use gridglobal
    
    implicit none

    integer              :: m,m1
    real(dp), intent(in) :: q(nfuncmode)
    real(dp)             :: q1(nmodes)
    
!----------------------------------------------------------------------
! Compute the ADT matrix at the grid point q
!----------------------------------------------------------------------
    ! Normal mode coordinates in the full space
    q1=0.0d0
    do m1=1,nfuncmode
       m=funcmode(m1)
       q1(m)=q(m1)
    enddo

    ! Get the ADT matrix via the diagonalisation of the model
    ! diabatic potential at the grid point q
    print*,"SORT THIS OUT!"
    stop
    
    
    return
    
  end subroutine adiabatic_projector
    
!######################################################################
  
end module func

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

    integer(8)            :: ntotal,i
    integer               :: m,m1,s1,s2
    integer               :: indx(nfuncmode),iswap(nfuncmode)
    real(dp)              :: q(nfuncmode)
    real(dp), allocatable :: func(:,:,:)
    
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
! Allocate arrays
!----------------------------------------------------------------------
    allocate(func(nsta,nsta,ntotal))
    func=0.0d0
    
!----------------------------------------------------------------------
! Compute the matrix representation of the function on the
! (sub) direct product grid
!----------------------------------------------------------------------
    ! Loop over grid points
    do i=1,ntotal

       ! Mode values at the current grid point
       call mode_values(i,q)

       ! Function value at the current grid point
       select case(ifunc)

          case(1) ! Projector onto an adiabatic state
             func(:,:,i)=adiabatic_projector(q)

          case(2) ! Adiabatic state excitation operator
             func(:,:,i)=adiabatic_excitation(q)
             
       end select
       
    enddo

!----------------------------------------------------------------------
! Write the elements of the matrix representation of the function
! to disk
!----------------------------------------------------------------------
    ! Loop over pairs of states
    do s1=1,nsta
       do s2=s1,nsta

          ! Write the current element of the matrix representation
          ! of the function to disk
          call wrfunc_1element(s1,s2,func(s1,s2,:),ntotal)

       enddo
    enddo
       
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(func)
    
    return
    
  end subroutine calc_func

!######################################################################

  subroutine mode_values(i,q)

    use constants
    use gridglobal
    
    implicit none

    integer(8), intent(in) :: i
    integer(8)             :: j,jj
    integer                :: m1,m,indx
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
       indx=j-jj*ndvr(m)+1
       q(m1)=grid(indx,m)
       j=jj
    enddo
    
    return
    
  end subroutine mode_values
  
!######################################################################

  function adiabatic_projector(q) result(proj)

    use constants
    use sysinfo
    use potfuncs
    use gridglobal
    
    implicit none

    integer              :: m,m1,i,j,iproj
    real(dp), intent(in) :: q(nfuncmode)
    real(dp)             :: proj(nsta,nsta)
    real(dp)             :: q1(nmodes)
    real(dp)             :: adt(nsta,nsta)
        
!----------------------------------------------------------------------
! Normal mode coordinates in the full space
!----------------------------------------------------------------------
    q1=0.0d0
    do m1=1,nfuncmode
       m=funcmode(m1)
       q1(m)=q(m1)
    enddo

!----------------------------------------------------------------------
! Compute the ADT matrix at the grid point q
!----------------------------------------------------------------------    
    adt=adtmatrix(q1)

!----------------------------------------------------------------------
! Diabatic representation of the projector onto the adiabatic state
! of interest
!----------------------------------------------------------------------
    ! Adiabatic state index
    iproj=funcsta(1)

    ! Form the matrix representation of the projector
    do i=1,nsta
       do j=1,nsta
          proj(j,i)=adt(j,iproj)*adt(i,iproj)
       enddo
    enddo
       
    return
    
  end function adiabatic_projector

!######################################################################

  function adiabatic_excitation(q) result(proj)

    use constants
    use sysinfo
    use potfuncs
    use gridglobal
    
    implicit none

    integer              :: m,m1,i,j,iexci1,iexci2
    real(dp), intent(in) :: q(nfuncmode)
    real(dp)             :: proj(nsta,nsta)
    real(dp)             :: q1(nmodes)
    real(dp)             :: adt(nsta,nsta)
        
!----------------------------------------------------------------------
! Normal mode coordinates in the full space
!----------------------------------------------------------------------
    q1=0.0d0
    do m1=1,nfuncmode
       m=funcmode(m1)
       q1(m)=q(m1)
    enddo

!----------------------------------------------------------------------
! Compute the ADT matrix at the grid point q
!----------------------------------------------------------------------    
    adt=adtmatrix(q1)

!----------------------------------------------------------------------
! Diabatic representation of the adiabatic state excitation operator
!----------------------------------------------------------------------
    ! Adiabatic state indices
    iexci1=funcsta(1)
    iexci2=funcsta(2)
    
    ! Form the matrix representation of the excitation operator
    do i=1,nsta
       do j=1,nsta
          proj(j,i)=adt(j,iexci1)*adt(i,iexci2) &
               +adt(i,iexci1)*adt(j,iexci2)
       enddo
    enddo
       
    return
    
  end function adiabatic_excitation
  
!######################################################################

  subroutine wrfunc_1element(s1,s2,func,ntotal)

    use constants
    use iomod
    use sysinfo
    use gridglobal
    
    implicit none

    integer, intent(in)    :: s1,s2
    integer(8), intent(in) :: ntotal
    integer                :: unit
    real(dp), intent(in)   :: func(ntotal)
    character(len=80)      :: filename,stem
    
!----------------------------------------------------------------------
! Open the output file
!----------------------------------------------------------------------
    call freeunit(unit)

    select case(ifunc)
    case(1) ! Projector onto an adiabatic state
       write(stem,'(a,i0)') 'adproj',funcsta(1)
    case(2) ! Adiabatic state excitation operator
       write(stem,'(a,2(i0))') 'adexci',funcsta(1),funcsta(2)
    end select

    write(filename,'(a,i0,a,i0,a)') trim(stem)//'_',s1,'_',s2,'.dat'

    open(unit,file=filename,form='unformatted',status='unknown')

!----------------------------------------------------------------------
! Write the function to file
!----------------------------------------------------------------------
    write(unit) func
    
!----------------------------------------------------------------------
! Close the output file
!----------------------------------------------------------------------
    close(unit)
    
    return
    
  end subroutine wrfunc_1element
    
!######################################################################
  
end module func

module extfunc

  implicit none

contains

!######################################################################

  subroutine getdim_extfunc(i)

    use constants
    use sysinfo
    use iomod
        
    implicit none

    integer, intent(in) :: i
    
    select case(i)

    case(2)
       ! BMA LVC Hamiltonian
       nmodes=2
       nsta=2

    case(3)
       ! Butatriene cation LVC Hamiltonian
       nmodes=2
       nsta=3
       
    case default
       errmsg='Error in nmodes_extfunc: unrecognised function index'
       call error_control
       
    end select
    
    return
    
  end subroutine getdim_extfunc

!######################################################################

  function extfunc_adtmatrix(q) result(adt)

    use constants
    use sysinfo
    use gridglobal
    use iomod
    
    implicit none

    real(dp), intent(in) :: q(nmodes)
    real(dp)             :: adt(nsta,nsta)

    integer              :: e2,error
    real(dp)             :: wmat(nsta,nsta)
    real(dp)             :: v(nsta)
    real(dp)             :: work(3*nsta)
    
!----------------------------------------------------------------------
! Compute the diabatic potential matrix
!----------------------------------------------------------------------
    select case(idiabfunc)

    case(2)
       ! BMA LVC Hamiltonian
       wmat=extfunc_bmalvc(q)

    case(3)
       ! Butatriene cation LVC Hamiltonian
       wmat=extfunc_c4h4lvc(q)
       
    end select

!----------------------------------------------------------------------
! Diagonalise the diabatic potential to get the ADT matrix
!----------------------------------------------------------------------
    e2=3*nsta
    adt=wmat
    call dsyev('V','U',nsta,adt,nsta,v,work,e2,error)

    if (error /= 0) then
       errmsg='error in extfunc_adtmatrix: ' &
            //'diagonalisation of wmat failed'
       call error_control
    endif

    return
    
  end function extfunc_adtmatrix
  
!######################################################################

  function extfunc_bmalvc(q) result(wmat)

    use constants
    use sysinfo
    use gridglobal
    
    implicit none

    ! Function input & output
    real(dp), intent(in) :: q(nmodes)
    real(dp)             :: wmat(nsta,nsta)

    ! Parameters of the BMA LVC model
    real(dp), parameter  :: delta=0.786435d0
    real(dp), parameter  :: omega1=0.210698d0
    real(dp), parameter  :: omega2=0.181772d0
    real(dp), parameter  :: kappa=-0.575674d0
    real(dp), parameter  :: lambda=0.026941d0

    !
    ! Zeroth-order Hamiltonian
    !
    wmat(1,1)=0.5d0*(omega1*q(1)**2+omega2*q(2)**2)

    wmat(2,2)=delta+0.5d0*(omega1*q(1)**2+omega2*q(2)**2)
    
    !
    ! Linear coupling
    !
    wmat(2,2)=wmat(2,2)+kappa*q(1)
    wmat(1,2)=lambda*q(2)
    wmat(2,1)=wmat(1,2)
    
    return
    
  end function extfunc_bmalvc
  
!######################################################################

  function extfunc_c4h4lvc(q) result(wmat)

    use constants
    use sysinfo
    use gridglobal
    
    implicit none

    ! Function input & output
    real(dp), intent(in) :: q(nmodes)
    real(dp)             :: wmat(nsta,nsta)

    ! Parameters of the butatriene cation LVC model
    real(dp), parameter  :: shift=9.0d0
    real(dp), parameter  :: delta=0.019728265d0
    real(dp), parameter  :: omega1=0.26014d0
    real(dp), parameter  :: omega2=0.09116d0
    real(dp), parameter  :: kappa0=-0.1644d0
    real(dp), parameter  :: kappa1=-0.51157432d0
    real(dp), parameter  :: lambda=0.28844084d0

    !
    ! Zeroth-order Hamiltonian
    !
    wmat(1,1)=0.5d0*(omega1*q(1)**2+omega2*q(2)**2)
    wmat(2,2)=-delta+0.5d0*(omega1*q(1)**2+omega2*q(2)**2)
    wmat(3,3)=delta+0.5d0*(omega1*q(1)**2+omega2*q(2)**2)
    
    !
    ! Linear couplings
    !
    wmat(1,1)=wmat(1,1)+kappa0*q(1)
    wmat(2,2)=wmat(2,2)+kappa1*q(1)
    wmat(2,3)=lambda*q(2)
    wmat(3,2)=wmat(2,3)
    
    !
    ! D0 & D1 shift
    !
    wmat(2,2)=wmat(2,2)+shift
    wmat(3,3)=wmat(3,3)+shift

    return
    
  end function extfunc_c4h4lvc

!######################################################################
  
end module extfunc

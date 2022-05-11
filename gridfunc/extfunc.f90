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
    real(dp), parameter  :: omega1=0.210616236d0
    real(dp), parameter  :: omega2=0.181772152d0
    real(dp), parameter  :: kappa1=0.289257182
    real(dp), parameter  :: kappa2=-0.289257182
    real(dp), parameter  :: lambda=0.013469643

    !
    ! Initialisation
    !
    wmat=0.0d0
    
    !
    ! W11
    !
    wmat(1,1)=wmat(1,1)+0.5*omega1*q(1)**2
    wmat(1,1)=wmat(1,1)+0.5*omega2*q(2)**2
    wmat(1,1)=wmat(1,1)+kappa1*q(1)

    !
    ! W22
    !
    wmat(2,2)=wmat(2,2)+0.5*omega1*q(1)**2
    wmat(2,2)=wmat(2,2)+0.5*omega2*q(2)**2
    wmat(2,2)=wmat(2,2)+kappa2*q(1)

    !
    ! W12
    !
    wmat(1,2)=wmat(1,2)+lambda*q(2)
    wmat(2,1)=wmat(1,2)
    
    return
    
  end function extfunc_bmalvc
  
!######################################################################
  
end module extfunc

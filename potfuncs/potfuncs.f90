module potfuncs

  implicit none

contains

!######################################################################

  function adiabaticpot(q) result(v)

    use constants
    use iomod
    use sysinfo
    
    implicit none

    integer                        :: e2,error
    real(dp), dimension(nmodes)    :: q
    real(dp), dimension(nsta)      :: v
    real(dp), dimension(nsta,nsta) :: w
    real(dp), dimension(3*nsta)    :: work

!----------------------------------------------------------------------
! Construct the model potential
!----------------------------------------------------------------------
    w=pot(q)

!----------------------------------------------------------------------
! Diagonalise the model potential to yield the model adiabatic
! potentials
!----------------------------------------------------------------------
    e2=3*nsta
    call dsyev('V','U',nsta,w,nsta,v,work,e2,error)

    if (error.ne.0) then
       write(6,'(/,2x,a,/)') 'Diagonalisation of the potential &
            matrix failed'
       stop
    endif
    
    return

  end function adiabaticpot

!######################################################################
  
  function diabaticpot(q) result(w)

    use constants
    use sysinfo

    implicit none

    integer                        :: i
    real(dp), dimension(nmodes)    :: q
    real(dp), dimension(nsta,nsta) :: w

!----------------------------------------------------------------------
! Construct the model potential
!----------------------------------------------------------------------
    w=pot(q)

    return

  end function diabaticpot

!######################################################################

  function adtmatrix(q) result(adt)

    use constants
    use iomod
    use sysinfo
    
    implicit none

    integer                        :: e2,error
    real(dp), dimension(nmodes)    :: q
    real(dp), dimension(nsta)      :: v
    real(dp), dimension(nsta,nsta) :: w
    real(dp), dimension(nsta,nsta) :: adt
    real(dp), dimension(3*nsta)    :: work
    
!----------------------------------------------------------------------
! Construct the model potential
!----------------------------------------------------------------------
    w=pot(q)

!----------------------------------------------------------------------
! Diagonalise the model potential
!----------------------------------------------------------------------
    e2=3*nsta
    call dsyev('V','U',nsta,w,nsta,v,work,e2,error)

    if (error.ne.0) then
       write(6,'(/,2x,a,/)') 'Diagonalisation of the potential &
            matrix failed'
       stop
    endif

!----------------------------------------------------------------------
! Return the eigenvectors of the model diabatic potential matrix
!----------------------------------------------------------------------
    adt=w
    
    return
    
  end function adtmatrix
    
!######################################################################

  function pot(q) result(w)

    use constants
    use sysinfo
    use symmetry
    use parameters

    implicit none

    integer                        :: m,m1,m2,s,s1,s2,n
    real(dp), dimension(nmodes)    :: q
    real(dp), dimension(nsta,nsta) :: w
    real(dp)                       :: fac,pre
    
!----------------------------------------------------------------------
! Initialisation of the model potential
!----------------------------------------------------------------------
    w=0.0d0

!----------------------------------------------------------------------
! Zeroth-order contributions
!----------------------------------------------------------------------
    ! Energies
    do s=1,nsta
       w(s,s)=w(s,s)+e0(s)
    enddo
    
    ! Harmonic potentials
    do s=1,nsta
       do m=1,nmodes
          w(s,s)=w(s,s)+0.5d0*freq(m)*q(m)**2
       enddo
    enddo

!----------------------------------------------------------------------
! One-mode contributions
!----------------------------------------------------------------------
    fac=1.0d0

    do n=1,order1

       fac=fac*n
       pre=1.0d0/fac
       
       do s2=1,nsta
          do s1=1,nsta
             do m=1,nmodes
                w(s1,s2)=w(s1,s2)+pre*coeff1(m,s1,s2,n)*q(m)**n
             enddo
          enddo
       enddo
       
    enddo

!----------------------------------------------------------------------
! Two-mode contributions (2nd-order only)
!----------------------------------------------------------------------
    do s2=1,nsta
       do s1=1,nsta
          do m2=1,nmodes
             do m1=1,nmodes
                w(s1,s2)=w(s1,s2)+0.5d0*coeff2(m1,m2,s1,s2)*q(m1)*q(m2)
             enddo
          enddo
       enddo
    enddo

!!----------------------------------------------------------------------
!! First-order contributions
!!----------------------------------------------------------------------
!    ! kappa
!    do s=1,nsta
!       do m=1,nmodes
!          if (kappa_mask(m,s).eq.0) cycle
!          w(s,s)=w(s,s)+kappa(m,s)*q(m)
!       enddo
!    enddo
!
!    ! lambda
!    do s1=1,nsta-1
!       do s2=s1+1,nsta
!          do m=1,nmodes
!             if (lambda_mask(m,s1,s2).eq.0) cycle
!             w(s1,s2)=w(s1,s2)+lambda(m,s1,s2)*q(m)
!             w(s2,s1)=w(s2,s1)+lambda(m,s1,s2)*q(m)
!          enddo
!       enddo
!    enddo
!
!!----------------------------------------------------------------------
!! Second-order contributions
!!----------------------------------------------------------------------
!    ! gamma
!    do s=1,nsta
!       do m1=1,nmodes
!          do m2=1,nmodes
!             if (gamma_mask(m1,m2,s).eq.0) cycle
!             w(s,s)=w(s,s)+0.5d0*gamma(m1,m2,s)*q(m1)*q(m2)
!          enddo
!       enddo
!    enddo
!
!    ! mu
!    do s1=1,nsta-1
!       do s2=s1+1,nsta
!          do m1=1,nmodes
!             do m2=1,nmodes
!                if (mu_mask(m1,m2,s1,s2).eq.0) cycle
!                w(s1,s2)=w(s1,s2)+0.5d0*mu(m1,m2,s1,s2)*q(m1)*q(m2)
!                w(s2,s1)=w(s2,s1)+0.5d0*mu(m1,m2,s1,s2)*q(m1)*q(m2)
!             enddo
!          enddo
!       enddo
!    enddo
!
!!----------------------------------------------------------------------
!! Third-order contributions
!!----------------------------------------------------------------------
!    ! iota
!    do s=1,nsta
!       do m=1,nmodes
!          if (iota_mask(m,s).eq.0) cycle
!          w(s,s)=w(s,s)+(1.0d0/6.0d0)*iota(m,s)*q(m)**3
!       enddo
!    enddo
!
!    ! tau
!    do s1=1,nsta-1
!       do s2=s1+1,nsta
!          do m=1,nmodes
!             if (tau_mask(m,s1,s2).eq.0) cycle
!             w(s1,s2)=w(s1,s2)+(1.0d0/6.0d0)*tau(m,s1,s2)*q(m)**3
!             w(s2,s1)=w(s2,s1)+(1.0d0/6.0d0)*tau(m,s1,s2)*q(m)**3
!          enddo
!       enddo
!    enddo
!
!!----------------------------------------------------------------------
!! Fourth-order contributions
!!----------------------------------------------------------------------
!    ! epsilon
!    do s=1,nsta
!       do m=1,nmodes
!          if (epsilon_mask(m,s).eq.0) cycle
!          w(s,s)=w(s,s)+(1.0d0/24.0d0)*epsilon(m,s)*q(m)**4
!       enddo
!    enddo
!
!    ! xi
!    do s1=1,nsta-1
!       do s2=s1+1,nsta
!          do m=1,nmodes
!             if (xi_mask(m,s1,s2).eq.0) cycle
!             w(s1,s2)=w(s1,s2)+(1.0d0/24.0d0)*xi(m,s1,s2)*q(m)**4
!             w(s2,s1)=w(s2,s1)+(1.0d0/24.0d0)*xi(m,s1,s2)*q(m)**4
!          enddo
!       enddo
!    enddo
    
    return

  end function pot

!######################################################################

  function diabdip(q,s1,s2) result(dip)

    use constants
    use sysinfo
    use symmetry
    use parameters

    implicit none

    integer, intent(in)                     :: s1,s2
    integer                                 :: m,m1,m2,c
    real(dp), dimension(nmodes), intent(in) :: q
    real(dp), dimension(3)                  :: dip

!----------------------------------------------------------------------
! Zeroth-order contribution
!----------------------------------------------------------------------
    dip(:)=dip0(s1,s2,:)

!----------------------------------------------------------------------
! First-order contributions
!----------------------------------------------------------------------
    do c=1,3
       do m=1,nmodes
          if (dip1_mask(m,s1,s2,c).eq.0) cycle
          dip(c)=dip(c)+dip1(m,s1,s2,c)*q(m)
       enddo
    enddo

!----------------------------------------------------------------------
! Second-order contributions
!----------------------------------------------------------------------
    do c=1,3
       do m1=1,nmodes
          do m2=1,nmodes
             if (dip2_mask(m1,m2,s1,s2,c).eq.0) cycle
             dip(c)=dip(c)+0.5d0*dip2(m1,m2,s1,s2,c)*q(m1)*q(m2)
          enddo
       enddo
    enddo

!----------------------------------------------------------------------
! Third-order contributions
!----------------------------------------------------------------------
    do c=1,3
       do m=1,nmodes
          if (dip3_mask(m,s1,s2,c).eq.0) cycle
          dip(c)=dip(c)+(1.0d0/6.0d0)*dip3(m,s1,s2,c)*q(m)**3
       enddo
    enddo

!----------------------------------------------------------------------
! Fourth-order contributions
!----------------------------------------------------------------------
    do c=1,3
       do m=1,nmodes
          if (dip4_mask(m,s1,s2,c).eq.0) cycle
          dip(c)=dip(c)+(1.0d0/24.0d0)*dip4(m,s1,s2,c)*q(m)**4
       enddo
    enddo
    
    return
    
  end function diabdip
  
!######################################################################

end module potfuncs

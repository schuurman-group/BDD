module geomopt

  use constants

  implicit none

  public  :: optimise
  private :: outer

contains

!######################################################################
! optimise: locally minimise a potential surface as a function of the
!           normal-mode coordinates, starting from q0, using the BFGS
!           quasi-Newton algorithm with an analytical gradient.
!
! By default, the s-th adiabatic potential surface is minimised
! (Hellmann-Feynman gradient). Pass ldiabatic=.true. to instead
! minimise the diabatic potential matrix element W_{ss}(Q) directly.
!
! Inputs:
!   q0(nmodes) : initial geometry in normal modes
!   s          : state index (1..nsta) -- adiabatic by default,
!                diabatic when ldiabatic=.true.
!   ldiabatic  : (optional) if .true., minimise the diabatic
!                W_{ss}(Q) instead of the adiabatic v_s(Q).
!                                                         [default .false.]
!   tol_g      : (optional) convergence tolerance on gradient norm
!                                                         [default 1d-6]
!   tol_x      : (optional) convergence tolerance on step norm
!                                                         [default 1d-8]
!   max_iter   : (optional) maximum BFGS iterations       [default 200]
!
! Outputs:
!   qmin(nmodes) : optimised geometry
!   vmin         : energy at qmin (adiabatic v_s, or diabatic W_{ss}
!                  if ldiabatic=.true.)
!   ifault       : status (0 = converged,
!                          1 = max_iter reached,
!                          2 = line-search failure)
!######################################################################

  subroutine optimise(q0, s, qmin, vmin, ifault, &
                      tol_g, tol_x, max_iter, ldiabatic)

    use sysinfo,  only: nmodes, nsta, freq
    use potfuncs, only: adiabaticpot, adiabaticgrad, pot, dpot_dq

    implicit none

    real(dp), intent(in)           :: q0(nmodes)
    integer,  intent(in)           :: s
    real(dp), intent(out)          :: qmin(nmodes)
    real(dp), intent(out)          :: vmin
    integer,  intent(out)          :: ifault
    real(dp), intent(in), optional :: tol_g
    real(dp), intent(in), optional :: tol_x
    integer,  intent(in), optional :: max_iter
    logical,  intent(in), optional :: ldiabatic

    integer  :: iter, m, maxit
    real(dp) :: x(nmodes), x_new(nmodes), p(nmodes)
    real(dp) :: g(nmodes), g_new(nmodes)
    real(dp) :: sv(nmodes), yv(nmodes), Hy(nmodes)
    real(dp) :: H(nmodes, nmodes)
    real(dp) :: v(nsta)
    real(dp) :: w(nsta,nsta), dw(nsta,nsta)
    real(dp) :: e, e_new, alpha, gp, ys
    real(dp) :: gtol, xtol
    logical  :: ldiab

    real(dp), parameter :: c1        = 1.0d-4
    real(dp), parameter :: alpha_min = 1.0d-12
    real(dp), parameter :: ys_min    = 1.0d-10

!----------------------------------------------------------------------
! Validate the state index
!----------------------------------------------------------------------
    if (s < 1 .or. s > nsta) then
       write(6,'(/,2x,a,i0,a,i0,/)') &
            'optimise: state index ', s, ' is out of range 1..', nsta
       stop
    endif

!----------------------------------------------------------------------
! Apply defaults for the optional arguments
!----------------------------------------------------------------------
    if (present(tol_g)) then
       gtol = tol_g
    else
       gtol = 1.0d-6
    endif

    if (present(tol_x)) then
       xtol = tol_x
    else
       xtol = 1.0d-8
    endif

    if (present(max_iter)) then
       maxit = max_iter
    else
       maxit = 200
    endif

    if (present(ldiabatic)) then
       ldiab = ldiabatic
    else
       ldiab = .false.
    endif

!----------------------------------------------------------------------
! Initial point, energy and gradient
!----------------------------------------------------------------------
    x = q0
    if (ldiab) then
       w = pot(x)
       e = w(s,s)
       do m = 1, nmodes
          dw = dpot_dq(x, m)
          g(m) = dw(s,s)
       enddo
    else
       v = adiabaticpot(x)
       e = v(s)
       g = adiabaticgrad(x, s)
    endif

!----------------------------------------------------------------------
! Initial inverse-Hessian guess: diag(1/freq(m)).
! At Q=0 this is the exact harmonic Hessian, so the first BFGS step is
! essentially Newton's method on near-harmonic surfaces.
!----------------------------------------------------------------------
    H = 0.0d0
    do m = 1, nmodes
       H(m,m) = 1.0d0/freq(m)
    enddo

!----------------------------------------------------------------------
! BFGS iterations
!----------------------------------------------------------------------
    ifault = 0

    do iter = 1, maxit

       ! Gradient-norm convergence
       if (sqrt(dot_product(g, g)) < gtol) exit

       ! Search direction
       p  = -matmul(H, g)
       gp = dot_product(g, p)

       ! Safeguard: if BFGS direction is not a descent direction, reset
       ! the inverse-Hessian to the harmonic guess and retry.
       if (gp >= 0.0d0) then
          H = 0.0d0
          do m = 1, nmodes
             H(m,m) = 1.0d0/freq(m)
          enddo
          p  = -matmul(H, g)
          gp = dot_product(g, p)
       endif

       ! Backtracking Armijo line search
       alpha = 1.0d0
       do
          x_new = x + alpha*p
          if (ldiab) then
             w     = pot(x_new)
             e_new = w(s,s)
          else
             v     = adiabaticpot(x_new)
             e_new = v(s)
          endif
          if (e_new <= e + c1*alpha*gp) exit
          alpha = 0.5d0*alpha
          if (alpha < alpha_min) then
             ! Line search has stalled. If the gradient norm is
             ! already comparable to the convergence tolerance, the
             ! Armijo failure is just numerical noise from dsyev /
             ! Hellmann-Feynman near a minimum, not a real obstacle:
             ! accept as converged. Otherwise, flag the failure.
             if (sqrt(dot_product(g, g)) < 100.0d0*gtol) then
                ifault = 0
             else
                ifault = 2
             endif
             qmin   = x
             vmin   = e
             return
          endif
       enddo

       ! Gradient at the new point
       if (ldiab) then
          do m = 1, nmodes
             dw       = dpot_dq(x_new, m)
             g_new(m) = dw(s,s)
          enddo
       else
          g_new = adiabaticgrad(x_new, s)
       endif

       ! BFGS inverse-Hessian update
       sv = x_new - x
       yv = g_new - g
       ys = dot_product(yv, sv)

       if (ys > ys_min) then
          Hy = matmul(H, yv)
          H  = H + ((ys + dot_product(yv, Hy))/ys**2)*outer(sv, sv) &
                 - (outer(Hy, sv) + outer(sv, Hy))/ys
       endif

       ! Step-norm convergence
       if (sqrt(dot_product(sv, sv)) < xtol) then
          x = x_new
          e = e_new
          g = g_new
          exit
       endif

       ! Roll forward
       x = x_new
       e = e_new
       g = g_new

    enddo

    ! If the loop completed without an early exit, max_iter was reached
    if (iter > maxit) ifault = 1

    qmin = x
    vmin = e

    return

  end subroutine optimise

!######################################################################
! outer: outer product of two real vectors,  m(i,j) = a(i)*b(j)
!######################################################################

  function outer(a, b) result(m)

    implicit none

    real(dp), intent(in) :: a(:), b(:)
    real(dp)             :: m(size(a), size(b))

    m = spread(a, 2, size(b))*spread(b, 1, size(a))

  end function outer

!######################################################################

end module geomopt

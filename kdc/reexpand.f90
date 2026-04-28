!**********************************************************************
! Re-expansion of the fitted diabatic potential matrix elements about
! the minimum of an adiabatic surface (Qmin), in normal-mode
! coordinates. Activated by the $reexpand keyword in the input file.
!
! The KDC model is a polynomial of total degree max(order1, 2) in Q.
! Substituting Q -> Qmin + dQ produces another polynomial of the same
! degree with the same monomial basis, so the shift is exact and
! analytical (no re-fitting). The bilinear coefficients coeff2 are
! unchanged. The 1-mode coefficients coeff1 transform via the Taylor
! shift formula. The new diagonal energies are written into e0(:);
! the new off-diagonal zeroth-order constants (which were zero in the
! original model) are written into e0_off(:,:).
!
! The on-diagonal harmonic baseline 0.5*freq(m)*Q_m^2 is folded into
! the effective second-order coefficient before the shift, then the
! same freq(m) is subtracted back from the new order-2 coefficient on
! the diagonal so that the storage convention coeff1(m,s,s,2) =
! (true second-derivative) - freq(m) (set in nmeqmod.f90) is preserved.
!**********************************************************************
module reexpand

  use constants

  implicit none

  public  :: reexpand_about_qmin
  private :: factorial, taylor_shift_diag, taylor_shift_offdiag, &
             write_reexpand_log, reexpand_diagnostics

contains

!######################################################################
! reexpand_about_qmin: top-level driver. Locates Qmin via the BFGS
! optimiser in geomopt, then re-expands all diabatic potential matrix
! elements about Qmin, updating e0, coeff1, and e0_off in place
! (coeff2 is invariant under the shift).
!######################################################################
  subroutine reexpand_about_qmin

    use channels
    use sysinfo,    only: nmodes, nsta, freq
    use parameters, only: e0, e0_off, coeff1, coeff2, order1
    use kdcglobal,  only: lreexpand, ireexpand_state, qmin_re
    use geomopt,    only: optimise

    implicit none

    integer  :: s, s1, s2, ifault
    real(dp) :: q0(nmodes), qmin(nmodes), vmin
    real(dp) :: new_e0(nsta,nsta)
    real(dp) :: new_coeff1(nmodes,nsta,nsta,order1)

!----------------------------------------------------------------------
! Locate Qmin: BFGS minimisation of the adiabatic potential v_s(Q) on
! the requested state, starting from the reference geometry Q=0.
!----------------------------------------------------------------------
    q0 = 0.0d0

    call optimise(q0, ireexpand_state, qmin, vmin=vmin, &
                  ifault=ifault, tol_g=1.0d-8)

    if (ifault /= 0) then
       write(ilog,'(/,72a)') ('+',s=1,72)
       write(ilog,'(2x,a)') 'Re-expansion of the diabatic potential'
       write(ilog,'(72a)') ('+',s=1,72)
       write(ilog,'(/,2x,a,i0,a)') &
            'WARNING: BFGS optimisation did not converge (ifault = ', &
            ifault,'); the diabatic potential will NOT be re-expanded.'
       lreexpand = .false.
       return
    endif

!----------------------------------------------------------------------
! Save Qmin so it can be persisted in the binary file (qmin_re is
! allocated and zeroed earlier in get_coefficients)
!----------------------------------------------------------------------
    qmin_re = qmin

!----------------------------------------------------------------------
! Diagnostics: print items 1-6 (gradient norms, residual coupling,
! eigenvalues, max off-diagonal coeff1/coeff2, adiabatic vs diabatic
! Qmin comparison) before the Taylor shift overwrites the polynomial
!----------------------------------------------------------------------
    call reexpand_diagnostics(qmin, ireexpand_state, .false.)

!----------------------------------------------------------------------
! Compute the new coefficients for each diabatic matrix element
!----------------------------------------------------------------------
    new_e0     = 0.0d0
    new_coeff1 = 0.0d0

    do s = 1, nsta
       call taylor_shift_diag(s, qmin, freq, e0(s), &
                              coeff1(:,s,s,:), coeff2(:,:,s,s), &
                              new_e0(s,s), new_coeff1(:,s,s,:))
    enddo

    do s1 = 1, nsta-1
       do s2 = s1+1, nsta
          call taylor_shift_offdiag(qmin, &
                                    coeff1(:,s1,s2,:), coeff2(:,:,s1,s2), &
                                    new_e0(s1,s2), new_coeff1(:,s1,s2,:))
          ! Symmetrise: coeff1 and e0_off are symmetric in (s1,s2)
          new_e0(s2,s1)         = new_e0(s1,s2)
          new_coeff1(:,s2,s1,:) = new_coeff1(:,s1,s2,:)
       enddo
    enddo

!----------------------------------------------------------------------
! Commit the new coefficients
!----------------------------------------------------------------------
    do s = 1, nsta
       e0(s) = new_e0(s,s)
    enddo

    do s1 = 1, nsta
       do s2 = 1, nsta
          if (s1 == s2) then
             e0_off(s1,s2) = 0.0d0
          else
             e0_off(s1,s2) = new_e0(s1,s2)
          endif
       enddo
    enddo

    coeff1 = new_coeff1

!----------------------------------------------------------------------
! Log the result
!----------------------------------------------------------------------
    call write_reexpand_log(qmin, vmin, ifault)

!----------------------------------------------------------------------
! Diagnostics: item 7 (max |e0_off|) after the Taylor shift
!----------------------------------------------------------------------
    call reexpand_diagnostics(qmin, ireexpand_state, .true.)

    return

  end subroutine reexpand_about_qmin

!######################################################################
! taylor_shift_diag: compute new diagonal Taylor coefficients for a
! single diabatic matrix element W_{ss} after the shift Q -> Qmin + dQ.
!
! Inputs:
!   s_idx   : state index (used only for clarity in messages)
!   q       : Qmin(nmodes)
!   freq    : harmonic frequencies in normal modes (eV)
!   e0_s    : current vertical excitation energy e0(s) (eV)
!   c1_s    : current 1-mode coefficients coeff1(m,s,s,1..order1)
!   c2_s    : current bilinear coefficients coeff2(m1,m2,s,s)
!             (symmetric storage; diagonal m1=m2 may also be non-zero)
!
! Outputs:
!   new_e0    : new W_{ss}(Qmin) - i.e. the new diagonal energy
!   new_c1    : new 1-mode coefficients (with the freq(m) baseline
!               subtracted back from order 2 to preserve the model's
!               storage convention)
!
! Algorithm (per mode m):
!   a_n = c1_s(m,n) for n != 2
!   a_2 = freq(m) + c1_s(m,2)              (fold harmonic in)
!   new_a_k(m) = sum_{n=k..order1} a_n / (n-k)! * q(m)^(n-k), k=1..order1
!   new_a_1(m) += sum_{m' != m} c2_s(m,m') * q(m')          (bilinear)
!   new_e0_contribution from mode m: sum_{n=1..order1} a_n / n! * q(m)^n
! Bilinear contribution to constant:
!   new_e0 += 0.5 * sum_{m1,m2} c2_s(m1,m2) * q(m1) * q(m2)
! Restore harmonic baseline at order 2:
!   new_c1(m,2) = new_a_2(m) - freq(m)
!######################################################################
  subroutine taylor_shift_diag(s_idx, q, freq, e0_s, c1_s, c2_s, &
                               new_e0, new_c1)

    use sysinfo,    only: nmodes
    use parameters, only: order1

    implicit none

    integer,  intent(in)  :: s_idx
    real(dp), intent(in)  :: q(nmodes)
    real(dp), intent(in)  :: freq(nmodes)
    real(dp), intent(in)  :: e0_s
    real(dp), intent(in)  :: c1_s(nmodes,order1)
    real(dp), intent(in)  :: c2_s(nmodes,nmodes)
    real(dp), intent(out) :: new_e0
    real(dp), intent(out) :: new_c1(nmodes,order1)

    integer  :: m, mp, n, k
    real(dp) :: a(order1), new_a(order1)
    real(dp) :: per_mode_const

    new_e0 = e0_s
    new_c1 = 0.0d0

    ! Loop over modes
    do m = 1, nmodes

       ! Build the effective per-mode Taylor coefficients (a_n) by
       ! folding the harmonic baseline freq(m) into a_2
       do n = 1, order1
          a(n) = c1_s(m,n)
       enddo
       if (order1 >= 2) a(2) = a(2) + freq(m)

       ! Apply the 1-mode Taylor shift:
       ! new_a_k = sum_{n=k..order1} a_n / (n-k)! * q(m)^(n-k)
       do k = 1, order1
          new_a(k) = 0.0d0
          do n = k, order1
             new_a(k) = new_a(k) &
                      + a(n) / real(factorial(n-k), dp) * q(m)**(n-k)
          enddo
       enddo

       ! Per-mode contribution to the new constant term (k=0):
       ! sum_{n=1..order1} a_n / n! * q(m)^n
       per_mode_const = 0.0d0
       do n = 1, order1
          per_mode_const = per_mode_const &
                         + a(n) / real(factorial(n), dp) * q(m)**n
       enddo
       new_e0 = new_e0 + per_mode_const

       ! Store the new 1-mode coefficients for this mode, restoring
       ! the storage convention (subtract freq(m) back at order 2)
       do k = 1, order1
          new_c1(m,k) = new_a(k)
       enddo
       if (order1 >= 2) new_c1(m,2) = new_c1(m,2) - freq(m)

    enddo

    ! Bilinear contribution to the constant:
    !   0.5 * sum_{m1,m2} c2(m1,m2) * q(m1) * q(m2)
    do m = 1, nmodes
       do mp = 1, nmodes
          new_e0 = new_e0 + 0.5d0 * c2_s(m,mp) * q(m) * q(mp)
       enddo
    enddo

    ! Bilinear contribution to the new linear (k=1) coefficient of
    ! mode m:
    !   sum_{mp != m} c2(m,mp) * q(mp)
    ! (the factor of 1/2 in pot() is cancelled by the symmetric storage
    ! c2(m,mp) = c2(mp,m): the (m,mp) and (mp,m) pair contributes
    ! c2(m,mp)*q(mp) + c2(mp,m)*q(mp) = 2*c2(m,mp)*q(mp), times 0.5.)
    !
    ! NB: the m=mp diagonal of coeff2 is assumed to be zero (the
    ! standard fit only populates m1<m2 pairs). If a future caller
    ! sets c2(m,m) /= 0, the order-1 and order-2 contributions from
    ! the m=mp term would also need to be folded in here.
    do m = 1, nmodes
       do mp = 1, nmodes
          if (mp == m) cycle
          new_c1(m,1) = new_c1(m,1) + c2_s(m,mp) * q(mp)
       enddo
    enddo

    return

  end subroutine taylor_shift_diag

!######################################################################
! taylor_shift_offdiag: compute new Taylor coefficients for a single
! off-diagonal matrix element W_{s1,s2} (s1 != s2) after the shift
! Q -> Qmin + dQ. There is no harmonic baseline on off-diagonals, so
! a_n is just c1(m,n).
!
! Inputs / outputs as in taylor_shift_diag, except no e0_s and no
! freq baseline. new_e0 here is the new off-diagonal constant
! (W_{s1,s2}(Qmin)).
!######################################################################
  subroutine taylor_shift_offdiag(q, c1_ij, c2_ij, new_e0, new_c1)

    use sysinfo,    only: nmodes
    use parameters, only: order1

    implicit none

    real(dp), intent(in)  :: q(nmodes)
    real(dp), intent(in)  :: c1_ij(nmodes,order1)
    real(dp), intent(in)  :: c2_ij(nmodes,nmodes)
    real(dp), intent(out) :: new_e0
    real(dp), intent(out) :: new_c1(nmodes,order1)

    integer  :: m, mp, n, k
    real(dp) :: a(order1), new_a(order1)
    real(dp) :: per_mode_const

    new_e0 = 0.0d0
    new_c1 = 0.0d0

    do m = 1, nmodes

       do n = 1, order1
          a(n) = c1_ij(m,n)
       enddo

       do k = 1, order1
          new_a(k) = 0.0d0
          do n = k, order1
             new_a(k) = new_a(k) &
                      + a(n) / real(factorial(n-k), dp) * q(m)**(n-k)
          enddo
       enddo

       per_mode_const = 0.0d0
       do n = 1, order1
          per_mode_const = per_mode_const &
                         + a(n) / real(factorial(n), dp) * q(m)**n
       enddo
       new_e0 = new_e0 + per_mode_const

       do k = 1, order1
          new_c1(m,k) = new_a(k)
       enddo

    enddo

    ! Bilinear contribution to the constant
    do m = 1, nmodes
       do mp = 1, nmodes
          new_e0 = new_e0 + 0.5d0 * c2_ij(m,mp) * q(m) * q(mp)
       enddo
    enddo

    ! Bilinear contribution to the new linear coefficient
    do m = 1, nmodes
       do mp = 1, nmodes
          if (mp == m) cycle
          new_c1(m,1) = new_c1(m,1) + c2_ij(m,mp) * q(mp)
       enddo
    enddo

    return

  end subroutine taylor_shift_offdiag

!######################################################################
! factorial: integer factorial. order1 is small (default 6, capped by
! input), so a plain loop suffices.
!######################################################################
  function factorial(n) result(f)

    implicit none

    integer, intent(in) :: n
    integer             :: f, i

    f = 1
    do i = 2, n
       f = f * i
    enddo

    return

  end function factorial

!######################################################################
! write_reexpand_log: emit a summary of the re-expansion to the log
! file. Reports Qmin in normal modes and Cartesians, the new diagonal
! energies e0(s), and the largest off-diagonal constant.
!######################################################################
  subroutine write_reexpand_log(qmin, vmin, ifault)

    use channels
    use sysinfo,    only: nmodes, nsta, natm, ncoo, ncoo, &
                          xcoo0, atlbl, nmcoo, nmlab
    use parameters, only: e0, e0_off
    use kdcglobal,  only: ireexpand_state
    use iomod,      only: freeunit

    implicit none

    real(dp), intent(in) :: qmin(nmodes)
    real(dp), intent(in) :: vmin
    integer,  intent(in) :: ifault

    integer  :: m, s, s1, s2, i, j, iout
    real(dp) :: x(ncoo), max_off
    integer  :: imax, jmax

!----------------------------------------------------------------------
! Section header
!----------------------------------------------------------------------
    write(ilog,'(/,72a)') ('+',i=1,72)
    write(ilog,'(2x,a,i0)') &
         'Re-expansion of the diabatic potential about the minimum '&
         //'of adiabatic state ', ireexpand_state
    write(ilog,'(72a)') ('+',i=1,72)

    if (ifault /= 0) then
       write(ilog,'(/,2x,a,i0)') &
            'WARNING: BFGS reported non-zero ifault = ', ifault
    endif

!----------------------------------------------------------------------
! Qmin in normal modes
!----------------------------------------------------------------------
    write(ilog,'(/,2x,a)') 'Qmin (normal modes):'
    write(ilog,'(2x,42a)') ('-',i=1,42)
    write(ilog,'(2x,a)')   '  Mode  Label        Q'
    write(ilog,'(2x,42a)') ('-',i=1,42)
    do m = 1, nmodes
       write(ilog,'(2x,i6,2x,a3,4x,F14.8)') m, nmlab(m), qmin(m)
    enddo
    write(ilog,'(2x,42a)') ('-',i=1,42)

!----------------------------------------------------------------------
! New diagonal energies
!----------------------------------------------------------------------
    write(ilog,'(/,2x,a)') &
         'New diagonal energies (eV) after re-expansion:'
    write(ilog,'(2x,28a)') ('-',i=1,28)
    write(ilog,'(2x,a)')   '  State        e0 (eV)'
    write(ilog,'(2x,28a)') ('-',i=1,28)
    do s = 1, nsta
       write(ilog,'(2x,i6,4x,F14.8)') s, e0(s)
    enddo
    write(ilog,'(2x,28a)') ('-',i=1,28)

!----------------------------------------------------------------------
! Largest off-diagonal zeroth-order constant
!----------------------------------------------------------------------
    max_off = 0.0d0
    imax    = 0
    jmax    = 0
    do s1 = 1, nsta-1
       do s2 = s1+1, nsta
          if (abs(e0_off(s1,s2)) > max_off) then
             max_off = abs(e0_off(s1,s2))
             imax    = s1
             jmax    = s2
          endif
       enddo
    enddo

    if (imax > 0) then
       write(ilog,'(/,2x,a,F14.8,a,2(x,i0))') &
            'Largest |e0_off| = ', max_off, &
            ' eV  for states', imax, jmax
    else
       write(ilog,'(/,2x,a)') &
            'All off-diagonal zeroth-order constants vanish.'
    endif

!----------------------------------------------------------------------
! Qmin in Cartesians, written to xyz file
!----------------------------------------------------------------------
    x = xcoo0/ang2bohr + matmul(nmcoo, qmin)

    call freeunit(iout)
    open(iout, file='reexpand_qmin.xyz', form='formatted', &
         status='unknown')
    write(iout,'(i0)') natm
    write(iout,'(a,i0)') &
         'Re-expansion centre Qmin: minimum of adiabatic state ', &
         ireexpand_state
    do i = 1, natm
       write(iout,'(a2,3(2x,F12.7))') atlbl(i), &
            (x(j), j=i*3-2, i*3)
    enddo
    close(iout)

    write(ilog,'(/,2x,a)') &
         'Qmin (Cartesians) written to reexpand_qmin.xyz'

    return

  end subroutine write_reexpand_log

!######################################################################
! reexpand_diagnostics: instrumented output to localise discrepancies
! between the diabatic and adiabatic optimiser paths, and to verify
! the analytical Taylor shift. Called twice from reexpand_about_qmin:
!   - lpost = .false. : before the shift; reports items 1-6 of the
!                       diagnostic plan (gradient norms, residual
!                       coupling, eigenvalues, max off-diagonal
!                       coefficients, Qmin_diab vs Qmin_adiab).
!   - lpost = .true.  : after the shift; reports item 7 (max e0_off
!                       in detail, broken down per state pair).
!######################################################################
  subroutine reexpand_diagnostics(qmin, s, lpost)

    use channels
    use sysinfo,    only: nmodes, nsta
    use parameters, only: e0_off, coeff1, coeff2, order1
    use potfuncs,   only: pot, dpot_dq, adiabaticgrad
    use geomopt,    only: optimise

    implicit none

    real(dp), intent(in) :: qmin(nmodes)
    integer,  intent(in) :: s
    logical,  intent(in) :: lpost

    integer  :: i, m, n, j, m1, m2, ifault_a
    real(dp) :: dw(nsta,nsta), wmat(nsta,nsta)
    real(dp) :: g_diab(nmodes), g_adiab(nmodes)
    real(dp) :: gn_diab, gn_adiab
    real(dp) :: max_off_pot, max_c1, max_c2, max_e0off
    integer  :: m_max1, j_max1, n_max1
    integer  :: m1_max2, m2_max2, j_max2
    integer  :: j_max_pot, j_max_e0off
    real(dp) :: ev(nsta), work(3*nsta)
    integer  :: info
    real(dp) :: q0(nmodes), qmin_a(nmodes), vmin_a, dq(nmodes)
    real(dp) :: qmin_diff_inf
    integer  :: m_qdiff
    real(dp), parameter :: tol_print = 1.0d-7

    if (lpost) then
       !----------------------------------------------------------------
       ! Item 7: max |e0_off(s,j)| after the Taylor shift
       !----------------------------------------------------------------
       write(ilog,'(/,72a)') ('+',i=1,72)
       write(ilog,'(2x,a)') 'Re-expansion diagnostics (post-shift)'
       write(ilog,'(72a)') ('+',i=1,72)

       max_e0off    = 0.0d0
       j_max_e0off  = 0
       do j = 1, nsta
          if (j == s) cycle
          if (abs(e0_off(s,j)) > max_e0off) then
             max_e0off   = abs(e0_off(s,j))
             j_max_e0off = j
          endif
       enddo
       write(ilog,'(/,2x,a,ES12.4,a,i0,a)') &
            'max |e0_off(s,j)|, j != s : ', max_e0off, &
            ' eV   (j = ', j_max_e0off, ')'

       return
    endif

    !-------------------------------------------------------------------
    ! Pre-shift diagnostics (items 1-6) at the converged qmin
    !-------------------------------------------------------------------
    write(ilog,'(/,72a)') ('+',i=1,72)
    write(ilog,'(2x,a)') 'Re-expansion diagnostics (pre-shift)'
    write(ilog,'(72a)') ('+',i=1,72)

    !-- Item 1: diabatic gradient norm |dW_{ss}/dQ| at qmin
    do m = 1, nmodes
       dw = dpot_dq(qmin, m)
       g_diab(m) = dw(s,s)
    enddo
    gn_diab = sqrt(dot_product(g_diab, g_diab))

    !-- Item 2: adiabatic gradient norm |dv_s/dQ| at qmin
    g_adiab = adiabaticgrad(qmin, s)
    gn_adiab = sqrt(dot_product(g_adiab, g_adiab))

    write(ilog,'(/,2x,a,ES12.4,a)') &
         '|grad W_{ss}(Qmin)|   = ', gn_diab,  '   (diabatic)'
    write(ilog,'(2x,a,ES12.4,a)')   &
         '|grad v_s(Qmin)|      = ', gn_adiab, '   (adiabatic)'

    !-- Item 3: max |W(s,j)(qmin)| for j != s
    wmat = pot(qmin)
    max_off_pot = 0.0d0
    j_max_pot   = 0
    do j = 1, nsta
       if (j == s) cycle
       if (abs(wmat(s,j)) > max_off_pot) then
          max_off_pot = abs(wmat(s,j))
          j_max_pot   = j
       endif
    enddo
    write(ilog,'(2x,a,ES12.4,a,i0,a)') &
         'max |W(s,j)(Qmin)|, j != s : ', max_off_pot, &
         ' eV   (j = ', j_max_pot, ')'

    !-- Item 4: eigenvalues of pot(qmin) and W(s,s)
    write(ilog,'(2x,a,F14.8,a)') &
         'W(s,s)(Qmin)          = ', wmat(s,s), ' eV'

    call dsyev('N','U',nsta,wmat,nsta,ev,work,3*nsta,info)
    if (info == 0) then
       write(ilog,'(2x,a)') 'Adiabatic eigenvalues v(:) at Qmin (eV):'
       do i = 1, nsta
          write(ilog,'(4x,a,i3,a,F14.8)') 'v(', i, ') = ', ev(i)
       enddo
    else
       write(ilog,'(2x,a,i0)') &
            'WARNING: dsyev failed at Qmin, info = ', info
    endif

    !-- Item 5: max |coeff1(m, s, j, n)| and max |coeff2(m1,m2, s, j)|
    !            for j != s, n = 1..order1
    max_c1  = 0.0d0
    m_max1  = 0
    j_max1  = 0
    n_max1  = 0
    do n = 1, order1
       do j = 1, nsta
          if (j == s) cycle
          do m = 1, nmodes
             if (abs(coeff1(m,s,j,n)) > max_c1) then
                max_c1 = abs(coeff1(m,s,j,n))
                m_max1 = m
                j_max1 = j
                n_max1 = n
             endif
          enddo
       enddo
    enddo

    if (max_c1 > 0.0d0) then
       write(ilog,'(2x,a,ES12.4,a,3(i0,a))') &
            'max |coeff1(m,s,j,n)|, j != s : ', max_c1, &
            ' eV   (m=', m_max1, ', j=', j_max1, ', n=', n_max1, ')'
    else
       write(ilog,'(2x,a)') &
            'max |coeff1(m,s,j,n)|, j != s : 0.0  (state s fully ' &
            //'decoupled in coeff1)'
    endif

    max_c2  = 0.0d0
    m1_max2 = 0
    m2_max2 = 0
    j_max2  = 0
    do j = 1, nsta
       if (j == s) cycle
       do m2 = 1, nmodes
          do m1 = 1, nmodes
             if (abs(coeff2(m1,m2,s,j)) > max_c2) then
                max_c2  = abs(coeff2(m1,m2,s,j))
                m1_max2 = m1
                m2_max2 = m2
                j_max2  = j
             endif
          enddo
       enddo
    enddo

    if (max_c2 > 0.0d0) then
       write(ilog,'(2x,a,ES12.4,a,3(i0,a))') &
            'max |coeff2(m1,m2,s,j)|, j != s : ', max_c2, &
            ' eV   (m1=', m1_max2, ', m2=', m2_max2, &
            ', j=', j_max2, ')'
    else
       write(ilog,'(2x,a)') &
            'max |coeff2(m1,m2,s,j)|, j != s : 0.0  (state s fully ' &
            //'decoupled in coeff2)'
    endif

    !-- Item 6: re-run optimise() with the diabatic path and compare
    !  to the adiabatic Qmin found by reexpand_about_qmin. Under the
    !  decoupled-state-s precondition, the two minima must coincide.
    q0 = 0.0d0
    call optimise(q0, s, qmin_a, vmin=vmin_a, ifault=ifault_a, &
                  tol_g=1.0d-8, ldiabatic=.true.)

    if (ifault_a /= 0) then
       write(ilog,'(2x,a,i0,a)') &
            'WARNING: diabatic BFGS did not converge (ifault = ', &
            ifault_a, ')'
    endif

    dq = qmin - qmin_a
    qmin_diff_inf = maxval(abs(dq))
    m_qdiff = maxloc(abs(dq), dim=1)

    write(ilog,'(2x,a,ES12.4,a,i0,a)') &
         '||Qmin_adiab - Qmin_diab||_inf = ', qmin_diff_inf, &
         '   (largest in mode ', m_qdiff, ')'

    if (qmin_diff_inf > tol_print) then
       write(ilog,'(2x,a)') &
            'Per-mode (Qmin_diab - Qmin_adiab):'
       do m = 1, nmodes
          if (abs(dq(m)) > tol_print) then
             write(ilog,'(4x,a,i3,a,ES14.6)') &
                  'm = ', m, ' :  dQ = ', dq(m)
          endif
       enddo
    endif

    return

  end subroutine reexpand_diagnostics

!######################################################################

end module reexpand

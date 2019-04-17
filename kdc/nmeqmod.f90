module nmeqmod

  implicit none

contains

!######################################################################

  subroutine get_coefficients_nmeq

    use constants
    use channels
    use iomod
    use sysinfo
    use kdcglobal
    
    implicit none

!----------------------------------------------------------------------
! Perform the fitting for the 1-mode terms
!----------------------------------------------------------------------
    call fit_1mode_terms

!----------------------------------------------------------------------
! Convert all coupling coefficients to units of eV
!----------------------------------------------------------------------
    kappa=kappa*eh2ev
    lambda=lambda*eh2ev
    gamma=gamma*eh2ev
    mu=mu*eh2ev
    iota=iota*eh2ev
    tau=tau*eh2ev
    epsilon=epsilon*eh2ev
    xi=xi*eh2ev
    
    return
    
  end subroutine get_coefficients_nmeq

!######################################################################

  subroutine fit_1mode_terms

    use constants
    use channels
    use iomod
    use sysinfo
    use kdcglobal
    
    implicit none

    integer               :: m,n,s1,s2,ndat
    integer               :: order
    real(dp), allocatable :: coeff(:),q(:),w(:)
    logical               :: lpseudo
    
!----------------------------------------------------------------------
! Determine the sets of files corresponding to single mode
! displacements
!----------------------------------------------------------------------
    call get_indices_1mode

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    ! Order of the polynomial to be fit (hard-wired for now)
    order=4

    ! Coefficient vector
    allocate(coeff(order+1))
    coeff=0.0d0

    ! Input coordinates and diabatic potential values
    allocate(q(maxfiles))
    allocate(w(maxfiles))
    q=0.0d0
    w=0.0d0
    
!----------------------------------------------------------------------
! Perform the fits of the 1-mode terms
!----------------------------------------------------------------------
    ! Loop over modes
    do m=1,nmodes

       ! Cycle if there are no points for the current mode
       ndat=nfiles1m(m)
       if (ndat.eq.0) cycle

       ! Loop over elements of the diabatic potential matrix
       do s1=1,nsta
          do s2=s1,nsta

             ! Fill in the coordinate and potential vectors to be
             ! sent to the fitting routine
             q=0.0d0
             w=0.0d0
             do n=1,ndat
                q(n)=qvec(m,findx1m(m,n))
                w(n)=diabpot(s1,s2,findx1m(m,n))
             enddo
             
             ! Perform the fitting for the current mode and
             ! diabatic potential matrix element
             call nmeq1d(order,coeff,ndat,q(1:ndat),w(1:ndat),&
                  lpseudo)
             
             ! Output a warning if the psedo-inverse was used in
             ! the fitting
             if (lpseudo) write(ilog,'(/,a,i0,a,i0,x,i0)') &
                  'WARNING: Pseudoinverse used in the fitting of the &
                  coefficients for mode ',m,' and states ',s1,s2
             
             ! Fill in the global coefficient arrays
             call fill_coeffs1d(coeff,order,m,s1,s2)
             
          enddo
       enddo
             
    enddo

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(coeff)
    deallocate(q)
    deallocate(w)
    
    return
    
  end subroutine fit_1mode_terms
    
!######################################################################

  subroutine get_indices_1mode

    use constants
    use channels
    use iomod
    use sysinfo
    use kdcglobal
    
    implicit none

    integer             :: m,n,mindx,ndisp
    real(dp), parameter :: thrsh=1e-4_dp
    
!----------------------------------------------------------------------
! Determine the no. files/geometries corresponding to single-mode
! displacements for each mode
!----------------------------------------------------------------------
    ! Allocate arrays
    allocate(nfiles1m(nmodes))
    nfiles1m=0.0d0

    ! Loop over displaced geometries
    do n=1,nfiles

       ! Determine the no. displaced modes
       ndisp=0
       do m=1,nmodes
          if (abs(qvec(m,n)).gt.thrsh) then
             ndisp=ndisp+1
             mindx=m
          endif
       enddo

       ! If only one mode is displaced, then update the
       ! corresponding element of nfiles1m
       if (ndisp.eq.1) nfiles1m(mindx)=nfiles1m(mindx)+1
          
    enddo

!----------------------------------------------------------------------
! Determine the sets of files corresponding to single mode
! displacements
!----------------------------------------------------------------------
    ! Allocate arrays
    maxfiles=maxval(nfiles1m)
    allocate(findx1m(nmodes,maxfiles))
    findx1m=0

    ! Reset the nfiles1m array to act as a counter (it will also be
    ! re-filled in in the following)
    nfiles1m=0
    
    ! Loop over displaced geometries
    do n=1,nfiles

       ! Determine the no. displaced modes
       ndisp=0
       do m=1,nmodes
          if (abs(qvec(m,n)).gt.thrsh) then
             ndisp=ndisp+1
             mindx=m
          endif
       enddo

       ! Fill in the findx1m array
       if (ndisp.eq.1) then
          nfiles1m(mindx)=nfiles1m(mindx)+1
          findx1m(mindx,nfiles1m(mindx))=n
       endif
       
    enddo

    return
    
  end subroutine get_indices_1mode

!######################################################################

  subroutine nmeq1d(order,coeff,npnts,x,y,lpseudo)

    use constants
    use channels
    use iomod
    use utils
    
    implicit none

    integer, intent(in)                    :: order,npnts
    integer                                :: i,j
    real(dp), dimension(order+1)           :: coeff
    real(dp), intent(in), dimension(npnts) :: x,y
    real(dp), dimension(npnts,order+1)     :: Xmat
    real(dp), dimension(order+1,order+1)   :: XTX,invXTX
    logical                                :: lpseudo
    
!----------------------------------------------------------------------
! Fill in the X-matrix
!----------------------------------------------------------------------
    do i=1,npnts
       do j=1,order+1
          Xmat(i,j)=x(i)**(j-1)
       enddo
    enddo

!----------------------------------------------------------------------
! Calculate the coefficients
!----------------------------------------------------------------------
    XTX=matmul(transpose(Xmat),Xmat)

    call invert_matrix(XTX,invXTX,order+1,lpseudo)

    coeff=matmul(invXTX,matmul(transpose(Xmat),y))
    
    return
    
  end subroutine nmeq1d

!######################################################################

  subroutine fill_coeffs1d(coeff,order,m,s1,s2)

    use constants
    use channels
    use iomod
    use sysinfo
    use kdcglobal
    
    implicit none

    integer, intent(in)                      :: order,m,s1,s2
    real(dp), intent(in), dimension(order+1) :: coeff

!----------------------------------------------------------------------
! Intrastate coupling coefficients
!----------------------------------------------------------------------
    if (s1.eq.s2) then
       if (order.gt.0) kappa(m,s1)=coeff(2)
       if (order.gt.1) gamma(m,m,s1)=2.0d0*coeff(3)-freq(m)/eh2ev
       if (order.gt.2) iota(m,s1)=6.0d0*coeff(4)
       if (order.gt.3) epsilon(m,s1)=24.0d0*coeff(5)
    endif

!----------------------------------------------------------------------
! Intrastate coupling coefficients
!----------------------------------------------------------------------
    if (s1.ne.s2) then
       if (order.gt.0) then
          lambda(m,s1,s2)=coeff(2)
          lambda(m,s2,s1)=lambda(m,s1,s2)
       endif
       if (order.gt.1) then
          mu(m,m,s1,s2)=2.0d0*coeff(3)
          mu(m,m,s2,s1)=mu(m,m,s1,s2)
       endif
       if (order.gt.2) then
          tau(m,s1,s2)=6.0d0*coeff(4)
          tau(m,s2,s1)=tau(m,s1,s2)
       endif
       if (order.gt.3) then
          xi(m,s1,s2)=24.0d0*coeff(5)
          xi(m,s2,s1)=xi(m,s1,s2)
       endif
    endif
    
    return
    
  end subroutine fill_coeffs1d
    
!######################################################################
  
end module nmeqmod

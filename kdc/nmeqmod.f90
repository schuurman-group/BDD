module nmeqmod

  implicit none

contains

!######################################################################

  subroutine get_coefficients_nmeq

    use constants
    use channels
    use iomod
    use sysinfo
    use parameters
    use kdcglobal
    
    implicit none

!----------------------------------------------------------------------
! Perform the fitting for the 1-mode terms
!----------------------------------------------------------------------
    call fit_1mode_terms
    
!----------------------------------------------------------------------
! Perform the fitting of the 2-mode terms
! Note that this has to be done after the 1-mode terms have been
! calculated
!----------------------------------------------------------------------
    call fit_2mode_terms
    
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

!----------------------------------------------------------------------
! Calculate the RMSD of the fit
!----------------------------------------------------------------------
    call calc_rmsd

    return
    
  end subroutine get_coefficients_nmeq

!######################################################################

  subroutine fit_1mode_terms

    use kdcglobal
    
    implicit none
    
!----------------------------------------------------------------------
! Determine the sets of files corresponding to single mode
! displacements
!----------------------------------------------------------------------
    call get_indices_1mode

!----------------------------------------------------------------------
! Fit the 1-mode terms entering into the vibronic coupling Hamiltonian
!----------------------------------------------------------------------
    call fit_1mode_terms_potential

!----------------------------------------------------------------------
! Fit the 1-mode terms entering into the expansion of the dipole
! matrix
!----------------------------------------------------------------------
    if (ldipfit) call fit_1mode_terms_dipole
    
    return
    
  end subroutine fit_1mode_terms

!######################################################################

  subroutine fit_1mode_terms_potential

    use constants
    use channels
    use iomod
    use sysinfo
    use parameters
    use kdcglobal
    
    implicit none

    integer               :: m,n,s1,s2,ndat
    integer               :: order
    real(dp), allocatable :: coeff(:),q(:),w(:)
    real(dp)              :: wfac1
    logical               :: lpseudo

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    ! Order of the polynomial to be fit (hard-wired for now)
    order=4

    ! Coefficient vector
    allocate(coeff(order))
    coeff=0.0d0

    ! Input coordinates and diabatic potential values
    allocate(q(maxfiles1m))
    allocate(w(maxfiles1m))
    q=0.0d0
    w=0.0d0
    
!----------------------------------------------------------------------
! Perform the fits of the 1-mode terms
!----------------------------------------------------------------------
    ! Loop over modes
    do m=1,nmodes

       ! Cycle if there are no points for the current mode
       ndat=ngeom1m(m)
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

             ! Subtract off the zeroth-order potential value
             ! Note that WIJ(Q0)=0 for I!=J
             if (s1.eq.s2) w(1:ndat)=w(1:ndat)-q0pot(s1)

             ! Perform the fitting for the current mode and
             ! diabatic potential matrix element
             if (s1.eq.s2) then
                wfac1=wfac
             else
                wfac1=0.0d0
             endif
             call nmeq1d(order,coeff,ndat,q(1:ndat),w(1:ndat),&
                  lpseudo,wfac1)
             
             ! Output a warning if the psedo-inverse was used in
             ! the fitting
             if (lpseudo) write(ilog,'(/,a,i0,a,i0,x,i0)') &
                  'WARNING: Pseudoinverse used in the fitting of the &
                  coefficients for mode ',m,' and states ',s1,s2
             
             ! Fill in the global coefficient arrays
             call fill_coeffs1d_potential(coeff,order,m,s1,s2)
             
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
    
  end subroutine fit_1mode_terms_potential

!######################################################################

  subroutine fit_1mode_terms_dipole

    use constants
    use channels
    use iomod
    use sysinfo
    use kdcglobal
    
    implicit none

    integer               :: m,n,s1,s2,c,ndat
    integer               :: order
    real(dp), allocatable :: coeff(:),q(:),d(:)
    logical               :: lpseudo

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    ! Order of the polynomial to be fit (hard-wired for now)
    order=4

    ! Coefficient vector
    allocate(coeff(order))
    coeff=0.0d0

    ! Input coordinates and diabatic dipole matrix element values
    allocate(q(maxfiles1m))
    allocate(d(maxfiles1m))
    q=0.0d0
    d=0.0d0

!----------------------------------------------------------------------
! Perform the fits of the 1-mode terms
!----------------------------------------------------------------------
    ! Loop over modes
    do m=1,nmodes
       
       ! Cycle if there are no points for the current mode
       ndat=ngeom1m(m)
       if (ndat.eq.0) cycle
       
       ! Loop over elements of the diabatic dipole matrix
       do s1=1,nsta
          do s2=s1,nsta
             
             ! Loop over components of the dipole
             do c=1,3
                
                ! Fill in the coordinate and dipole matrix element
                ! vectors to be sent to the fitting routine
                q=0.0d0
                d=0.0d0
                do n=1,ndat
                   q(n)=qvec(m,findx1m(m,n))
                   d(n)=diabdip(s1,s2,c,findx1m(m,n))
                enddo

                ! Perform the fitting for the current mode and
                ! diabatic dipole matrix element
                call nmeq1d(order,coeff,ndat,q(1:ndat),d(1:ndat),&
                     lpseudo,0.0d0)

                ! Output a warning if the psedo-inverse was used in
                ! the fitting
                if (lpseudo) write(ilog,'(/,a,i0,a,i0,x,i0,a,i0)') &
                     'WARNING: Pseudoinverse used in the fitting of &
                     the coefficients for mode ',m,', states ',s1,s2,&
                     ' and dipole component ',c

                ! Fill in the global coefficient arrays
                call fill_coeffs1d_dipole(coeff,order,m,s1,s2,c)
                
             enddo
          enddo
                
       enddo
          
    enddo
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(coeff)
    deallocate(q)
    deallocate(d)
    
    return
    
  end subroutine fit_1mode_terms_dipole
    
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
    allocate(ngeom1m(nmodes))
    ngeom1m=0

    ! Loop over displaced geometries
    do n=1,ngeom

       ! Determine the no. displaced modes
       ndisp=0
       do m=1,nmodes
          if (abs(qvec(m,n)).gt.thrsh) then
             ndisp=ndisp+1
             mindx=m
          endif
       enddo

       ! If only one mode is displaced, then update the
       ! corresponding element of ngeom1m
       if (ndisp.eq.1) ngeom1m(mindx)=ngeom1m(mindx)+1
          
    enddo

!----------------------------------------------------------------------
! Determine the sets of files corresponding to single mode
! displacements
!----------------------------------------------------------------------
    ! Allocate arrays
    maxfiles1m=maxval(ngeom1m)
    allocate(findx1m(nmodes,maxfiles1m))
    findx1m=0

    ! Reset the ngeom1m array to act as a counter (it will also be
    ! re-filled in in the following)
    ngeom1m=0
    
    ! Loop over displaced geometries
    do n=1,ngeom

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
          ngeom1m(mindx)=ngeom1m(mindx)+1
          findx1m(mindx,ngeom1m(mindx))=n
       endif
       
    enddo

    return
    
  end subroutine get_indices_1mode

!######################################################################

  subroutine nmeq1d(order,coeff,npnts,x,y,lpseudo,wfac)

    use constants
    use channels
    use iomod
    use utils
    
    implicit none

    integer, intent(in)              :: order,npnts
    integer                          :: i,j
    real(dp), dimension(order)       :: coeff
    real(dp), dimension(npnts)       :: x,y
    real(dp), dimension(npnts,order) :: Xmat
    real(dp), dimension(order,order) :: XTX,invXTX
    real(dp), intent(in)             :: wfac
    logical                          :: lpseudo

    real(dp)                   :: miny
    real(dp), dimension(npnts) :: wsq

!----------------------------------------------------------------------
! Fill in the X-matrix
!----------------------------------------------------------------------
    do i=1,npnts
       do j=1,order
          Xmat(i,j)=x(i)**j
       enddo
    enddo

!----------------------------------------------------------------------
! Weights
!----------------------------------------------------------------------
    miny=minval(y)
    do i=1,npnts
       wsq(i)=exp(-wfac*(y(i)-miny))**2
    enddo

    do i=1,npnts
       y(i)=y(i)*wsq(i)
       Xmat(i,:)=Xmat(i,:)*wsq(i)
    enddo
    
!----------------------------------------------------------------------
! Calculate the coefficients
!----------------------------------------------------------------------
    XTX=matmul(transpose(Xmat),Xmat)

    call invert_matrix(XTX,invXTX,order,lpseudo)

    coeff=matmul(invXTX,matmul(transpose(Xmat),y))
    
    return
    
  end subroutine nmeq1d

!######################################################################

  subroutine fill_coeffs1d_potential(coeff,order,m,s1,s2)

    use constants
    use channels
    use iomod
    use sysinfo
    use parameters
    use kdcglobal
    
    implicit none

    integer, intent(in)                    :: order,m,s1,s2
    real(dp), intent(in), dimension(order) :: coeff

!----------------------------------------------------------------------
! Intrastate coupling coefficients
!----------------------------------------------------------------------
    if (s1.eq.s2) then
       if (order.gt.0) kappa(m,s1)=coeff(1)
       if (order.gt.1) gamma(m,m,s1)=2.0d0*coeff(2)-freq(m)/eh2ev
       if (order.gt.2) iota(m,s1)=6.0d0*coeff(3)
       if (order.gt.3) epsilon(m,s1)=24.0d0*coeff(4)
    endif

!----------------------------------------------------------------------
! Intrastate coupling coefficients
!----------------------------------------------------------------------
    if (s1.ne.s2) then
       if (order.gt.0) then
          lambda(m,s1,s2)=coeff(1)
          lambda(m,s2,s1)=lambda(m,s1,s2)
       endif
       if (order.gt.1) then
          mu(m,m,s1,s2)=2.0d0*coeff(2)
          mu(m,m,s2,s1)=mu(m,m,s1,s2)
       endif
       if (order.gt.2) then
          tau(m,s1,s2)=6.0d0*coeff(3)
          tau(m,s2,s1)=tau(m,s1,s2)
       endif
       if (order.gt.3) then
          xi(m,s1,s2)=24.0d0*coeff(4)
          xi(m,s2,s1)=xi(m,s1,s2)
       endif
    endif
    
    return
    
  end subroutine fill_coeffs1d_potential

!######################################################################
  
  subroutine fill_coeffs1d_dipole(coeff,order,m,s1,s2,c)

    use constants
    use channels
    use iomod
    use sysinfo
    use parameters
    use kdcglobal
    
    implicit none

    integer, intent(in)                      :: order,m,s1,s2,c
    real(dp), intent(in), dimension(order+1) :: coeff

    ! 1st-order terms
    if (order.gt.0) then
       dip1(m,s1,s2,c)=coeff(2)
       dip1(m,s2,s1,c)=dip1(m,s1,s2,c)
    endif

    ! 2nd-order terms
    if (order.gt.1) then
       dip2(m,m,s1,s2,c)=2.0d0*coeff(3)
       dip2(m,m,s2,s1,c)=dip2(m,m,s1,s2,c)
    endif

    ! 3rd-order terms
    if (order.gt.2) then
       dip3(m,s1,s2,c)=6.0d0*coeff(4)
       dip3(m,s2,s1,c)=dip3(m,s1,s2,c)
    endif

    ! 4th-order terms
    if (order.gt.3) then
       dip4(m,s1,s2,c)=24.0d0*coeff(5)
       dip4(m,s2,s1,c)=dip4(m,s1,s2,c)
    endif
    
    return
    
  end subroutine fill_coeffs1d_dipole
    
!######################################################################

  subroutine fit_2mode_terms

    use constants
    use channels
    use iomod
    use sysinfo
    use kdcglobal
    
    implicit none

    integer               :: n,m1,m2,s1,s2,ndat
    real(dp), allocatable :: q(:,:),w(:)
    real(dp)              :: coeff
    logical               :: present

!----------------------------------------------------------------------
! Determine the sets of files corresponding to single mode
! displacements
!----------------------------------------------------------------------
    call get_indices_2modes(present)

!----------------------------------------------------------------------
! Return here if no files containing two-mode displacements were found
!----------------------------------------------------------------------
    if (.not.present) return

!----------------------------------------------------------------------
! Fit the 2-mode terms entering into the vibronic coupling Hamiltonian
!----------------------------------------------------------------------
    !call fit_2mode_terms_potential

    call fit_2mode_terms_potential_new
    
!----------------------------------------------------------------------
! Fit the 2-mode terms entering into the expansion of the dipole
! matrix
!----------------------------------------------------------------------
    if (ldipfit) then

       print*,'2-mode dipole matrix fitting needs updating'
       stop
       
       call fit_2mode_terms_dipole

    endif
    
    return
    
  end subroutine fit_2mode_terms

!######################################################################

  subroutine fit_2mode_terms_potential

    use constants
    use channels
    use iomod
    use sysinfo
    use kdcglobal
    
    implicit none

    integer               :: n,m1,m2,s1,s2,ndat,mask
    real(dp), allocatable :: q(:,:),w(:)
    real(dp)              :: coeff
    logical               :: present

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    ! Input coordinates and diabatic potential values
    allocate(q(2,maxfiles2m))
    allocate(w(maxfiles2m))
    q=0.0d0
    w=0.0d0

!----------------------------------------------------------------------
! Perform the fits of the 2-mode terms
!----------------------------------------------------------------------
    ! Loop over pairs of modes
    do m1=1,nmodes-1
       do m2=m1+1,nmodes

          ! Cycle if there are no points for the current pair of modes
          ndat=ngeom2m(m1,m2)
          if (ndat.eq.0) cycle

          ! Loop over elements of the diabatic potential matrix
          do s1=1,nsta
             do s2=s1,nsta

                ! Fill in the coordinate and potential vectors to be
                ! sent to the fitting routine
                q=0.0d0
                w=0.0d0
                do n=1,ndat
                   q(1,n)=qvec(m1,findx2m(m1,m2,n))
                   q(2,n)=qvec(m2,findx2m(m1,m2,n))
                   w(n)=diabpot(s1,s2,findx2m(m1,m2,n))
                enddo

                ! Subtract off the one-mode contributions to the
                ! potential to get the correlated contribution
                call get_2mode_contrib_potential(w(1:ndat),q(:,1:ndat),&
                     ndat,m1,m2,s1,s2)

                ! Skip the fitting if the current bi-linear parameter is
                ! zero by symmetry
                mask=mask_pot_2mode(m1,m2,s1,s2)
                if (mask.eq.0) cycle
                
                ! Perform the fitting for the current pair of modes and
                ! diabatic potential matrix element
                call nmeq_bilinear(coeff,ndat,q(:,1:ndat),w(1:ndat))

                ! Fill in the global coefficient arrays
                call fill_coeffs2d_potential(coeff,m1,m2,s1,s2)
                
             enddo
          enddo
                
       enddo
    enddo

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(q)
    deallocate(w)
    
    return
    
  end subroutine fit_2mode_terms_potential

!######################################################################

  subroutine fit_2mode_terms_potential_new

    use constants
    use channels
    use iomod
    use sysinfo
    use kdcglobal

    use parameters
    
    implicit none

    integer                  :: n,m1,m2,s1,s2,ndat,mask
    integer                  :: order
    real(dp), allocatable    :: coeff(:),q(:),w(:)
    real(dp)                 :: wfac1,val
    real(dp), dimension(2)   :: qi,qf
    real(dp), dimension(2,2) :: U
    logical                  :: present,lpseudo

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    ! Order of the polynomial to be fit (hard-wired for now)
    order=4

    ! Coefficient vector
    allocate(coeff(order))
    coeff=0.0d0

    ! Input coordinates and diabatic potential values
    allocate(q(maxfiles2m))
    allocate(w(maxfiles2m))
    q=0.0d0
    w=0.0d0

!----------------------------------------------------------------------
! Transformation matrix
!----------------------------------------------------------------------
    U(1,1)=1.0d0/sqrt(2.0d0)
    U(2,1)=1.0d0/sqrt(2.0d0)
    U(1,2)=-1.0d0/sqrt(2.0d0)
    U(2,2)=1.0d0/sqrt(2.0d0)
    
!----------------------------------------------------------------------
! Perform the fits of the 2-mode terms
!----------------------------------------------------------------------
    ! Loop over pairs of modes
    do m1=1,nmodes-1
       do m2=m1+1,nmodes

          ! Cycle if there are no points for the current pair of modes
          ndat=ngeom2m(m1,m2)
          if (ndat.eq.0) cycle

          ! Loop over elements of the diabatic potential matrix
          do s1=1,nsta
             do s2=s1,nsta

                ! Fill in the coordinate and potential vectors to be
                ! sent to the fitting routine
                q=0.0d0
                w=0.0d0
                do n=1,ndat
                   ! Normal mode coordinates in the original basis
                   qi(1)=qvec(m1,findx2m(m1,m2,n))
                   qi(2)=qvec(m2,findx2m(m1,m2,n))
                   ! Normal mode coordinates in the rotated basis
                   qf=matmul(U,qi)
                   ! Fill in the coordinate and potential vectors
                   q(n)=qf(2)
                   w(n)=diabpot(s1,s2,findx2m(m1,m2,n))
                enddo

                ! Subtract off the zeroth-order potential value
                ! Note that WIJ(Q0)=0 for I!=J
                if (s1.eq.s2) w(1:ndat)=w(1:ndat)-q0pot(s1)

                ! Perform the fitting for the current mode and
                ! diabatic potential matrix element
                if (s1.eq.s2) then
                   wfac1=wfac
                else
                   wfac1=0.0d0
                endif

                call nmeq1d(order,coeff,ndat,q(1:ndat),w(1:ndat),&
                     lpseudo,wfac1)

                ! Check on the 1st-order coupling coefficient
                if (s1 == s2) then
                   val=(kappa(m1,s1)+kappa(m2,s2))/sqrt(2.0d0)
                else
                   val=(lambda(m1,s1,s2)+lambda(m2,s1,s2))/sqrt(2.0d0)
                endif
                if (abs(val-coeff(1)) > 0.01d0) then
                   write(6,'(/,x,a)') &
                        'Warning: inconsistency in the 2-mode fits'
                endif

                ! Fill in the global coefficient arrays
                call fill_coeffs2d_potential_new(coeff,order,m1,m2,s1,s2)
                
             enddo
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
    
  end subroutine fit_2mode_terms_potential_new
    
!######################################################################

  subroutine fit_2mode_terms_dipole

    use constants
    use channels
    use iomod
    use sysinfo
    use kdcglobal
    
    implicit none

    integer               :: n,m1,m2,s1,s2,c,ndat
    real(dp), allocatable :: q(:,:),d(:)
    real(dp)              :: coeff
    logical               :: present

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    ! Input coordinates and diabatic dipole matrix element values
    allocate(q(2,maxfiles2m))
    allocate(d(maxfiles2m))
    q=0.0d0
    d=0.0d0

!----------------------------------------------------------------------
! Perform the fits of the 2-mode terms
!----------------------------------------------------------------------
    ! Loop over pairs of modes
    do m1=1,nmodes-1
       do m2=m1+1,nmodes

          ! Cycle if there are no points for the current pair of modes
          ndat=ngeom2m(m1,m2)
          if (ndat.eq.0) cycle

          ! Loop over elements of the diabatic dipole matrix
          do s1=1,nsta
             do s2=s1,nsta

                ! Loop over components of the dipole
                do c=1,3
                
                   ! Fill in the coordinate and dipole matrix element
                   ! vectors to be sent to the fitting routine
                   q=0.0d0
                   d=0.0d0
                   do n=1,ndat
                      q(1,n)=qvec(m1,findx2m(m1,m2,n))
                      q(2,n)=qvec(m2,findx2m(m1,m2,n))
                      d(n)=diabdip(s1,s2,c,findx2m(m1,m2,n))
                   enddo
                
                   ! Subtract off the one-mode contributions to the
                   ! potential to get the correlated contribution
                   call get_2mode_contrib_dipole(d(1:ndat),q(:,1:ndat),&
                        ndat,m1,m2,s1,s2,c)
                
                   ! Perform the fitting for the current pair of modes
                   ! and diabatic dipole matrix element
                   call nmeq_bilinear(coeff,ndat,q(:,1:ndat),d(1:ndat))
                   
                   ! Fill in the global coefficient arrays
                   call fill_coeffs2d_dipole(coeff,m1,m2,s1,s2,c)

                enddo
                   
             enddo
          enddo

       enddo
    enddo
          
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(q)
    deallocate(d)
    
    return
    
  end subroutine fit_2mode_terms_dipole
    
!######################################################################

  subroutine get_indices_2modes(present)

    use constants
    use channels
    use iomod
    use sysinfo
    use kdcglobal
    
    implicit none

    integer                     :: m,n,mindx1,mindx2,ndisp
    real(dp), parameter         :: thrsh=1e-4_dp
    real(dp), dimension(nmodes) :: qunit
    real(dp)                    :: theta,inner,norm,diff
    logical                     :: present

!----------------------------------------------------------------------
! Determine the no. files/geometries corresponding to two-mode
! displacements for each mode
!----------------------------------------------------------------------
    ! Allocate arrays
    allocate(ngeom2m(nmodes,nmodes))
    ngeom2m=0

    ! Loop over displaced geometries
    do n=1,ngeom

       ! Determine the no. displaced modes
       ndisp=0
       do m=1,nmodes
          if (abs(qvec(m,n)).gt.thrsh) then
             if (ndisp.eq.0) then
                mindx1=m                
             else
                mindx2=m
             endif
             ndisp=ndisp+1
          endif
       enddo

       ! If only two modes are displaced diagonally, then update the
       ! corresponding element of ngeom2m
       if (ndisp.eq.2) then

          ! Unit vector pointing along the diagonal of the two modes
          qunit=0.0d0
          qunit(mindx1)=1.0d0/sqrt(2.0d0)
          qunit(mindx2)=qunit(mindx1)

          ! Angle between the above unit vector and the current
          ! normal mode coordinate vector
          inner=dot_product(qvec(:,n),qunit)
          norm=sqrt(dot_product(qvec(:,n),qvec(:,n)))
          theta=acos(inner/norm)

          ! Are we at a diagonal point?
          diff=min(abs(theta),abs(theta-pi))
          if (diff < thrsh) then
             ngeom2m(mindx1,mindx2)=ngeom2m(mindx1,mindx2)+1
             ngeom2m(mindx2,mindx1)=ngeom2m(mindx1,mindx2)
          endif

          !ngeom2m(mindx1,mindx2)=ngeom2m(mindx1,mindx2)+1
          !ngeom2m(mindx2,mindx1)=ngeom2m(mindx1,mindx2)
          
       endif
          
    enddo

!----------------------------------------------------------------------
! Return if no files containing two-mode displacements were found
!----------------------------------------------------------------------
    if (maxval(ngeom2m).gt.0) then
       present=.true.
    else
       present=.false.
       return
    endif
       
!----------------------------------------------------------------------
! Determine the sets of files corresponding to two-mode displacements
!----------------------------------------------------------------------
    ! Allocate arrays
    maxfiles2m=maxval(ngeom2m)
    allocate(findx2m(nmodes,nmodes,maxfiles2m))
    findx2m=0

    ! Reset the ngeom1m array to act as a counter (it will also be
    ! re-filled in in the following)
    ngeom2m=0

    ! Loop over displaced geometries
    do n=1,ngeom

       ! Determine the no. displaced modes
       ndisp=0
       do m=1,nmodes
          if (abs(qvec(m,n)).gt.thrsh) then
             if (ndisp.eq.0) then
                mindx1=m                
             else
                mindx2=m
             endif
             ndisp=ndisp+1
          endif
       enddo
       
       ! Fill in the findx2m array
       if (ndisp.eq.2) then

          ! Unit vector pointing along the diagonal of the two modes
          qunit=0.0d0
          qunit(mindx1)=1.0d0/sqrt(2.0d0)
          qunit(mindx2)=qunit(mindx1)

          ! Angle between the above unit vector and the current
          ! normal mode coordinate vector
          inner=dot_product(qvec(:,n),qunit)
          norm=sqrt(dot_product(qvec(:,n),qvec(:,n)))
          theta=acos(inner/norm)

          ! Are we at a diagonal point?
          diff=min(abs(theta),abs(theta-pi))
          if (diff < thrsh) then
             ngeom2m(mindx1,mindx2)=ngeom2m(mindx1,mindx2)+1
             ngeom2m(mindx2,mindx1)=ngeom2m(mindx1,mindx2)
             findx2m(mindx1,mindx2,ngeom2m(mindx1,mindx2))=n
             findx2m(mindx2,mindx1,ngeom2m(mindx1,mindx2))=n
          endif
          
          !ngeom2m(mindx1,mindx2)=ngeom2m(mindx1,mindx2)+1
          !ngeom2m(mindx2,mindx1)=ngeom2m(mindx1,mindx2)
          !findx2m(mindx1,mindx2,ngeom2m(mindx1,mindx2))=n
          !findx2m(mindx2,mindx1,ngeom2m(mindx1,mindx2))=n
       endif
       
    enddo
       
    return
    
  end subroutine get_indices_2modes

!######################################################################

  subroutine get_2mode_contrib_potential(w,q,ndat,m1,m2,s1,s2)

    use constants
    use channels
    use iomod
    use sysinfo
    use parameters
    use kdcglobal
    
    implicit none

    integer, intent(in)                     :: ndat,m1,m2,s1,s2
    integer                                 :: n
    real(dp), dimension(ndat)               :: w
    real(dp), dimension(2,ndat), intent(in) :: q

!----------------------------------------------------------------------
! On-diagonal diabatic potential matrix element
!----------------------------------------------------------------------
    if (s1.eq.s2) then

       ! Loop over points
       do n=1,ndat

          ! Zeroth-order contributions
          w(n)=w(n)-q0pot(s1)
          w(n)=w(n)-0.5d0*(freq(m1)*q(1,n)**2+freq(m2)*q(2,n)**2)/eh2ev
          
          ! First-order contributions
          w(n)=w(n)-kappa(m1,s1)*q(1,n)-kappa(m2,s1)*q(2,n)

          ! Second-order contributions
          w(n)=w(n)-0.5d0*(gamma(m1,m1,s1)*q(1,n)**2 &
               +gamma(m2,m2,s1)*q(2,n)**2)

          ! Third-order contributions
          w(n)=w(n)-(1.0d0/6.0d0)*(iota(m1,s1)*q(1,n)**3 &
              +iota(m2,s1)*q(2,n)**3)

          ! Fourth-order contributions
          w(n)=w(n)-(1.0d0/24.0d0)*(epsilon(m1,s1)*q(1,n)**4 &
              +epsilon(m2,s1)*q(2,n)**4)
          
       enddo
       
    endif

!----------------------------------------------------------------------
! Off-diagonal diabatic potential matrix element
!----------------------------------------------------------------------
    if (s1.ne.s2) then

       ! Loop over points
       do n=1,ndat
          
          ! First-order contributions
          w(n)=w(n)-lambda(m1,s1,s2)*q(1,n)-lambda(m2,s1,s2)*q(2,n)

          ! Second-order contributions
          w(n)=w(n)-0.5d0*(mu(m1,m1,s1,s2)*q(1,n)**2 &
               +mu(m2,m2,s1,s2)*q(2,n)**2)
          
          ! Third-order contributions
          w(n)=w(n)-(1.0d0/6.0d0)*(tau(m1,s1,s2)*q(1,n)**3 &
               +tau(m2,s1,s2)*q(2,n)**3)
          
          ! Fourth-order contributions
          w(n)=w(n)-(1.0d0/24.0d0)*(xi(m1,s1,s2)*q(1,n)**4 &
               +xi(m2,s1,s2)*q(2,n)**4)

       enddo
          
    endif
       
    return
    
  end subroutine get_2mode_contrib_potential

!######################################################################

  subroutine get_2mode_contrib_dipole(d,q,ndat,m1,m2,s1,s2,c)

    use constants
    use channels
    use iomod
    use sysinfo
    use parameters
    use kdcglobal
    
    implicit none

    integer, intent(in)                     :: ndat,m1,m2,s1,s2,c
    integer                                 :: n
    real(dp), dimension(ndat)               :: d
    real(dp), dimension(2,ndat), intent(in) :: q

    ! Loop over points
    do n=1,ndat

       ! Zeroth-order contributions
       d(n)=d(n)-dip0(s1,s2,c)

       ! First-order contributions
       d(n)=d(n)-dip1(m1,s1,s2,c)*q(1,n)
       d(n)=d(n)-dip1(m2,s1,s2,c)*q(2,n)
       
       ! Second-order contributions
       d(n)=d(n)-0.5d0*dip2(m1,m1,s1,s2,c)*q(1,n)**2
       d(n)=d(n)-0.5d0*dip2(m2,m2,s1,s2,c)*q(2,n)**2

       ! Third-order contributions
       d(n)=d(n)-(1.0d0/6.0d0)*dip3(m1,s1,s2,c)*q(1,n)**3
       d(n)=d(n)-(1.0d0/6.0d0)*dip3(m2,s1,s2,c)*q(2,n)**3
       
       ! Fourth-order contributions
       d(n)=d(n)-(1.0d0/24.0d0)*dip4(m1,s1,s2,c)*q(1,n)**4
       d(n)=d(n)-(1.0d0/24.0d0)*dip4(m2,s1,s2,c)*q(2,n)**4
       
    enddo
    
    return
    
  end subroutine get_2mode_contrib_dipole
  
!######################################################################

  subroutine nmeq_bilinear(coeff,ndat,q,w)

    use constants
    use channels
    use iomod
    use utils
    
    implicit none

    integer, intent(in)                     :: ndat
    integer                                 :: i
    real(dp)                                :: coeff
    real(dp), dimension(2,ndat), intent(in) :: q
    real(dp), dimension(ndat), intent(in)   :: w
    real(dp)                                :: numer,denom
    
!----------------------------------------------------------------------
! Calculate the bi-linear coefficient
!----------------------------------------------------------------------
    numer=sum(w)

    denom=0.0d0
    do i=1,ndat
       denom=denom+q(1,i)*q(2,i)
    enddo

    coeff=numer/denom
    
    return
    
  end subroutine nmeq_bilinear

!######################################################################

  subroutine fill_coeffs2d_potential(coeff,m1,m2,s1,s2)

    use constants
    use channels
    use iomod
    use sysinfo
    use parameters
    use kdcglobal
    
    implicit none

    integer, intent(in)  :: m1,m2,s1,s2
    real(dp), intent(in) :: coeff
    
!----------------------------------------------------------------------
! Intrastate coupling coefficients
!----------------------------------------------------------------------
    if (s1.eq.s2) then
       gamma(m1,m2,s1)=coeff
       gamma(m2,m1,s1)=gamma(m1,m2,s1)
    endif 

!----------------------------------------------------------------------
! Interstate coupling coefficients
!----------------------------------------------------------------------
    if (s1.ne.s2) then
       mu(m1,m2,s1,s2)=coeff
       mu(m2,m1,s1,s2)=mu(m1,m2,s1,s2)
       mu(m2,m1,s2,s1)=mu(m1,m2,s1,s2)
    endif
    
    return
    
  end subroutine fill_coeffs2d_potential

!######################################################################

  subroutine fill_coeffs2d_potential_new(coeff,order,m1,m2,s1,s2)

    use constants
    use channels
    use iomod
    use sysinfo
    use parameters
    use kdcglobal
    
    implicit none

    integer, intent(in)                    :: order,m1,m2,s1,s2
    real(dp), intent(in), dimension(order) :: coeff
    real(dp)                               :: c1,c2,c12

!----------------------------------------------------------------------
! Intrastate coupling coefficients
!----------------------------------------------------------------------
    if (s1.eq.s2) then

       c1=gamma(m1,m1,s1)+freq(m1)/eh2ev
       c2=gamma(m2,m2,s1)+freq(m2)/eh2ev

       c12=2.0d0*coeff(2)

       print*,'here'
       stop
       
       gamma(m1,m2,s1)=0.5d0*(2.0d0*c12-c1-c2)

       gamma(m2,m1,s1)=gamma(m1,m2,s1)
    endif

!----------------------------------------------------------------------
! Interstate coupling coefficients
!----------------------------------------------------------------------
    !if (s1.ne.s2) then
    !   mu(m1,m2,s1,s2)=2.0d0*coeff(2)-mu(m1,m1,s1,s2)-mu(m2,m2,s1,s2)
    !   mu(m1,m2,s1,s2)=0.5d0*mu(m1,m2,s1,s2)
    !   mu(m2,m1,s1,s2)=mu(m1,m2,s1,s2)
    !   mu(m2,m1,s2,s1)=mu(m1,m2,s1,s2)
    !endif
    
    return
    
  end subroutine fill_coeffs2d_potential_new
    
!######################################################################

   subroutine fill_coeffs2d_dipole(coeff,m1,m2,s1,s2,c)

    use constants
    use channels
    use iomod
    use sysinfo
    use parameters
    use kdcglobal
    
    implicit none

    integer, intent(in)  :: m1,m2,s1,s2,c
    real(dp), intent(in) :: coeff

    dip2(m1,m2,s1,s2,c)=coeff
    dip2(m2,m1,s1,s2,c)=dip2(m1,m2,s1,s2,c)
    dip2(m2,m2,s2,s1,c)=dip2(m1,m2,s1,s2,c)
    
    return
    
  end subroutine fill_coeffs2d_dipole

!######################################################################

  function mask_pot_2mode(m1,m2,s1,s2) result(mask)
    
    use symmetry
    
    implicit none

    integer, intent(in) :: m1,m2,s1,s2
    integer             :: mask

!----------------------------------------------------------------------
! Determine the mask value for the given mode and state indices
!----------------------------------------------------------------------
    ! Intrastate coupling coefficient
    if (s1.eq.s2) mask=gamma_mask(m1,m2,s1)

    ! Interstate coupling coefficient
    if (s1.ne.s2) mask=mu_mask(m1,m2,s1,s2)
        
    return
    
  end function mask_pot_2mode
  
!######################################################################

  subroutine calc_rmsd

    use constants
    use channels
    use parameters
    use potfuncs
    use sysinfo
    use utils
    use kdcglobal
    
    implicit none

    integer                        :: m,m1,m2,s1,s2,s,n,ndat,indx,&
                                      ndat_tot
    real(dp), dimension(nsta,nsta) :: wmod,wact
    real(dp), dimension(nsta)      :: vmod,vact
    real(dp), dimension(nsta,nsta) :: eigvec
    
!----------------------------------------------------------------------
! Allocate and initialise the RMSD arrays
!----------------------------------------------------------------------
    allocate(rmsd1m(nmodes))
    allocate(rmsd2m(nmodes,nmodes))
    
    rmsd=0
    rmsd1m=0.0d0
    rmsd2m=0.0d0
    
!----------------------------------------------------------------------
! Contribution to the RMSD from the geometries with one mode
! displaced
!----------------------------------------------------------------------
    ! Loop over modes
    do m=1,nmodes

       ! Cycle if there are no points for the current mode
       ndat=ngeom1m(m)
       if (ndat.eq.0) cycle

       ! Loop over geometries with the current mode displaced
       do n=1,ndat

          ! File index
          indx=findx1m(m,n)
          
          ! Model potential (in eV)
          wmod=diabaticpot(qvec(:,indx))
          vmod=adiabaticpot(qvec(:,indx))
          
          ! Actual potential (in eV)
          do s1=1,nsta
             wact(s1,s1)=(diabpot(s1,s1,indx)-q0pot(s1))*eh2ev
          enddo
          do s1=1,nsta-1
             do s2=s1+1,nsta
                wact(s1,s2)=diabpot(s1,s2,indx)*eh2ev
                wact(s2,s1)=wact(s1,s2)
             enddo
          enddo
          call diag_matrix(diabpot(:,:,indx),vact,nsta)
          vact=(vact-q0pot(1))*eh2ev
          
          ! RMSD contribution
          do s=1,nsta
             rmsd1m(m)=rmsd1m(m)+(vmod(s)-vact(s))**2
          enddo
          
       enddo
       
    enddo

!----------------------------------------------------------------------
! Contribution to the RMSD from the geometries with two modes
! displaced
!----------------------------------------------------------------------
    ! Loop over pairs of modes
    do m1=1,nmodes-1
       do m2=m1+1,nmodes

          ! Cycle if there are no points for the current pair of modes
          ndat=ngeom2m(m1,m2)
          if (ndat.eq.0) cycle

          ! Loop over geometries with the current pair of modes
          ! displaced
          do n=1,ndat

             ! File index
             indx=findx2m(m1,m2,n)

             ! Model potential (in eV)
             wmod=diabaticpot(qvec(:,indx))
             vmod=adiabaticpot(qvec(:,indx))
             
             ! Actual potential (in eV)
             do s1=1,nsta
                wact(s1,s1)=(diabpot(s1,s1,indx)-q0pot(s1))*eh2ev
             enddo
             do s1=1,nsta-1
                do s2=s1+1,nsta
                   wact(s1,s2)=diabpot(s1,s2,indx)*eh2ev
                   wact(s2,s1)=wact(s1,s2)
                enddo
             enddo
             call diag_matrix(diabpot(:,:,indx),vact,nsta)
             vact=(vact-q0pot(1))*eh2ev
             
          enddo

          ! RMSD contribution
          do s=1,nsta
             rmsd2m(m1,m2)=rmsd2m(m1,m2)+(vmod(s)-vact(s))**2
          enddo
          
       enddo
    enddo

!----------------------------------------------------------------------
! Total RMSD
!----------------------------------------------------------------------
    ! Total number of displaced geometries
    ndat_tot=sum(ngeom1m(:))
    do m1=1,nmodes-1
       do m2=m1+1,nmodes
          ndat_tot=ndat_tot+ngeom2m(m1,m2)
       enddo
    enddo

    ! Total RMSD
    rmsd=sum(rmsd1m)
    do m1=1,nmodes-1
       do m2=m1+1,nmodes
          rmsd=rmsd+rmsd2m(m1,m2)
       enddo
    enddo
    rmsd=sqrt(rmsd/(ndat_tot*nsta))

!----------------------------------------------------------------------
! One-mode RMSDs
!----------------------------------------------------------------------
    do m=1,nmodes
       ndat=ngeom1m(m)
       if (ndat.eq.0) cycle
       rmsd1m(m)=sqrt(rmsd1m(m)/(ndat*nsta))
    enddo

!----------------------------------------------------------------------
! Two-mode RMSDs
!----------------------------------------------------------------------
    do m1=1,nmodes-1
       do m2=m1+1,nmodes
          ndat=ngeom2m(m1,m2)
          if (ndat.eq.0) cycle
          rmsd2m(m1,m2)=rmsd2m(m1,m2)/(ndat*nsta)
       enddo
    enddo
    
    return
    
  end subroutine calc_rmsd
    
!######################################################################
  
end module nmeqmod

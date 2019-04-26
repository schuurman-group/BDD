module fdmod

  implicit none

contains

!######################################################################
  
  subroutine get_coefficients_fd

    implicit none

!----------------------------------------------------------------------
! Determine the displacement information for each blockdiag log file
!----------------------------------------------------------------------
    call getdispinfo
  
!----------------------------------------------------------------------
! Calculate the coupling coefficients
!----------------------------------------------------------------------
    call calc_coefficients
    
    return

  end subroutine get_coefficients_fd

!######################################################################

  subroutine getdispinfo

    use constants
    use channels
    use iomod
    use sysinfo
    use symmetry
    use kdcglobal
    
    implicit none

    integer             :: n,m,m1,m2,ndisp,n1,n2,n3,n4,nmissing
    real(dp)            :: sum1,sum2
    real(dp), parameter :: thrsh=1e-4_dp

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(coefftyp(nfiles))
    coefftyp=0

    allocate(findx1m(nmodes,2))
    findx1m=0

    allocate(findx2m(nmodes,nmodes,4))
    findx2m=0
    
!----------------------------------------------------------------------
! Determine the displacement type for each file/geometry
!
! coefftyp(n) = 1 <-> linear or quadratic terms
!                     (kappa_a, lambda_a, gamma_aa, mu_aa)
!
! coefftyp(n) = 2 <-> bi-linear terms
!                     (gamma_ab, mu_ab)
!----------------------------------------------------------------------
    lbilinear=.false.

    ! Loop over displaced geometries
    do n=1,nfiles

       ! Determine the no. displaced modes
       ndisp=0
       do m=1,nmodes
          if (abs(qvec(m,n)).gt.thrsh) ndisp=ndisp+1
       enddo

       ! Set coefftyp(n)
       if (ndisp.eq.1) then
          coefftyp(n)=1
       else if (ndisp.eq.2) then
          coefftyp(n)=2
          lbilinear=.true.
       else
          errmsg='Too many displaced modes in file: '&
               //trim(bdfiles(n))
          call error_control
       endif
       
    enddo

!----------------------------------------------------------------------
! Set the mapping for file <-> mode displacement
!----------------------------------------------------------------------
    ! Linear and quadratic terms
    do n=1,nfiles

       ! Skip if not a file for linear and quadratic terms
       if (coefftyp(n).ne.1) cycle

       ! Determine which mode is displaced and in which direction
       do m=1,nmodes
          if (abs(qvec(m,n)).gt.thrsh) then
             if (qvec(m,n).lt.0.0d0) then
                ! Negatice displacement along mode m
                findx1m(m,1)=n
             else
                ! Positive displacement along mode m
                findx1m(m,2)=n
             endif
          endif
       enddo
       
    enddo

    ! Bi-linear terms
    if (lbilinear) then

       do n=1,nfiles

          ! Skip if not a file for bi-linear terms
          if (coefftyp(n).ne.2) cycle

          ! Determine which modes are displaced and in which directions
          do m1=1,nmodes-1
             do m2=m1+1,nmodes
                if (abs(qvec(m1,n)).gt.thrsh &
                     .and.abs(qvec(m2,n)).gt.thrsh) then
                   if (qvec(m1,n).gt.0.0d0.and.qvec(m2,n).gt.0.0d0) then
                      ! + +
                      findx2m(m1,m2,1)=n
                   else if (qvec(m1,n).lt.0.0d0.and.qvec(m2,n).lt.0.0d0) then
                      ! - -
                      findx2m(m1,m2,2)=n
                   else if (qvec(m1,n).gt.0.0d0.and.qvec(m2,n).lt.0.0d0) then
                      ! + -
                      findx2m(m1,m2,3)=n
                   else if (qvec(m1,n).lt.0.0d0.and.qvec(m2,n).gt.0.0d0) then
                      ! - +
                      findx2m(m1,m2,4)=n
                   endif
                endif
             enddo
          enddo
          
       enddo

    endif

!----------------------------------------------------------------------
! Check that consistent displacements were used
!----------------------------------------------------------------------
    ! Linear and quadratic terms
    do m=1,nmodes

       n1=findx1m(m,1)
       n2=findx1m(m,2)
       
       ! Cycle if we are missing information for this mode
       if (n1.eq.0.or.n2.eq.0) cycle

       ! Check that the positive and negative displacements are the
       ! same
       if (abs(abs(qvec(m,n1))-abs((qvec(m,n2)))).gt.thrsh) then
          write(errmsg,'(a,x,i3)') 'Inconsistent positive/negative &
               displacement sizes for mode',m
          call error_control
       endif
       
    enddo
    
    ! Bi-linear terms
    ! Write this part of the code later...
    if (lbilinear) then

       do m1=1,nmodes-1
          do m2=m1+1,nmodes

             ! Skip if there are no symmetry-allowed coefficients for
             ! the current mode pair
             if (cut_mask(m1,m2).eq.0) cycle
             
             n1=findx2m(m1,m2,1)
             n2=findx2m(m1,m2,2)
             n3=findx2m(m1,m2,3)
             n4=findx2m(m1,m2,4)

             ! Cycle if we are missing information for this pair
             ! of modes
             if (n1.eq.0.or.n2.eq.0.or.n3.eq.0.or.n4.eq.0) cycle

             ! Check that all displacements are the same for each mode
             sum1=(qvec(m1,n1)+qvec(m1,n2)+qvec(m1,n3)+qvec(m1,n4))/4.0d0
             sum2=(qvec(m2,n1)+qvec(m2,n2)+qvec(m2,n3)+qvec(m2,n4))/4.0d0
             if (abs(sum1).gt.thrsh) then
                write(errmsg,'(a,x,i0,x,a,x,i0,a,i0)') 'Inconsistent &
                     positive/negative displacement sizes for mode',m1,&
                     'for the mode pair',m1,'/',m2
                call error_control
             endif
             if (abs(sum2).gt.thrsh) then
                write(errmsg,'(a,x,i0,x,a,x,i0,a,i0)') 'Inconsistent &
                     positive/negative displacement sizes for mode',m2,&
                     'for the mode pair',m1,'/',m2
                call error_control
             endif

          enddo
       enddo
       
    endif

!----------------------------------------------------------------------
! Write out some information about missing data
!----------------------------------------------------------------------
    ! Linear and quadratic terms
    do m=1,nmodes

       n1=findx1m(m,1)
       n2=findx1m(m,2)
       
       ! Completely missing
       if (n1.eq.0.and.n2.eq.0) then
          write(ilog,'(/,2x,a,2x,i3)') 'Linear and quadratic terms will &
               not be computed for mode',m
          write(ilog,'(2x,a)') 'Reason: two missing files'
       endif

       ! One missing file
       if (n1.eq.0.and.n2.ne.0 &
            .or.n2.eq.0.and.n1.ne.0) then
          write(ilog,'(/,2x,a,2x,i3)') 'Linear and quadratic terms will &
               not be computed for mode',m
          write(ilog,'(2x,a)') 'Reason: one missing file'
       endif
       
    enddo

    ! Bi-linear terms
    ! Write this part of the code later...
    if (lbilinear) then

       do m1=1,nmodes-1
          do m2=m1+1,nmodes

             ! Skip if there are no symmetry-allowed coefficients for
             ! the current mode pair
             if (cut_mask(m1,m2).eq.0) cycle

             ! No. missing files
             nmissing=0
             n1=findx2m(m1,m2,1)
             if (n1.eq.0) nmissing=nmissing+1
             n2=findx2m(m1,m2,2)
             if (n2.eq.0) nmissing=nmissing+1
             n3=findx2m(m1,m2,3)
             if (n3.eq.0) nmissing=nmissing+1
             n4=findx2m(m1,m2,4)
             if (n4.eq.0) nmissing=nmissing+1
                          
             ! Missing files
             if (nmissing.ne.0) then
                write(ilog,'(/,2x,a,x,i0,x,a,x,i0)') 'Bi-linear terms &
                     will not be computed for modes',m1,'and',m2
                write(ilog,'(2x,a,x,i0,x,a)') &
                     'Reason:',nmissing,'missing file(s)'
             endif
                          
          enddo
       enddo
             
    endif
       
    return
    
  end subroutine getdispinfo

!######################################################################

  subroutine calc_coefficients

    use constants
    use channels
    use iomod
    use sysinfo
    use parameters
    use kdcglobal
    
    implicit none

    integer  :: m1,m2,s1,s2,n1,n2,n3,n4
    real(dp) :: dq1,dq2

!**********************************************************************
! Note that we are here assuming that the adiabatic and diabatic
! representations are being taken to be equal at the reference point.
!
! Specifically, in the calculation of the second-order terms, we are
! taking the Q0 adiabatic and diabatic potential to be equal.
!**********************************************************************
    
!----------------------------------------------------------------------
! First-order and quadratic terms: kappa_a, lambda_a, gamma_aa, mu_aa
!----------------------------------------------------------------------
    do m1=1,nmodes

       ! File index for the negative displacement
       n1=findx1m(m1,1)

       ! File index for the positive displacement
       n2=findx1m(m1,2)

       ! Cycle if we are missing information for this mode
       if (n1.eq.0.or.n2.eq.0) cycle

       ! Step size
       dq1=(abs(qvec(m1,n1))+abs(qvec(m1,n2)))/2.0d0
       
       ! kappa
       do s1=1,nsta
          kappa(m1,s1)=(diabpot(s1,s1,n2)-diabpot(s1,s1,n1))/(2.0d0*dq1)
       enddo

       ! lambda
       do s1=1,nsta-1
          do s2=s1+1,nsta
             lambda(m1,s1,s2)=(diabpot(s1,s2,n2)-diabpot(s1,s2,n1))/(2.0d0*dq1)
             lambda(m1,s2,s1)=lambda(m1,s1,s2)
          enddo
       enddo

       ! Quadratic gamma terms
       do s1=1,nsta
          gamma(m1,m1,s1)=(diabpot(s1,s1,n2)+diabpot(s1,s1,n1)&
               -2.0d0*q0pot(s1))/dq1**2
          gamma(m1,m1,s1)=gamma(m1,m1,s1)-freq(m1)/eh2ev
       enddo

       ! Quadratic mu terms
       do s1=1,nsta-1
          do s2=s1+1,nsta
             mu(m1,m1,s1,s2)=&
                  (diabpot(s1,s2,n2)+diabpot(s1,s2,n1))/dq1**2
             mu(m1,m1,s2,s1)=mu(m1,m1,s1,s2)
          enddo
       enddo
       
    enddo
    
!----------------------------------------------------------------------
! Bi-linear terms: gamma_ab and mu_ab
!----------------------------------------------------------------------
    if (lbilinear) then

       do m1=1,nmodes-1
          do m2=m1+1,nmodes

             ! File index for the + + displacement
             n1=findx2m(m1,m2,1)
             
             ! File index for the - - displacement
             n2=findx2m(m1,m2,2)

             ! File index for the + - displacement
             n3=findx2m(m1,m2,3)

             ! File index for the - + displacement
             n4=findx2m(m1,m2,4)

             ! Step sizes
             dq1=(abs(qvec(m1,n1))+abs(qvec(m1,n2)) &
                  +abs(qvec(m1,n3))+abs(qvec(m1,n4)))/4.0d0
             dq2=(abs(qvec(m2,n1))+abs(qvec(m2,n2)) &
                  +abs(qvec(m2,n3))+abs(qvec(m2,n4)))/4.0d0

             ! Bi-linear gamma terms
             do s1=1,nsta
                gamma(m1,m2,s1)=(diabpot(s1,s1,n1)+diabpot(s1,s1,n2) &
                     -diabpot(s1,s1,n3)-diabpot(s1,s1,n4))&
                     /(4.0d0*dq1*dq2)
                gamma(m2,m1,s1)=gamma(m1,m2,s1)
             enddo

             ! Bi-linear mu terms
             do s1=1,nsta-1
                do s2=s1+1,nsta
                   mu(m1,m2,s1,s2)=&
                        (diabpot(s1,s2,n1)+diabpot(s1,s2,n2) &
                        -diabpot(s1,s2,n3)-diabpot(s1,s2,n4))&
                        /(4.0d0*dq1*dq2)
                   mu(m2,m1,s1,s2)=mu(m1,m2,s1,s2)
                   mu(m1,m2,s2,s1)=mu(m1,m2,s1,s2)
                   mu(m2,m1,s2,s1)=mu(m1,m2,s1,s2)
                enddo
             enddo
             
          enddo
       enddo
             
    endif
       
!----------------------------------------------------------------------
! Convert all coupling coefficients to units of eV
!----------------------------------------------------------------------
    kappa=kappa*eh2ev
    lambda=lambda*eh2ev
    gamma=gamma*eh2ev
    mu=mu*eh2ev
    
    return
    
  end subroutine calc_coefficients
  
!######################################################################
  
end module fdmod

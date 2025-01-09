!######################################################################
! parinfo: routines to output some useful information about the QVC
!          Hamiltonian parameters to the log file
!######################################################################

module parinfo

  use constants
  
  implicit none

  save

  real(dp), allocatable :: efreq(:,:)
  real(dp), allocatable :: wkappa(:,:)
  real(dp), allocatable :: wlambda(:,:,:,:)
  
contains

!######################################################################  

  subroutine get_parinfo

    implicit none

!----------------------------------------------------------------------
! Allocate and initialise arrays
!----------------------------------------------------------------------
    call initialise_parinfo
    
!----------------------------------------------------------------------
! Calculate and output the effective frequencies
!----------------------------------------------------------------------
    call effective_freqs

!----------------------------------------------------------------------
! Calculate and output the effective frequency-weighted first-order
! coupling coefficients
!----------------------------------------------------------------------
    call weighted_pars

!----------------------------------------------------------------------
! Ouput an estimation of the relative spectroscopic importance of
! the normal modes
!----------------------------------------------------------------------
    call rank_modes
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    call finalise_parinfo
    
    return
    
  end subroutine get_parinfo

!######################################################################  

  subroutine initialise_parinfo

    use constants
    use sysinfo
    
    implicit none

    ! Effective frequencies
    allocate(efreq(nmodes,nsta))
    efreq=0.0d0

    ! Effective frequency-weighted first-order coupling coefficients
    allocate(wkappa(nmodes,nsta))
    wkappa=0.0d0
    allocate(wlambda(nmodes,nsta,nsta,2))
    wlambda=0.0d0
    
    return
    
  end subroutine initialise_parinfo

!######################################################################  

  subroutine finalise_parinfo

    use constants
    use sysinfo
    
    implicit none

    deallocate(efreq)
    deallocate(wkappa)
    deallocate(wlambda)
    
    return
    
  end subroutine finalise_parinfo
  
!######################################################################  

  subroutine effective_freqs

    use constants
    use channels
    use iomod
    use sysinfo
    use parameters
    use kdcglobal
    
    implicit none

    integer                          :: s,m,k
    real(dp)                         :: fac
    character(len=60)                :: string
    logical, dimension(nmodes,nsta)  :: imaginary
   
!----------------------------------------------------------------------
! Calculate the effective frequencies
!----------------------------------------------------------------------
    ! Flag for imaginary effective frequencies
    imaginary=.false.
    
    do s=1,nsta
       do m=1,nmodes
          fac=freq(m)**2+gamma(m,m,s)*freq(m)
          if (fac.lt.0.0d0) imaginary(m,s)=.true.
          efreq(m,s)=sqrt(abs(fac))
       enddo
    enddo

!----------------------------------------------------------------------
! Write the effective frequencies to the log file
!----------------------------------------------------------------------
    ! Section Header
    write(ilog,'(/,72a)') ('+',k=1,72)
    write(ilog,'(2x,a)') 'Effective Frequencies (eV)'
    write(ilog,'(72a)') ('+',k=1,72)

    ! Table header
    write(ilog,'(/,40a)') ('-',k=1,40)
    write(ilog,'(a)') ' Mode  |  State  | Effective Frequency'
    write(ilog,'(40a)') ('-',k=1,40)

    ! Table of effective frequencies
    do s=1,nsta
       do m=1,nmodes
          write(string,'(i3,3x,a,2x,i2,5x,a,2x,F6.4)') m,'|',s,'|',&
               efreq(m,s)
          if (imaginary(m,s)) string=trim(string)//' (Imaginary)'
          write(ilog,'(x,a)') trim(string)
       enddo
       write(ilog,'(72a)') ('-',k=1,40)
    enddo
    
    return
    
  end subroutine effective_freqs
    
!######################################################################

  subroutine weighted_pars

    use constants
    use channels
    use iomod
    use sysinfo
    use symmetry
    use parameters
    use kdcglobal
    
    implicit none

    integer :: m,s,s1,s2,k
    
!----------------------------------------------------------------------
! Calculate the effective frequency-weighted first-order coupling
! coefficients
!----------------------------------------------------------------------
    ! kappa
    do m=1,nmodes
       do s=1,nsta
          wkappa(m,s)=kappa(m,s)/efreq(m,s)
       enddo
    enddo

    ! lambda
    do m=1,nmodes
       do s1=1,nsta-1
          do s2=s1+1,nsta
             wlambda(m,s1,s2,1)=lambda(m,s1,s2)/efreq(m,s1)
             wlambda(m,s1,s2,2)=lambda(m,s1,s2)/efreq(m,s2)
             wlambda(m,s2,s1,1)=wlambda(m,s1,s2,1)
             wlambda(m,s2,s1,2)=wlambda(m,s1,s2,2)
          enddo
       enddo
    enddo

!----------------------------------------------------------------------
! Write the effective frequency-weighted first-order coupling
! coefficients to the log file
!----------------------------------------------------------------------
    ! Section Header
    write(ilog,'(/,72a)') ('+',k=1,72)
    write(ilog,'(2x,a)') 'Effective Frequency-Weighted Coupling &
         Coefficients'
    write(ilog,'(72a)') ('+',k=1,72)

    ! kappa
    write(ilog,'(/,40a)') ('-',k=1,40)
    write(ilog,'(a)') ' Mode  |  State  |  w-kappa'
    write(ilog,'(40a)') ('-',k=1,40)
    do s=1,nsta
       do m=1,nmodes
          if (kappa_mask(m,s).eq.0) cycle
          write(ilog,'(x,i3,3x,a,2x,i2,5x,a,2x,F7.4)') m,'|',s,'|',&
               wkappa(m,s)
       enddo
       write(ilog,'(44a)') ('-',k=1,40)
    enddo

    ! lambda, lower state weighting
    write(ilog,'(/,44a)') ('-',k=1,44)
    write(ilog,'(a)') ' Mode  |  State i  |  State j  |  wi-lambda'
    write(ilog,'(44a)') ('-',k=1,44)
    do s1=1,nsta-1
       do s2=s1+1,nsta
          do m=1,nmodes
             if (lambda_mask(m,s1,s2).eq.0) cycle
             write(ilog,'(x,i3,3x,a,2x,i2,7x,a,2x,i2,7x,a,2x,F7.4)') &
                  m,'|',s1,'|',s2,'|',wlambda(m,s1,s2,1)
          enddo
          write(ilog,'(44a)') ('-',k=1,44)
       enddo
    enddo

    ! lambda, upper state weighting
    write(ilog,'(/,44a)') ('-',k=1,44)
    write(ilog,'(a)') ' Mode  |  State i  |  State j  |  wj-lambda'
    write(ilog,'(44a)') ('-',k=1,44)
    do s1=1,nsta-1
       do s2=s1+1,nsta
          do m=1,nmodes
             if (lambda_mask(m,s1,s2).eq.0) cycle
             write(ilog,'(x,i3,3x,a,2x,i2,7x,a,2x,i2,7x,a,2x,F7.4)') &
                  m,'|',s1,'|',s2,'|',wlambda(m,s1,s2,2)
          enddo
          write(ilog,'(44a)') ('-',k=1,44)
       enddo
    enddo
    
    return
    
  end subroutine weighted_pars

!######################################################################

  subroutine rank_modes

    use constants
    use channels
    use iomod
    use sysinfo
    use symmetry
    use kdcglobal
    use utils
    
    implicit none

    integer                          :: m,s,s1,s2,k
    integer, dimension(nmodes)       :: indx
    real(dp), dimension(nmodes,nsta) :: maxpar
    real(dp)                         :: fwp
    real(dp), parameter              :: thrsh=1e-4_dp
    
!----------------------------------------------------------------------
! Determine which modes give rise to the largest effective
! frequency-weigted first-order coupling coefficients for each state
!----------------------------------------------------------------------
    maxpar=0.0d0

    ! kappa
    do m=1,nmodes
       do s=1,nsta
          if (kappa_mask(m,s).eq.0) cycle
          fwp=abs(wkappa(m,s))
          if (fwp.ge.maxpar(m,s)) maxpar(m,s)=fwp
       enddo
    enddo

    ! lambda
    do m=1,nmodes
       do s1=1,nsta-1
          do s2=s1+1,nsta
             if (lambda_mask(m,s1,s2).eq.0) cycle
             fwp=max(abs(wlambda(m,s1,s2,1)),abs(wlambda(m,s1,s2,2)))
             if (fwp.ge.maxpar(m,s1)) maxpar(m,s1)=fwp
             if (fwp.ge.maxpar(m,s2)) maxpar(m,s2)=fwp
          enddo
       enddo
    enddo

!----------------------------------------------------------------------    
! Estimation of the relative spectroscopic importance of the normal
! modes
!----------------------------------------------------------------------
    ! Section header
    write(ilog,'(/,72a)') ('+',k=1,72)
    write(ilog,'(x,a)') 'Estimation of the relative spectroscopic &
         importance of the normal modes'
    write(ilog,'(72a)') ('+',k=1,72)

    ! Table header
    write(ilog,'(/,30a)') ('-', k=1,30)
    write(ilog,'(x,a)') 'State  |  Mode  |  Max. Val.'
    write(ilog,'(30a)') ('-', k=1,30)

    ! Table
    do s=1,nsta
       call dsortindxa1('D',nmodes,maxpar(:,s),indx)
       do m=1,nmodes
          if (maxpar(indx(m),s).lt.thrsh) cycle
          write(ilog,'(x,i2,5x,a,2x,i3,3x,a,3x,F7.4)') &
               s,'|',indx(m),'|',maxpar(indx(m),s)
       enddo
       write(ilog,'(30a)') ('-', k=1,30)
    enddo

    return
    
  end subroutine rank_modes
  
!######################################################################
  
end module parinfo

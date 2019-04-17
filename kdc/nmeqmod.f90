module nmeqmod

  implicit none

contains

!######################################################################

  subroutine get_coefficients_nmeq

    implicit none

!----------------------------------------------------------------------
! Perform the fitting for the 1-mode terms
!----------------------------------------------------------------------
    call fit_1mode_terms
    
    return
    
  end subroutine get_coefficients_nmeq

!######################################################################

  subroutine fit_1mode_terms

    use constants
    use channels
    use iomod
    use sysinfo
    use symmetry
    use kdcglobal
    
    implicit none

    integer               :: m,n,s1,s2,ndat
    integer               :: order
    real(dp), allocatable :: coeff(:),q(:),w(:)
    
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
    allocate(coeff(order))
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
             call nmeq1d(order,coeff,ndat,q(1:ndat),w(1:ndat))

          enddo
       enddo
             
    enddo

    STOP
    
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
    use symmetry
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

  subroutine nmeq1d(order,coeff,npnts,x,y)

    use constants
    use channels
    use iomod
    
    implicit none

    integer, intent(in)                    :: order,npnts
    real(dp), dimension(order)             :: coeff
    real(dp), intent(in), dimension(npnts) :: x,y

    
    
    return
    
  end subroutine nmeq1d
  
!######################################################################
  
end module nmeqmod

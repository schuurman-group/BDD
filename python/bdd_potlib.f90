!######################################################################
! bdd_potlib: Fortran wrapper subroutines for f2py interface to the
!             BDD potential evaluation routines
!######################################################################

!######################################################################
! init: read the binary parameter file produced by KDC
!######################################################################

  subroutine init(filename)

    use constants
    use sysinfo
    use parameters
    use symmetry
    use iomod
    use pltdata

    implicit none

    character(len=*), intent(in) :: filename
    integer                      :: iunit

!----------------------------------------------------------------------
! Open the parameter file
!----------------------------------------------------------------------
    call freeunit(iunit)
    open(iunit,file=filename,form='unformatted',status='old',err=999)

!----------------------------------------------------------------------
! System dimensions
!----------------------------------------------------------------------
    read(iunit) nmodes
    read(iunit) nsta
    read(iunit) ncoo
    read(iunit) natm

!----------------------------------------------------------------------
! No. ab initio diabatic potential values
!----------------------------------------------------------------------
    read(iunit) ndat

!----------------------------------------------------------------------
! Order of the one-mode expansions
!----------------------------------------------------------------------
    read(iunit) order1

!----------------------------------------------------------------------
! Diabatic dipole flag
!----------------------------------------------------------------------
    read(iunit) ldip

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    ! Reference geometry
    allocate(xcoo0(ncoo))
    allocate(atnum(natm))

    ! Coordinate transformations
    allocate(nmcoo(nmodes,ncoo))
    allocate(coonm(ncoo,nmodes))

    ! Model diabatic potential arrays
    allocate(e0(nsta))
    allocate(freq(nmodes))
    allocate(coeff1(nmodes,nsta,nsta,order1))
    allocate(coeff2(nmodes,nmodes,nsta,nsta))
    allocate(coeff1_mask(nmodes,nsta,nsta,order1))
    allocate(coeff2_mask(nmodes,nmodes,nsta,nsta))

    ! Ab initio diabatic potential arrays
    allocate(wq0(nsta))
    allocate(wdisp(nsta,nsta,ndat))
    allocate(qvec(nmodes,ndat))

    ! Diabatic dipole arrays
    if (ldip) then
       ! Coefficients
       allocate(dip0(nsta,nsta,3))
       allocate(dip1(nmodes,nsta,nsta,3))
       allocate(dip2(nmodes,nmodes,nsta,nsta,3))
       allocate(dip3(nmodes,nsta,nsta,3))
       allocate(dip4(nmodes,nsta,nsta,3))
       ! Masks
       allocate(dip1_mask(nmodes,nsta,nsta,3))
       allocate(dip2_mask(nmodes,nmodes,nsta,nsta,3))
       allocate(dip3_mask(nmodes,nsta,nsta,3))
       allocate(dip4_mask(nmodes,nsta,nsta,3))
       ! Ab initio dipole values
       allocate(ddisp(nsta,nsta,3,ndat))
    endif

!----------------------------------------------------------------------
! Reference geometry
!----------------------------------------------------------------------
    read(iunit) xcoo0
    read(iunit) atnum

!----------------------------------------------------------------------
! Coordinate transformations
!----------------------------------------------------------------------
    read(iunit) nmcoo
    read(iunit) coonm

!----------------------------------------------------------------------
! Vertical excitation energies
!----------------------------------------------------------------------
    read(iunit) e0

!----------------------------------------------------------------------
! Frequencies
!----------------------------------------------------------------------
    read(iunit) freq

!----------------------------------------------------------------------
! Coupling coefficients
!----------------------------------------------------------------------
    read(iunit) coeff1
    read(iunit) coeff2

!----------------------------------------------------------------------
! Masks
!----------------------------------------------------------------------
    read(iunit) coeff1_mask
    read(iunit) coeff2_mask

!----------------------------------------------------------------------
! Ab initio diabatic potential values
!----------------------------------------------------------------------
    read(iunit) wq0
    read(iunit) wdisp
    read(iunit) qvec

!----------------------------------------------------------------------
! Diabatic dipole surfaces
!----------------------------------------------------------------------
    if (ldip) then
       ! Coefficients
       read(iunit) dip0
       read(iunit) dip1
       read(iunit) dip2
       read(iunit) dip3
       read(iunit) dip4
       ! Masks
       read(iunit) dip1_mask
       read(iunit) dip2_mask
       read(iunit) dip3_mask
       read(iunit) dip4_mask
       ! Ab initio diabatic dipoles
       read(iunit) ddisp
    endif

!----------------------------------------------------------------------
! Close the parameter file
!----------------------------------------------------------------------
    close(iunit)

    return

999 continue
    write(6,'(/,2x,a,/)') 'Error opening file: '//trim(filename)
    stop

  end subroutine init

!######################################################################
! get_dimensions: return system dimensions to Python
!######################################################################

  subroutine get_dimensions(nmodes_out,nsta_out,ndat_out,ldip_out)

    use sysinfo
    use pltdata

    implicit none

    integer, intent(out) :: nmodes_out,nsta_out,ndat_out
    logical, intent(out) :: ldip_out

    !f2py intent(out) nmodes_out,nsta_out,ndat_out,ldip_out

    nmodes_out=nmodes
    nsta_out=nsta
    ndat_out=ndat
    ldip_out=ldip

    return

  end subroutine get_dimensions

!######################################################################
! get_e0: return vertical excitation energies
!######################################################################

  subroutine get_e0(nsta_in,e0_out)

    use constants
    use parameters

    implicit none

    integer, intent(in)                    :: nsta_in
    real(dp), intent(out)                  :: e0_out(nsta_in)

    !f2py intent(in) nsta_in
    !f2py intent(out) e0_out
    !f2py depend(nsta_in) e0_out

    e0_out=e0

    return

  end subroutine get_e0

!######################################################################
! get_freq: return normal mode frequencies
!######################################################################

  subroutine get_freq(nmodes_in,freq_out)

    use constants
    use sysinfo

    implicit none

    integer, intent(in)                    :: nmodes_in
    real(dp), intent(out)                  :: freq_out(nmodes_in)

    !f2py intent(in) nmodes_in
    !f2py intent(out) freq_out
    !f2py depend(nmodes_in) freq_out

    freq_out=freq

    return

  end subroutine get_freq

!######################################################################
! calc_adiab_1d: compute adiabatic potentials along a 1D cut
!######################################################################

  subroutine calc_adiab_1d(mode,nsta_in,npnts,qi,qf,qgrid,vpot)

    use constants
    use sysinfo
    use potfuncs

    implicit none

    integer, intent(in)   :: mode,nsta_in,npnts
    real(dp), intent(in)  :: qi,qf
    real(dp), intent(out) :: qgrid(npnts)
    real(dp), intent(out) :: vpot(npnts,nsta_in)

    !f2py intent(in) mode,nsta_in,npnts,qi,qf
    !f2py intent(out) qgrid,vpot
    !f2py depend(npnts) qgrid
    !f2py depend(npnts,nsta_in) vpot

    integer               :: i
    real(dp)              :: dq
    real(dp)              :: q(nmodes)

    dq=(qf-qi)/(npnts-1)

    do i=1,npnts
       q=0.0d0
       q(mode)=qi+(i-1)*dq
       qgrid(i)=q(mode)
       vpot(i,1:nsta)=adiabaticpot(q)
    enddo

    return

  end subroutine calc_adiab_1d

!######################################################################
! calc_adiab_diag: compute adiabatic potentials along a diagonal cut
!######################################################################

  subroutine calc_adiab_diag(mode1,mode2,nsta_in,npnts,qi,qf, &
       qgrid,vpot)

    use constants
    use sysinfo
    use potfuncs

    implicit none

    integer, intent(in)   :: mode1,mode2,nsta_in,npnts
    real(dp), intent(in)  :: qi,qf
    real(dp), intent(out) :: qgrid(npnts)
    real(dp), intent(out) :: vpot(npnts,nsta_in)

    !f2py intent(in) mode1,mode2,nsta_in,npnts,qi,qf
    !f2py intent(out) qgrid,vpot
    !f2py depend(npnts) qgrid
    !f2py depend(npnts,nsta_in) vpot

    integer               :: i
    real(dp)              :: dq,coord
    real(dp)              :: q(nmodes)

    dq=(qf-qi)/(npnts-1)

    do i=1,npnts
       q=0.0d0
       coord=qi+(i-1)*dq
       q(mode1)=coord/sqrt(2.0d0)
       q(mode2)=coord/sqrt(2.0d0)
       qgrid(i)=coord
       vpot(i,1:nsta)=adiabaticpot(q)
    enddo

    return

  end subroutine calc_adiab_diag

!######################################################################
! calc_diab_1d: compute diabatic potentials along a 1D cut
!######################################################################

  subroutine calc_diab_1d(mode,nsta_in,npnts,qi,qf,qgrid,vpot)

    use constants
    use sysinfo
    use potfuncs

    implicit none

    integer, intent(in)   :: mode,nsta_in,npnts
    real(dp), intent(in)  :: qi,qf
    real(dp), intent(out) :: qgrid(npnts)
    real(dp), intent(out) :: vpot(npnts,nsta_in)

    !f2py intent(in) mode,nsta_in,npnts,qi,qf
    !f2py intent(out) qgrid,vpot
    !f2py depend(npnts) qgrid
    !f2py depend(npnts,nsta_in) vpot

    integer                        :: i,s
    real(dp)                       :: dq
    real(dp)                       :: q(nmodes)
    real(dp), dimension(nsta,nsta) :: wmat

    dq=(qf-qi)/(npnts-1)

    do i=1,npnts
       q=0.0d0
       q(mode)=qi+(i-1)*dq
       qgrid(i)=q(mode)
       wmat=diabaticpot(q)
       do s=1,nsta
          vpot(i,s)=wmat(s,s)
       enddo
    enddo

    return

  end subroutine calc_diab_1d

!######################################################################
! calc_diab_diag: compute diabatic potentials along a diagonal cut
!######################################################################

  subroutine calc_diab_diag(mode1,mode2,nsta_in,npnts,qi,qf, &
       qgrid,vpot)

    use constants
    use sysinfo
    use potfuncs

    implicit none

    integer, intent(in)   :: mode1,mode2,nsta_in,npnts
    real(dp), intent(in)  :: qi,qf
    real(dp), intent(out) :: qgrid(npnts)
    real(dp), intent(out) :: vpot(npnts,nsta_in)

    !f2py intent(in) mode1,mode2,nsta_in,npnts,qi,qf
    !f2py intent(out) qgrid,vpot
    !f2py depend(npnts) qgrid
    !f2py depend(npnts,nsta_in) vpot

    integer                        :: i,s
    real(dp)                       :: dq,coord
    real(dp)                       :: q(nmodes)
    real(dp), dimension(nsta,nsta) :: wmat

    dq=(qf-qi)/(npnts-1)

    do i=1,npnts
       q=0.0d0
       coord=qi+(i-1)*dq
       q(mode1)=coord/sqrt(2.0d0)
       q(mode2)=coord/sqrt(2.0d0)
       qgrid(i)=coord
       wmat=diabaticpot(q)
       do s=1,nsta
          vpot(i,s)=wmat(s,s)
       enddo
    enddo

    return

  end subroutine calc_diab_diag

!######################################################################
! calc_diabcp_1d: compute a diabatic coupling element along a 1D cut
!######################################################################

  subroutine calc_diabcp_1d(mode,s1,s2,npnts,qi,qf,qgrid,vcoup)

    use constants
    use sysinfo
    use potfuncs

    implicit none

    integer, intent(in)   :: mode,s1,s2,npnts
    real(dp), intent(in)  :: qi,qf
    real(dp), intent(out) :: qgrid(npnts)
    real(dp), intent(out) :: vcoup(npnts)

    !f2py intent(in) mode,s1,s2,npnts,qi,qf
    !f2py intent(out) qgrid,vcoup
    !f2py depend(npnts) qgrid,vcoup

    integer                        :: i
    real(dp)                       :: dq
    real(dp)                       :: q(nmodes)
    real(dp), dimension(nsta,nsta) :: wmat

    dq=(qf-qi)/(npnts-1)

    do i=1,npnts
       q=0.0d0
       q(mode)=qi+(i-1)*dq
       qgrid(i)=q(mode)
       wmat=diabaticpot(q)
       vcoup(i)=wmat(s1,s2)
    enddo

    return

  end subroutine calc_diabcp_1d

!######################################################################
! calc_diabcp_diag: compute a diabatic coupling along a diagonal cut
!######################################################################

  subroutine calc_diabcp_diag(mode1,mode2,s1,s2,npnts,qi,qf, &
       qgrid,vcoup)

    use constants
    use sysinfo
    use potfuncs

    implicit none

    integer, intent(in)   :: mode1,mode2,s1,s2,npnts
    real(dp), intent(in)  :: qi,qf
    real(dp), intent(out) :: qgrid(npnts)
    real(dp), intent(out) :: vcoup(npnts)

    !f2py intent(in) mode1,mode2,s1,s2,npnts,qi,qf
    !f2py intent(out) qgrid,vcoup
    !f2py depend(npnts) qgrid,vcoup

    integer                        :: i
    real(dp)                       :: dq,coord
    real(dp)                       :: q(nmodes)
    real(dp), dimension(nsta,nsta) :: wmat

    dq=(qf-qi)/(npnts-1)

    do i=1,npnts
       q=0.0d0
       coord=qi+(i-1)*dq
       q(mode1)=coord/sqrt(2.0d0)
       q(mode2)=coord/sqrt(2.0d0)
       qgrid(i)=coord
       wmat=diabaticpot(q)
       vcoup(i)=wmat(s1,s2)
    enddo

    return

  end subroutine calc_diabcp_diag

!######################################################################
! calc_dip_1d: compute diabatic dipole along a 1D cut
!######################################################################

  subroutine calc_dip_1d(mode,s1,s2,npnts,qi,qf,qgrid,vdip)

    use constants
    use sysinfo
    use potfuncs

    implicit none

    integer, intent(in)   :: mode,s1,s2,npnts
    real(dp), intent(in)  :: qi,qf
    real(dp), intent(out) :: qgrid(npnts)
    real(dp), intent(out) :: vdip(npnts,3)

    !f2py intent(in) mode,s1,s2,npnts,qi,qf
    !f2py intent(out) qgrid,vdip
    !f2py depend(npnts) qgrid,vdip

    integer               :: i
    real(dp)              :: dq
    real(dp)              :: q(nmodes)

    dq=(qf-qi)/(npnts-1)

    do i=1,npnts
       q=0.0d0
       q(mode)=qi+(i-1)*dq
       qgrid(i)=q(mode)
       vdip(i,:)=diabdip(q,s1,s2)
    enddo

    return

  end subroutine calc_dip_1d

!######################################################################
! calc_dip_diag: compute diabatic dipole along a diagonal cut
!######################################################################

  subroutine calc_dip_diag(mode1,mode2,s1,s2,npnts,qi,qf,qgrid,vdip)

    use constants
    use sysinfo
    use potfuncs

    implicit none

    integer, intent(in)   :: mode1,mode2,s1,s2,npnts
    real(dp), intent(in)  :: qi,qf
    real(dp), intent(out) :: qgrid(npnts)
    real(dp), intent(out) :: vdip(npnts,3)

    !f2py intent(in) mode1,mode2,s1,s2,npnts,qi,qf
    !f2py intent(out) qgrid,vdip
    !f2py depend(npnts) qgrid,vdip

    integer               :: i
    real(dp)              :: dq,coord
    real(dp)              :: q(nmodes)

    dq=(qf-qi)/(npnts-1)

    do i=1,npnts
       q=0.0d0
       coord=qi+(i-1)*dq
       q(mode1)=coord/sqrt(2.0d0)
       q(mode2)=coord/sqrt(2.0d0)
       qgrid(i)=coord
       vdip(i,:)=diabdip(q,s1,s2)
    enddo

    return

  end subroutine calc_dip_diag

!######################################################################
! get_abinit_pot: extract ab initio potential data for a given cut
!######################################################################

  subroutine get_abinit_pot(mode,mode2,ldiag,surftyp,s1,s2, &
       ndat_in,nsta_in,nab,qab,vab)

    use constants
    use sysinfo
    use parameters
    use pltdata

    implicit none

    integer, intent(in)   :: mode,mode2,surftyp,s1,s2
    integer, intent(in)   :: ndat_in,nsta_in
    logical, intent(in)   :: ldiag
    integer, intent(out)  :: nab
    real(dp), intent(out) :: qab(ndat_in+1)
    real(dp), intent(out) :: vab(max(3,nsta_in),max(3,nsta_in),ndat_in+1)

    !f2py intent(in) mode,mode2,ldiag,surftyp,s1,s2,ndat_in,nsta_in
    !f2py intent(out) nab,qab,vab
    !f2py depend(ndat_in) qab
    !f2py depend(ndat_in,nsta_in) vab

    integer               :: m,n,ndisp,s,ss1,ss2,c
    integer               :: mindx(2)
    integer               :: e2,error
    real(dp), parameter   :: thrsh=1e-4_dp
    real(dp)              :: v(nsta)
    real(dp)              :: w(nsta,nsta)
    real(dp)              :: work(3*nsta)

    nab=0
    qab=0.0d0
    vab=0.0d0

!----------------------------------------------------------------------
! Q=0 point (reference geometry)
!----------------------------------------------------------------------
    nab=1
    qab(1)=0.0d0

    if (surftyp.eq.1) then
       do s=1,nsta
          vab(s,s,1)=e0(s)
       enddo
    else if (surftyp.eq.2.or.surftyp.eq.4) then
       do s=1,nsta
          vab(s,s,1)=e0(s)
       enddo
    else if (surftyp.eq.3) then
       if (ldip) then
          do c=1,3
             vab(c,c,1)=dip0(s1,s2,c)
          enddo
       endif
    endif

!----------------------------------------------------------------------
! Displaced geometry points
!----------------------------------------------------------------------
    do n=1,ndat

       ndisp=0
       do m=1,nmodes
          if (abs(qvec(m,n)).gt.thrsh) then
             ndisp=ndisp+1
             if (ndisp.le.2) mindx(ndisp)=m
          endif
       enddo

       if (ldiag) then
          if (ndisp.ne.2) cycle
          if (mindx(1).ne.mode.or.mindx(2).ne.mode2) cycle
       else
          if (ndisp.ne.1) cycle
          if (mindx(1).ne.mode) cycle
       endif

       nab=nab+1

       if (ldiag) then
          qab(nab)=sqrt(2.0d0)*qvec(mode,n)
       else
          qab(nab)=qvec(mode,n)
       endif

       if (surftyp.eq.1) then
          e2=3*nsta
          w=wdisp(:,:,n)
          call dsyev('V','U',nsta,w,nsta,v,work,e2,error)
          do s=1,nsta
             vab(s,s,nab)=e0(s)+(v(s)-wq0(s))*eh2ev
          enddo
       endif

       if (surftyp.eq.2.or.surftyp.eq.4) then
          do s=1,nsta
             vab(s,s,nab)=e0(s)+(wdisp(s,s,n)-wq0(s))*eh2ev
          enddo
          do ss1=1,nsta-1
             do ss2=ss1+1,nsta
                vab(ss1,ss2,nab)=wdisp(ss1,ss2,n)*eh2ev
                vab(ss2,ss1,nab)=vab(ss1,ss2,nab)
             enddo
          enddo
       endif

       if (surftyp.eq.3) then
          if (ldip) then
             do c=1,3
                vab(c,c,nab)=ddisp(s1,s2,c,n)
             enddo
          endif
       endif

    enddo

    return

  end subroutine get_abinit_pot

!######################################################################

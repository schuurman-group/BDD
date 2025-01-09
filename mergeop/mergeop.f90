!######################################################################
! mergeop: a program to merge two KDC Hamiltonians
!######################################################################

program mergeop

  use channels
  use iomod
  use opermod
  
  implicit none

!----------------------------------------------------------------------
! Read the command line arguments
!----------------------------------------------------------------------
  call rdinp_mergeop

!----------------------------------------------------------------------
! Read the parameter binary files
!----------------------------------------------------------------------
  call rdbinfiles

!----------------------------------------------------------------------
! Merge the models
!----------------------------------------------------------------------
  call merge_models

!----------------------------------------------------------------------
! Write the merged operator file
!----------------------------------------------------------------------
  call freeunit(iop)
  open(iop,file='merged.op',form='formatted',status='unknown')
  call wroper
  
contains

!######################################################################  

  subroutine rdinp_mergeop

    use constants
    use channels
    use iomod
    use mergeglobal
    
    implicit none

!----------------------------------------------------------------------
! Set default values
!----------------------------------------------------------------------
    ! Parameter file names
    abinA=''
    abinB=''
    
!----------------------------------------------------------------------
! Exit if the incorrect number of arguments have been given
!----------------------------------------------------------------------
    if (iargc() /= 2) then
       write(6,'(/,2x,a,/)') 'Incorrect number of arguments given:' &
            // ' two expected'
       STOP
    endif

!----------------------------------------------------------------------
! Get the parameter binary file names
!----------------------------------------------------------------------
    call getarg(1,abinA)
    call getarg(2,abinB)

!----------------------------------------------------------------------
! Sanity checks
!----------------------------------------------------------------------
    ! The two files cannot be the same
    if (abinA == abinB) then
       errmsg='Error: the two files cannot be the same'
       call error_control
    endif
        
    return
    
  end subroutine rdinp_mergeop

!######################################################################

  subroutine rdbinfiles

    use constants
    use channels
    use iomod
    use mergeglobal
    
    implicit none

    ! Binary file unit numbers
    integer               :: ibinA,ibinB

    ! Various arrays that are only needed here
    real(dp), allocatable :: xcoo0(:),nmcoo(:,:),coonm(:,:)
    integer, allocatable  :: atnum(:)
    real(dp), allocatable :: wq0A(:),wq0B(:),wq0(:)
    real(dp), allocatable :: wdispA(:,:,:),wdispB(:,:,:)
    real(dp), allocatable :: qvecA(:,:),qvecB(:,:)
    real(dp), allocatable :: ddispA(:,:,:,:),ddispB(:,:,:,:)

    ! Everything else
    integer               :: i,ndatA,ndatB
    real(dp)              :: emin,eminA,eminB
    
!----------------------------------------------------------------------
! Open the parameter files
!----------------------------------------------------------------------
    call freeunit(ibinA)
    open(ibinA,file=abinA,form='unformatted',status='old')

    call freeunit(ibinB)
    open(ibinB,file=abinB,form='unformatted',status='old')

!----------------------------------------------------------------------
! System dimensions
!----------------------------------------------------------------------
    read(ibinA) nmodesA
    read(ibinA) nstaA
    read(ibinA) ncooA
    read(ibinA) natmA
    
    read(ibinB) nmodesB
    read(ibinB) nstaB
    read(ibinB) ncooB
    read(ibinB) natmB

!----------------------------------------------------------------------
! Check that the dimensions are compatible
!----------------------------------------------------------------------
    if (nmodesA /= nmodesB) then
       errmsg='Error: unequal normal mode coordinate dimensions'
       call error_control
    endif

    if (ncooA /= ncooB) then
       errmsg='Error: unequal Cartesian coordinate dimensions'
       call error_control
    endif

    if (natmA /= natmB) then
       errmsg='Error: unequal numbers of atoms'
       call error_control
    endif
    
!----------------------------------------------------------------------
! No. ab initio diabatic potential values
!----------------------------------------------------------------------
    read(ibinA) ndatA
    read(ibinB) ndatB

!----------------------------------------------------------------------
! Order of the one-mode expansions
!----------------------------------------------------------------------
    read(ibinA) order1A
    read(ibinB) order1B

!----------------------------------------------------------------------
! Check that the one-mode expansion orders are equal
!----------------------------------------------------------------------
    if (order1A /= order1B) then
       errmsg='Error: unequal one-mode expansion orders'
       call error_control
    endif
    
!----------------------------------------------------------------------
! Diabatic dipole flag
!----------------------------------------------------------------------
    read(ibinA) ldipA
    read(ibinB) ldipB

!----------------------------------------------------------------------
! Check that the dipole fitting flags are equal
!----------------------------------------------------------------------
    if (ldipA .neqv. ldipB) then
       errmsg='Error: unequal dipole fitting flags'
       call error_control
    endif
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    ! Reference geometry
    allocate(xcoo0(ncooA))
    allocate(atnum(natmA))

    ! Coordinate transformations
    allocate(nmcoo(nmodesA,ncooA))
    allocate(coonm(ncooA,nmodesA))
    
    ! Model diabatic potential arrays
    allocate(e0A(nstaA))
    allocate(freqA(nmodesA))
    allocate(coeff1A(nmodesA,nstaA,nstaA,order1A))
    allocate(coeff2A(nmodesA,nmodesA,nstaA,nstaA))
    allocate(coeff1_maskA(nmodesA,nstaA,nstaA,order1A))
    allocate(coeff2_maskA(nmodesA,nmodesA,nstaA,nstaA))
    allocate(e0B(nstaB))
    allocate(freqB(nmodesB))
    allocate(coeff1B(nmodesB,nstaB,nstaB,order1B))
    allocate(coeff2B(nmodesB,nmodesB,nstaB,nstaB))
    allocate(coeff1_maskB(nmodesB,nstaB,nstaB,order1B))
    allocate(coeff2_maskB(nmodesB,nmodesB,nstaB,nstaB))

    ! Ab initio diabatic potential arrays
    allocate(wq0A(nstaA))    
    allocate(wdispA(nstaA,nstaA,ndatA))
    allocate(qvecA(nmodesA,ndatA))
    allocate(wq0B(nstaB))
    allocate(wdispB(nstaB,nstaB,ndatB))
    allocate(qvecB(nmodesB,ndatB))

    ! Diabatic dipole arrays
    if (ldipA) then
       ! Coefficients
       allocate(dip0A(nstaA,nstaA,3))
       allocate(dip1A(nmodesA,nstaA,nstaA,3))
       allocate(dip2A(nmodesA,nmodesA,nstaA,nstaA,3))
       allocate(dip3A(nmodesA,nstaA,nstaA,3))
       allocate(dip4A(nmodesA,nstaA,nstaA,3))
       allocate(dip0B(nstaB,nstaB,3))
       allocate(dip1B(nmodesB,nstaB,nstaB,3))
       allocate(dip2B(nmodesB,nmodesB,nstaB,nstaB,3))
       allocate(dip3B(nmodesB,nstaB,nstaB,3))
       allocate(dip4B(nmodesB,nstaB,nstaB,3))
       ! Masks
       allocate(dip1_maskA(nmodesA,nstaA,nstaA,3))
       allocate(dip2_maskA(nmodesA,nmodesA,nstaA,nstaA,3))
       allocate(dip3_maskA(nmodesA,nstaA,nstaA,3))
       allocate(dip4_maskA(nmodesA,nstaA,nstaA,3))
       allocate(dip1_maskB(nmodesB,nstaB,nstaB,3))
       allocate(dip2_maskB(nmodesB,nmodesB,nstaB,nstaB,3))
       allocate(dip3_maskB(nmodesB,nstaB,nstaB,3))
       allocate(dip4_maskB(nmodesB,nstaB,nstaB,3))
       ! Ab initio dipole values
       allocate(ddispA(nstaA,nstaA,3,ndatA))
       allocate(ddispB(nstaB,nstaB,3,ndatB))
    endif

!----------------------------------------------------------------------
! Reference geometry
!----------------------------------------------------------------------
    read(ibinA) xcoo0
    read(ibinA) atnum

    read(ibinB) xcoo0
    read(ibinB) atnum

!----------------------------------------------------------------------
! Coordinate transformations
!----------------------------------------------------------------------
    read(ibinA) nmcoo
    read(ibinA) coonm

    read(ibinB) nmcoo
    read(ibinB) coonm

!----------------------------------------------------------------------
! Vertical excitation energies
!----------------------------------------------------------------------
    read(ibinA) e0A

    read(ibinB) e0B

!----------------------------------------------------------------------
! Frequencies
!----------------------------------------------------------------------
    read(ibinA) freqA
    
    read(ibinB) freqB

!----------------------------------------------------------------------
! Check that the frequencies match
!----------------------------------------------------------------------
    do i=1,nmodesA
       if (abs(freqA(i)-freqB(i)) > 1e-6_dp) then
          errmsg='Error: different frequencies found'
          call error_control
       endif
    enddo
    
!----------------------------------------------------------------------
! Coupling coefficients
!----------------------------------------------------------------------
    read(ibinA) coeff1A
    read(ibinA) coeff2A

    read(ibinB) coeff1B
    read(ibinB) coeff2B

!----------------------------------------------------------------------
! Masks
!----------------------------------------------------------------------
    read(ibinA) coeff1_maskA
    read(ibinA) coeff2_maskA

    read(ibinB) coeff1_maskB
    read(ibinB) coeff2_maskB
    
!----------------------------------------------------------------------
! Ab inito diabatic potential values
!----------------------------------------------------------------------
    read(ibinA) wq0A
    read(ibinA) wdispA
    read(ibinA) qvecA

    read(ibinB) wq0B
    read(ibinB) wdispB
    read(ibinB) qvecB

!----------------------------------------------------------------------
! Diabatic dipole surfaces
!----------------------------------------------------------------------
    if (ldipA) then

       ! Coefficients
       read(ibinA) dip0A
       read(ibinA) dip1A
       read(ibinA) dip2A
       read(ibinA) dip3A
       read(ibinA) dip4A

       read(ibinB) dip0B
       read(ibinB) dip1B
       read(ibinB) dip2B
       read(ibinB) dip3B
       read(ibinB) dip4B
       
       ! Masks
       read(ibinA) dip1_maskA
       read(ibinA) dip2_maskA
       read(ibinA) dip3_maskA
       read(ibinA) dip4_maskA

       read(ibinB) dip1_maskB
       read(ibinB) dip2_maskB
       read(ibinB) dip3_maskB
       read(ibinB) dip4_maskB
       
       ! Ab initio diabatic dipoles
       read(ibinA) ddispA

       read(ibinB) ddispB
       
    endif

!----------------------------------------------------------------------
! Take all energies relative to the lowest-lying state across both
! models
!----------------------------------------------------------------------
    allocate(wq0(nstaA+nstaB))
    wq0(1:nstaA)=wq0A
    wq0(nstaA+1:nstaA+nstaB)=wq0B

    emin=minval(wq0)

    e0A=(wq0A-emin)*eh2ev
    e0B=(wq0B-emin)*eh2ev

!----------------------------------------------------------------------
! Close the parameter files
!----------------------------------------------------------------------
    close(ibinA)
    
    close(ibinB)
    
    return
    
  end subroutine rdbinfiles
  
!######################################################################

  subroutine merge_models

    use constants
    use iomod
    use sysinfo
    use mergeglobal
    use symmetry
    use parameters
    use kdcglobal, only: ldipfit
    
    implicit none

!----------------------------------------------------------------------
! Dimensions
!----------------------------------------------------------------------
    nmodes=nmodesA
    nsta=nstaA+nstaB

!----------------------------------------------------------------------
! One-mode term expansion order
!----------------------------------------------------------------------
    order1=order1A
    
!----------------------------------------------------------------------
! Frequencies
!----------------------------------------------------------------------
    allocate(freq(nmodes))
    freq=freqA

!----------------------------------------------------------------------
! Vertical excitation energies
!----------------------------------------------------------------------
    allocate(e0(nsta))
    e0(1:nstaA)=e0A
    e0(nstaA+1:nsta)=e0B

!----------------------------------------------------------------------
! One-mode coupling coefficients
!----------------------------------------------------------------------
    allocate(coeff1(nmodes,nsta,nsta,order1))
    coeff1=0.0d0
    coeff1(:,1:nstaA,1:nstaA,:)=coeff1A
    coeff1(:,nstaA+1:nsta,nstaA+1:nsta,:)=coeff1B

    allocate(coeff1_mask(nmodes,nsta,nsta,order1))
    coeff1_mask=0
    coeff1_mask(:,1:nstaA,1:nstaA,:)=coeff1_maskA
    coeff1_mask(:,nstaA+1:nsta,nstaA+1:nsta,:)=coeff1_maskB

!----------------------------------------------------------------------
! Two-mode coupling coefficients
!----------------------------------------------------------------------
    allocate(coeff2(nmodes,nmodes,nsta,nsta))
    coeff2=0.0d0
    coeff2(:,:,1:nstaA,1:nstaA)=coeff2A
    coeff2(:,:,nstaA+1:nsta,nstaA+1:nsta)=coeff2B

    allocate(coeff2_mask(nmodes,nmodes,nsta,nsta))
    coeff2_mask=0
    coeff2_mask(:,:,1:nstaA,1:nstaA)=coeff2_maskA
    coeff2_mask(:,:,nstaA+1:nsta,nstaA+1:nsta)=coeff2_maskB

!----------------------------------------------------------------------
! Dipole expansion coefficients
!----------------------------------------------------------------------
    ldipfit=ldipA

    if (ldipfit) then

       allocate(dip0(nsta,nsta,3))
       dip0=0.0d0
       dip0(1:nstaA,1:nstaA,:)=dip0A
       dip0(nstaA+1:nsta,nstaA+1:nsta,:)=dip0B

       allocate(dip1(nmodes,nsta,nsta,3))
       dip1=0.0d0
       dip1(:,1:nstaA,1:nstaA,:)=dip1A
       dip1(:,nstaA+1:nsta,nstaA+1:nsta,:)=dip1B

       allocate(dip2(nmodes,nmodes,nsta,nsta,3))
       dip2=0.0d0
       dip2(:,:,1:nstaA,1:nstaA,:)=dip2A
       dip2(:,:,nstaA+1:nsta,nstaA+1:nsta,:)=dip2B

       allocate(dip3(nmodes,nsta,nsta,3))
       dip3=0.0d0
       dip3(:,1:nstaA,1:nstaA,:)=dip3A
       dip3(:,nstaA+1:nsta,nstaA+1:nsta,:)=dip3B

       allocate(dip4(nmodes,nsta,nsta,3))
       dip4=0.0d0
       dip4(:,1:nstaA,1:nstaA,:)=dip4A
       dip4(:,nstaA+1:nsta,nstaA+1:nsta,:)=dip4B
       
    endif
    
    return
    
  end subroutine merge_models
  
!######################################################################
  
end program mergeop

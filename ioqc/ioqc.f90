module ioqc

  implicit none

  save

  integer :: idum
  logical :: ldum
  
contains

!######################################################################
! freqtype: determines the quantum chemistry program used for the
!           frequency calculation, and sets freqtyp accordingly:
!          
!           freqtyp = 1 <-> G98
!                     2 <-> CFOUR
!                     3 <-> Hessian file
!                     4 <-> Turbomole, aoforce
!######################################################################

  subroutine freqtype

    use constants
    use sysinfo
    use iomod

    implicit none
    
    integer            :: unit
    character(len=120) :: string
    logical            :: dir

!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    freqtyp=-1

!----------------------------------------------------------------------
! Determine the quantum chemistry program used for the frequency
! calculation
!----------------------------------------------------------------------
    if (isg98(freqfile)) then
       ! G98
       freqtyp=1
    else if (iscfour(freqfile)) then
       ! CFOUR
       freqtyp=2
    else if (ishessian(freqfile)) then
       ! Hessian
       freqtyp=3
    else if (isaoforce(freqfile)) then
       ! Turbomole, aoforce
       freqtyp=4
    endif

!----------------------------------------------------------------------
! Check that the quantum chemistry program used is supported
!----------------------------------------------------------------------
    if (freqtyp.eq.-1) then
       errmsg='The quantum chemistry program used for the &
            frequency calculation is not supported.'
       call error_control
    endif

    return

  end subroutine freqtype

!######################################################################

  function isg98(filename) result(found)

    use constants
    use iomod

    implicit none

    integer            :: unit
    character(len=*)   :: filename
    character(len=120) :: string
    logical            :: found,dir

!----------------------------------------------------------------------
! First determine whether freqfile is actually a file.
! This is necessary as for certain programs, the name of a directory
! will actually be passed instead
!----------------------------------------------------------------------
    inquire(file=trim(filename)//'/.',exist=dir)

    if (dir) then
       found=.false.
       return
    endif

!----------------------------------------------------------------------
! Open file
!----------------------------------------------------------------------
    call freeunit(unit)
    open(unit,file=filename,form='formatted',status='old')

!----------------------------------------------------------------------
! Check whether the calculation was performed using G98
!----------------------------------------------------------------------
    read(unit,'(a)') string
    if (index(string,'Entering Gaussian System').ne.0) then
       found=.true.
    else
       found=.false.
    endif

!----------------------------------------------------------------------
! Close file
!----------------------------------------------------------------------
    close(unit)

    return

  end function isg98

!######################################################################

  function iscfour(filename) result(found)

    use constants
    use iomod

    implicit none

    integer            :: unit
    character(len=*)   :: filename
    character(len=120) :: string
    logical            :: found,dir

!----------------------------------------------------------------------
! First determine whether freqfile is actually a file.
! This is necessary as for certain programs, the name of a directory
! will actually be passed instead
!----------------------------------------------------------------------
    inquire(file=trim(filename)//'/.',exist=dir)

    if (.not.dir) then
       found=.false.
       return
    endif

!----------------------------------------------------------------------
! Check to see if the cfour FCM file exists
!----------------------------------------------------------------------
    inquire(file=trim(filename)//'/FCM',exist=found)

    return

  end function iscfour

!######################################################################
  
  function ishessian(filename) result(found)

    use constants
    use iomod

    implicit none

    integer            :: unit
    character(len=*)   :: filename
    character(len=120) :: string
    logical            :: found,dir

!----------------------------------------------------------------------
! First determine whether freqfile is actually a file.
! This is necessary as for certain programs, the name of a directory
! will actually be passed instead
!----------------------------------------------------------------------
    inquire(file=trim(filename)//'/.',exist=dir)

    if (dir) then
       found=.false.
       return
    endif

!----------------------------------------------------------------------
! Open file
!----------------------------------------------------------------------
    call freeunit(unit)
    open(unit,file=filename,form='formatted',status='old')

!----------------------------------------------------------------------
! Check whether the file is a Hessian file
!----------------------------------------------------------------------
    read(unit,'(a)') string
    if (index(string,'Hessian').ne.0) then
       found=.true.
    else
       found=.false.
    endif

!----------------------------------------------------------------------
! Close file
!----------------------------------------------------------------------
    close(unit)

    return

  end function ishessian

!######################################################################

  function isaoforce (filename) result(found)

    use constants
    use iomod

    implicit none

    integer            :: unit
    character(len=*)   :: filename
    character(len=120) :: string
    logical            :: found,dir

!----------------------------------------------------------------------
! First determine whether freqfile is actually a file.
! This is necessary as for certain programs, the name of a directory
! will actually be passed instead
!----------------------------------------------------------------------
    inquire(file=trim(filename)//'/.',exist=dir)

    if (dir) then
       found=.false.
       return
    endif

!----------------------------------------------------------------------
! Open file
!----------------------------------------------------------------------
    call freeunit(unit)
    open(unit,file=filename,form='formatted',status='old')

!----------------------------------------------------------------------
! Check whether the calculation was performed using the aoforce
! module in Turbomole
!----------------------------------------------------------------------
    found=.false.

5   read(unit,'(a)',end=10) string
    if (index(string,'a o f o r c e - program').ne.0) then
       found=.true.
    else
       goto 5
    endif

10  continue
    
!----------------------------------------------------------------------
! Close file
!----------------------------------------------------------------------
    close(unit)
    
    return
    
  end function isaoforce
  
!######################################################################

  subroutine getnatm

    use sysinfo
    use iomod
    
    implicit none

    if (freqtyp.eq.1) then
       ! G98
       call getnatm_g98
    else if (freqtyp.eq.2) then
       ! CFOUR
       errmsg='The interface for CFOUR still needs writing!'
       call error_control
    else if (freqtyp.eq.3) then
       ! Hessian file
       errmsg='The interface for Hessian files still needs writing!'       
       call error_control
    else if (freqtyp.eq.4) then
       ! Turbomole, aoforce
       call getnatm_aoforce
    endif
    
    return
    
  end subroutine getnatm

!######################################################################

  subroutine getnatm_g98

    use constants
    use channels
    use sysinfo
    use iomod
    
    implicit none

    integer            :: unit,i
    character(len=120) :: string
    
!----------------------------------------------------------------------
! Open the normal mode file
!----------------------------------------------------------------------
    call freeunit(unit)
    open(unit,file=freqfile,form='formatted',status='old')

!----------------------------------------------------------------------
! Determine the number of atoms
!----------------------------------------------------------------------
    ! Read to the Cartesian coordinates section (standard orientation)
5   read(unit,'(a)',end=100) string
    if (index(string,'Standard orientation').eq.0) goto 5
    do i=1,4
       read(unit,*)
    enddo

    ! Determine the number of atoms
    natm=0
10  read(unit,'(a)',end=200) string
    if (string(2:3).ne.'--') then
       natm=natm+1
       goto 10
    endif

!----------------------------------------------------------------------
! Close the normal mode file
!----------------------------------------------------------------------
    close(unit)
    
    return

100 continue
    errmsg='The Cartesian coordinates could not be found in: '&
         //trim(freqfile)
    call error_control
    
200 continue
    errmsg='End of file reached reading the Cartesian coordinate &
         section in: '//trim(freqfile)
    call error_control
    
  end subroutine getnatm_g98

!######################################################################

  subroutine getnatm_aoforce

    use constants
    use channels
    use sysinfo
    use iomod
    
    implicit none

    integer            :: unit,i
    character(len=120) :: string
    
!----------------------------------------------------------------------
! Open the normal mode file
!----------------------------------------------------------------------
    call freeunit(unit)
    open(unit,file=freqfile,form='formatted',status='old')

!----------------------------------------------------------------------
! Determine the number of atoms
!----------------------------------------------------------------------
    ! Read to the Cartesian coordinates section
5   read(unit,'(a)',end=100) string
    if (index(string,'atomic coordinates').eq.0) goto 5

    ! Determine the number of atoms
    natm=-2
10  read(unit,'(a)',end=200) string
    if (index(string,'center of nuclear mass').eq.0) then
       natm=natm+1
       goto 10
    endif

!----------------------------------------------------------------------
! Close the normal mode file
!----------------------------------------------------------------------
    close(unit)
    
    return

100 continue
    errmsg='The Cartesian coordinates could not be found in: '&
         //trim(freqfile)
    call error_control
    
200 continue
    errmsg='End of file reached reading the Cartesian coordinate &
         section in: '//trim(freqfile)
    call error_control
    
  end subroutine getnatm_aoforce
  
!######################################################################

    subroutine getxcoo0

    use constants
    use sysinfo
    use iomod

    implicit none

    real(dp), dimension(ncoo) :: xcoo
    
!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    ldum=.false.
    idum=0

!----------------------------------------------------------------------
! Read in the reference Cartesian coordinates
!----------------------------------------------------------------------
    if (freqtyp.eq.1) then
       ! G98
       call getxcoo_g98(xcoo,freqfile)
    else if (freqtyp.eq.2) then
       ! CFOUR
       call getxcoo_cfour(xcoo,freqfile)
    else if (freqtyp.eq.3) then
       ! Hessian
       call getxcoo_hessian(xcoo,freqfile)
    else if (freqtyp.eq.4) then
       ! Turbomole, aoforce
       call getxcoo_aoforce(xcoo,freqfile)
    endif

!----------------------------------------------------------------------
! Fill in the xcoo0 array
!----------------------------------------------------------------------
    xcoo0=xcoo
    
    return

  end subroutine getxcoo0

!######################################################################
  
  subroutine getxcoo_g98(xcoo,filename)

    use constants
    use sysinfo
    use iomod

    implicit none
  
    integer                   :: unit,i,j
    real(dp), dimension(ncoo) :: xcoo
    character(len=*)          :: filename
    character(len=120)        :: string

!----------------------------------------------------------------------
! Open file
!----------------------------------------------------------------------
    call freeunit(unit)
    open(unit,file=filename,form='formatted',status='old')

!----------------------------------------------------------------------
! Read the Cartesian coordinates
!----------------------------------------------------------------------
5   read(unit,'(a)',end=999) string
    if (index(string,'Standard orientation:').eq.0) goto 5

    do i=1,4
       read(unit,*)
    enddo

    do i=1,natm
       read(unit,'(14x,i2,18x,3F12.6)') atnum(i),(xcoo(j), j=i*3-2,i*3)
       atlbl(i)=num2lbl(atnum(i))
       mass(i*3-2:i*3)=num2mass(atnum(i))
    enddo

    ! Convert to Bohr
    xcoo=xcoo*ang2bohr

!----------------------------------------------------------------------
! Close file
!----------------------------------------------------------------------
    close(unit)

    return

999 continue
    errmsg='The Cartesian coordinates could not be found in: '&
         //trim(filename)
    call error_control

  end subroutine getxcoo_g98

!######################################################################

  subroutine getxcoo_cfour(xcoo,filename)

    use constants
    use iomod
    use sysinfo

    implicit none

    integer                   :: unit,i,n
    real(dp), dimension(ncoo) :: xcoo
    logical                   :: found
    character(len=*)          :: filename
    character(len=120)        :: string
    character(len=2)          :: atmp

!----------------------------------------------------------------------
! First check to make sure that the cfour.log file exists
!----------------------------------------------------------------------
    inquire(file=trim(filename)//'/cfour.log',exist=found)

    if (.not.found) then
       errmsg='The CFOUR log file is assumed to be named cfour.log. &
            This file could not be found in: '//trim(filename)
       call error_control
    endif

!----------------------------------------------------------------------
! Open the CFOUR log file
!----------------------------------------------------------------------
    call freeunit(unit)
    open(unit,file=trim(filename)//'/cfour.log',form='formatted',&
         status='old')

!----------------------------------------------------------------------
! Read the reference Cartesian coordinates from the CFOUR log file
!----------------------------------------------------------------------
5   read(unit,'(a)',end=999) string
    if (index(string,'Coordinates (in bohr)').eq.0) goto 5

    do i=1,2
       read(unit,*)
    enddo

    n=0
10  read(unit,'(a)') string
    if (index(string,'---').eq.0) then
       if (string(6:6).ne.'X') then
          n=n+1
          read(string,'(5x,a2,8x,i2,3x,3(F15.8))') &
               atlbl(n),atnum(n),(xcoo(i), i=n*3-2,n*3)
          mass(n*3-2:n*3)=num2mass(atnum(n))
       else
          ldum=.true.
          idum=n+1
       endif
       goto 10
    endif

!----------------------------------------------------------------------
! Close the CFOUR log file
!----------------------------------------------------------------------
    close(unit)

    return

999 continue
    errmsg='The Cartesian coordinates could not be found in: '&
         //trim(filename)//'/cfour.log'
    call error_control

  end subroutine getxcoo_cfour

!######################################################################

  subroutine getxcoo_hessian(xcoo,filename)

    use constants
    use iomod
    use sysinfo

    implicit none

    integer                   :: unit,i,j
    real(dp), dimension(ncoo) :: xcoo
    character(len=*)          :: filename

!----------------------------------------------------------------------
! Open file
!----------------------------------------------------------------------
    call freeunit(unit)
    open(unit,file=filename,form='formatted',status='old')

!----------------------------------------------------------------------
! Read the Cartesian coordinates
!----------------------------------------------------------------------
    read(unit,*)
    do i=1,natm
       read(unit,*) atlbl(i),(xcoo(j), j=i*3-2,i*3)
       mass(i*3-2:i*3)=lbl2mass(atlbl(i))
    enddo

    ! Convert to Bohr
    xcoo=xcoo*ang2bohr
    
!----------------------------------------------------------------------
! Close file
!----------------------------------------------------------------------
    close(unit)

    return

  end subroutine getxcoo_hessian

!######################################################################
  
  subroutine getxcoo_aoforce(xcoo,filename)

    use constants
    use sysinfo
    use iomod

    implicit none
  
    integer                   :: unit,i,j
    real(dp), dimension(ncoo) :: xcoo
    real(dp)                  :: ftmp
    character(len=*)          :: filename
    character(len=120)        :: string

!----------------------------------------------------------------------
! Open file
!----------------------------------------------------------------------
    call freeunit(unit)
    open(unit,file=filename,form='formatted',status='old')

!----------------------------------------------------------------------
! Read the Cartesian coordinates
!----------------------------------------------------------------------
    ! Read to the coordinates section
5   read(unit,'(a)',end=999) string
    if (index(string,'atomic coordinates').eq.0) goto 5

    ! Read the coordinates (in Bohr)
    do i=1,natm
       read(unit,'(x,3(2x,F12.8),4x,a2,9x,F6.4)') &
            (xcoo(j), j=i*3-2,i*3),atlbl(i),ftmp
       atnum(i)=int(ftmp)
       mass(i*3-2:i*3)=num2mass(atnum(i))
       call uppercase(atlbl(i)(1:1))
    enddo

!----------------------------------------------------------------------
! Close file
!----------------------------------------------------------------------
    close(unit)

    return

999 continue
    errmsg='The Cartesian coordinates could not be found in: '&
         //trim(filename)
    call error_control

  end subroutine getxcoo_aoforce

!######################################################################

  subroutine uppercase(string)
      
    implicit none
    
    integer*8    ::  i,j,length
    character(*) :: string
    
    length=len(string)
    
    do j = 1,length
       if(string(j:j).ge."a".and.string(j:j).le."z")&
            string(j:j)=achar(iachar(string(j:j))-32)
    enddo
    
    return

  end subroutine uppercase
  
!######################################################################

  function num2lbl(num) result(lbl)

    use constants
    use iomod

    implicit none

    integer          :: num
    character(len=2) :: lbl

    if (num.eq.1) then
       lbl='H'
    else if (num.eq.6) then
       lbl='C'
    else if (num.eq.7) then
       lbl='N'
    else if (num.eq.8) then
       lbl='O'
    else
       write(errmsg,'(a,x,i2,x,a)') 'The atomic number',num,&
            'is not supported. See function num2lbl.'
       call error_control
    endif

    return

  end function num2lbl

!######################################################################

  function num2mass(num) result(mass)

    use constants
    use iomod

    implicit none

    integer  :: num
    real(dp) :: mass

    if (num.eq.1) then
       mass=1.00794d0
    else if (num.eq.6) then
       mass=12.0107d0
    else if (num.eq.7) then
       mass=14.0067d0
    else if (num.eq.8) then
       mass=15.9994d0
    else
       write(errmsg,'(a,x,i2,x,a)') 'The atomic number',num,&
            'is not supported. See function num2mass.'
       call error_control
    endif

    return

  end function num2mass

!######################################################################
  
  function lbl2mass(lbl) result(mass)

    use constants
    use iomod

    implicit none

    real(dp)         :: mass
    character(len=*) :: lbl

    if (lbl.eq.'H') then
       mass=1.00794d0
    else if (lbl.eq.'C') then
       mass=12.0107d0
    else if (lbl.eq.'N') then
       mass=14.0067d0
    else if (lbl.eq.'O') then
       mass=15.9994d0
    else
       errmsg='The atom type '//trim(lbl)//' is not supported.'&
            //' See function lbl2mass'
    endif

    return

  end function lbl2mass

!######################################################################

  subroutine getdim

    use constants
    use sysinfo
    use iomod
    
    implicit none

!----------------------------------------------------------------------
! Determine the no. atoms
!----------------------------------------------------------------------
    call getnatm

!----------------------------------------------------------------------
! Set ncoo
!----------------------------------------------------------------------
    ncoo=3*natm
    
!----------------------------------------------------------------------
! Allocate arrays associated with the reference geometry
!----------------------------------------------------------------------
    ! Reference coordinates
    allocate(xcoo0(ncoo))
    xcoo0=0.0d0

    ! Atomic numbers
    allocate(atnum(natm))
    atnum=0

    ! Atom labels
    allocate(atlbl(natm))
    atlbl=''

    ! Atomic masses
    allocate(mass(ncoo))
    mass=0.0d0
    
    return
    
  end subroutine getdim
  
!######################################################################

  subroutine getnmodes

    use constants
    use sysinfo
    use iomod

    implicit none

    integer                  :: i,error,nzero
    real(dp), dimension(3)   :: iteig
    real(dp), dimension(3,3) :: itold
    real(dp), dimension(9)   :: work
    real(dp), parameter      :: thrsh=1e-10_dp
    
!-----------------------------------------------------------------------
! Calculate the moment of inertia tensor
!-----------------------------------------------------------------------
    itold=0.0d0
    do i=1,natm
       itold(1,1)=itold(1,1)+mass(i*3)*xcoo0(i*3-1)**2+xcoo0(i*3)**2
       itold(2,2)=itold(2,2)+mass(i*3)*xcoo0(i*3-2)**2+xcoo0(i*3)**2
       itold(3,3)=itold(3,3)+mass(i*3)*xcoo0(i*3-2)**2+xcoo0(i*3-1)**2
       itold(1,2)=itold(1,2)-mass(i*3)*xcoo0(i*3-2)*xcoo0(i*3-1)
       itold(1,3)=itold(1,3)-mass(i*3)*xcoo0(i*3-2)*xcoo0(i*3)
       itold(2,3)=itold(2,3)-mass(i*3)*xcoo0(i*3-1)*xcoo0(i*3)
    enddo
    itold(2,1)=itold(1,2)
    itold(3,1)=itold(1,3)
    itold(3,2)=itold(2,3)

!-----------------------------------------------------------------------
! Diagonalise the moment of inertia tensor
!-----------------------------------------------------------------------
    call dsyev('V','U',3,itold,3,iteig,work,9,error)

   if (error.ne.0) then
      errmsg='Diagonalisation of the moment of inertia tensor in &
           getnmodes failed'
      call error_control
   endif

!-----------------------------------------------------------------------
! Determine the number of normal modes
!-----------------------------------------------------------------------
   nzero=0
   do i=1,3
      if (abs(iteig(i)).lt.thrsh) nzero=nzero+1
   enddo

   if (nzero.eq.0) then
      ! Non-linear molecule
      nmodes=ncoo-6
   else if (nzero.eq.1) then
      ! Linear molecule
      nmodes=ncoo-5      
   else
      errmsg='Something has gone terribly wrong in getnmodes...'
      call error_control
   endif

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
   ! Transformation matrices
    allocate(nmcoo(ncoo,nmodes))
    allocate(coonm(nmodes,ncoo))
    nmcoo=0.0d0
    coonm=0.0d0

   ! Normal mode labels
   allocate(nmlab(nmodes))
   nmlab=''
   
   ! Normal mode frequencies
   allocate(freq(nmodes))
   freq=0.0d0
   
   return
    
  end subroutine getnmodes

!######################################################################

  subroutine getmodes
  
    use constants
    use sysinfo
    use iomod
    
    implicit none
    
    integer :: i,j

!----------------------------------------------------------------------
! Read the normal modes from file
!----------------------------------------------------------------------
    if (freqtyp.eq.1) then
       ! G98
       call getmodes_g98
    else if (freqtyp.eq.2) then
       ! CFOUR
       call getmodes_cfour
    else if (freqtyp.eq.3) then
       ! HESSIAN
       call getmodes_hessian
    else if (freqtyp.eq.4) then
       ! Turbomole, aoforce
       call getmodes_aoforce
    endif

!----------------------------------------------------------------------
! Orthogonalise the normal mode vectors
!----------------------------------------------------------------------
    call nmortho

!----------------------------------------------------------------------
! Write the normal modes to file for inspection
!----------------------------------------------------------------------
    call wrmodes
    
    return
    
  end subroutine getmodes

!######################################################################

  subroutine getmodes_g98

    use constants
    use sysinfo
    use iomod
    
    implicit none
    
    integer                          :: unit
    integer                          :: i,j,n,nm,nrd,ncurr
    real(dp), dimension(ncoo,ncoo)   :: matrix
    real(dp), dimension(10)          :: xdum
    real(dp)                         :: s2
    logical                          :: found
    character(len=120)               :: string
    character(len=16), dimension(10) :: buffer

!----------------------------------------------------------------------
! Open the frequency file
!----------------------------------------------------------------------
    call freeunit(unit)
    open(unit,file=freqfile,form='formatted',status='old')

!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    matrix=0.0d0

!----------------------------------------------------------------------
! Get the frequencies
!----------------------------------------------------------------------
100 continue

    read(unit,'(a)',end=999) string

    if (string(2:30) .ne. 'Harmonic frequencies (cm**-1)') go to 100

    nrd = ((nmodes-1)/5) + 1

110 continue
    read(unit,'(a)',end=900) string
    if (string(8:22) .ne. 'Frequencies ---') go to 110

    found=.true.
    backspace(unit)
    backspace(unit)
    backspace(unit)

    nm=0
    do n=1,nrd
       ncurr=min(5,nmodes-5*(n-1))
       read(unit,*)
       read(unit,'(26x,5a10)') (buffer(j),j=1,ncurr)
       read(unit,'(23x,5f10.4)',err=900) (xdum(j),j=1,ncurr)       
       do j=1,ncurr
          nm=nm+1
          nmlab(nm)=buffer(j)
          freq(nm)=xdum(j)
       enddo
120    continue
       read(unit,'(a)',end=999) string
       if (string(2:20) .ne. 'Coord Atom Element:') go to 120
       do j=1,ncoo
          read(unit,'(a)',err=900) string
          read(string,20) (xdum(i),i=1,ncurr)
          nm=nm-ncurr
          do i=1,ncurr
             nm=nm+1
             matrix(j,nm)=xdum(i)
          enddo
       enddo
    enddo
20  format(23x,5f10.5)

!----------------------------------------------------------------------
! Close the frequency file
!----------------------------------------------------------------------
    close(unit)

!-----------------------------------------------------------------------
! Scale to mass-weighted x to obtain true normal mode vectors
!
! N.B. This has to be done as Gaussian prints the normal mode vectors
!      in terms of non-mass-weighted Cartesians
!-----------------------------------------------------------------------
    do nm=1,nmodes
       do j=1,ncoo
          matrix(j,nm)=matrix(j,nm)*sqrt(mass(j))
       enddo
    enddo

!-----------------------------------------------------------------------
! Normalise the normal mode vectors
!-----------------------------------------------------------------------
    do i=1,nmodes
       matrix(:,i)=matrix(:,i) &
            /sqrt(dot_product(matrix(:,i),matrix(:,i)))
    enddo

!-----------------------------------------------------------------------
! Convert frequencies to eV
!-----------------------------------------------------------------------
    do i=1,nmodes
       freq(i)=freq(i)*invcm2ev
    enddo

!-----------------------------------------------------------------------
! Save the normal modes in the nmcoo array
!-----------------------------------------------------------------------
    nmcoo=matrix(1:ncoo,1:nmodes)

    return

900 continue
    errmsg='Incorrect normal mode format. Use IOP(7/8=11)'
    call error_control

999 continue
    errmsg='Problem reading the normal modes in: '//trim(freqfile)
    call error_control

  end subroutine getmodes_g98

!######################################################################

  subroutine getmodes_cfour

    use constants
    use iomod
    use sysinfo

    implicit none

    integer                            :: unit,i,j,k,lim1,lim2
    real(dp), dimension(ncoo,ncoo)     :: hess
    real(dp), dimension(ncoo+3,ncoo+3) :: dumhess
    real(dp), dimension(nmodes)        :: cffreq
    character(len=120)                 :: string


!**********************************************************************
! Note that cfour outputs the normal mode vectors to only a few
! decimal places. Therefore, we construct the normal modes ourself
! from the Hessian, which is written to a decent level of precision
! in the FCMFINAL file.
!**********************************************************************


!----------------------------------------------------------------------
! Read the Hessian from the FCMFINAL file
!----------------------------------------------------------------------
    ! Open the FCMFINAL file
    call freeunit(unit)
    open(unit,file=trim(freqfile)//'/FCMFINAL',form='formatted',&
         status='old')

    ! Read the FCMFINAL file
    read(unit,*)

    if (ldum) then
       ! Dummy atom present
       ! (1) Read the Hessian including the dummy atom
       do i=1,ncoo+3
          do j=1,natm+1
             read(unit,'(3F20.10)') (dumhess(i,k), k=j*3-2,j*3)
          enddo
       enddo
       ! (2) Remove the dummy atom contributions
       lim1=idum*3-2
       lim2=idum*3
       do i=1,ncoo+3
          if (i.ge.lim1.and.i.le.lim2) cycle
          do j=1,ncoo+3
             if (j.ge.lim1.and.j.le.lim2) cycle
             if (i.lt.lim1.and.j.lt.lim1) then
                hess(i,j)=dumhess(i,j)
             else if (i.lt.lim1.and.j.ge.lim1) then
                hess(i,j-3)=dumhess(i,j)
             else if (i.ge.lim1.and.j.lt.lim1) then
                hess(i-3,j)=dumhess(i,j)
             else if (i.ge.lim1.and.j.ge.lim1) then
                hess(i-3,j-3)=dumhess(i,j)
             endif
          enddo
       enddo
    else
       ! No dummy atom
       do i=1,ncoo
          do j=1,natm
             read(unit,'(3F20.10)') (hess(i,k), k=j*3-2,j*3)
          enddo
       enddo
    endif

    ! Close the FCMFINAL file
    close(unit)

!----------------------------------------------------------------------
! Calculate the normal mode vectors from the Hessian
!----------------------------------------------------------------------
    call hess2nm(hess)

!----------------------------------------------------------------------
! Consistency check: read the CFOUR log file frequencies and make sure
! that our frequencies match
!----------------------------------------------------------------------    
    call freeunit(unit)
    open(unit,file=trim(freqfile)//'/cfour.log')

5   read(unit,'(a)',end=888) string
    if (index(string,'rotational projection').eq.0) goto 5

    do i=1,7
       read(unit,*)
    enddo
       
    do i=1,nmodes
       read(unit,'(24x,F10.4)') cffreq(i)
       if (abs(cffreq(i)-freq(i)/invcm2ev).gt.1.0d0) then
          errmsg='Mismatch between calculated frequencies:'
          write(errmsg(len_trim(errmsg):),'(2(x,F10.4))') &
               cffreq(i),freq(i)/invcm2ev
          call error_control
       endif
    enddo

    close(unit)

!----------------------------------------------------------------------
! Read the symmetry labels from the CFOUR log file
!----------------------------------------------------------------------
    call freeunit(unit)
    open(unit,file=trim(freqfile)//'/cfour.log')

10  read(unit,'(a)',end=999) string
    if (index(string,'Normal Coordinate Analysis').eq.0) goto 10

    do i=1,12
       read(unit,*)
    enddo

    do i=1,nmodes
       read(unit,'(8x,a3)') nmlab(i)
    enddo

    close(unit)

    return

888 continue
    errmsg='The frequencies could not be found in: '&
         //trim(freqfile)//'/cfour.log'
    call error_control

999 continue
    errmsg='The Normal Coordinate Analysis section couldn''t &
         be found in: '//trim(freqfile)//'/cfour.log'
    call error_control

  end subroutine getmodes_cfour

!######################################################################

  subroutine getmodes_hessian

    use constants
    use iomod
    use sysinfo

    implicit none
    
    integer                        :: unit,i,j
    real(dp), dimension(ncoo,ncoo) :: hess

!----------------------------------------------------------------------
! Read the Hessian
!----------------------------------------------------------------------
    ! Open the Hessian file
    call freeunit(unit)
    open(unit,file=freqfile,form='formatted',status='old')

    ! Read the Hessian
    do i=1,natm+1
       read(unit,*)
    enddo
    do i=1,ncoo
       read(unit,*) (hess(i,j),j=1,ncoo)
    enddo

    ! Close the Hessian file
    close(unit)

!----------------------------------------------------------------------
! Calculate the normal mode vectors from the Hessian
!----------------------------------------------------------------------
    call hess2nm(hess)

!----------------------------------------------------------------------
! Assign the symmetry labels: only C1 symmetry is supported for now
!----------------------------------------------------------------------
    nmlab(1:nmodes)='A'
    
    return

  end subroutine getmodes_hessian

!######################################################################

  subroutine getmodes_aoforce

    use constants
    use sysinfo
    use iomod

    implicit none

    integer                        :: unit
    integer                        :: i,j,nm,nblock,indx1,indx2
    real(dp), dimension(ncoo,ncoo) :: matrix
    real(dp), dimension(6)         :: ftmp
    character(len=3), dimension(6) :: atmp
    character(len=120)             :: string

!----------------------------------------------------------------------
! Open the frequency file
!----------------------------------------------------------------------
    call freeunit(unit)
    open(unit,file=freqfile,form='formatted',status='old')

!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    matrix=0.0d0

!----------------------------------------------------------------------
! Read in the normal modes
!----------------------------------------------------------------------
    ! Read to the start of the non-zero frequency modes
5   read(unit,'(a)',end=100) string
    if (index(string,'reduced mass(g/mol)').eq.0) goto 5

    ! No. blocks of normal modes in the output
    nblock=ceiling(real(nmodes)/6)

    ! Read blocks of normal modes
    do i=1,nblock

       ! Indices of the first and last mode in the block
       indx1=(i-1)*6+1
       indx2=min(i*6,nmodes)

       ! Frequencies
       do j=1,4
          read(unit,*)
       enddo
       read(unit,'(20x,6(F9.2))') freq(indx1:indx2)

       ! Symmetry labels
       read(unit,*)
       read(unit,'(19x,6(6x,a3))') (nmlab(j),j=indx1,indx2)

       ! Normal mode vectors
       do j=1,8
          read(unit,*)
       enddo
       do j=1,ncoo
          read(unit,'(20x,6(F9.5))') matrix(j,indx1:indx2)
       enddo

       ! Skip to the end of this block
       read(unit,*)
       read(unit,*)
       
    enddo

!----------------------------------------------------------------------
! Close the frequency file
!----------------------------------------------------------------------
    close(unit)

!-----------------------------------------------------------------------
! Scale to mass-weighted x to obtain true normal mode vectors
!
! N.B. This has to be done as aoforce prints the normal mode vectors
!      in terms of non-mass-weighted Cartesians
!-----------------------------------------------------------------------
    do nm=1,nmodes
       do j=1,ncoo
          matrix(j,nm)=matrix(j,nm)*sqrt(mass(j))
       enddo
    enddo

!-----------------------------------------------------------------------
! Convert frequencies to eV
!-----------------------------------------------------------------------
    do i=1,nmodes
       freq(i)=freq(i)*invcm2ev
    enddo

!-----------------------------------------------------------------------
! Save the normal modes in the nmcoo array
!-----------------------------------------------------------------------
    nmcoo=matrix(1:ncoo,1:nmodes)    
    
    return

100 continue
    errmsg='The normal modes could not be found in the file: '&
         //trim(freqfile)
    
  end subroutine getmodes_aoforce
  
!######################################################################

  subroutine hess2nm(hess)

    use constants
    use iomod
    use sysinfo
    use utils

    implicit none

    integer                        :: i,j,k,l,e2,error,unit
    integer, dimension(ncoo)       :: indx
    real(dp), dimension(ncoo,ncoo) :: hess,proj,phess,eigvec
    real(dp), dimension(ncoo)      :: eigval
    real(dp), dimension(3*ncoo)    :: work

!----------------------------------------------------------------------
! Mass-weight the Hessian
!----------------------------------------------------------------------
    do i=1,ncoo
       do j=1,ncoo
          hess(i,j)=hess(i,j)/sqrt(mass(i)*mass(j))
       enddo
    enddo

!----------------------------------------------------------------------
! Project out the rotational and translational DOFs
!----------------------------------------------------------------------
    ! Construct the projector onto the complement of the
    ! translation-rotation subspace
    call rtproj(proj)

    ! Project out the rotational and translational DOFs
    phess=matmul(proj,matmul(hess,proj))

!----------------------------------------------------------------------
! Diagonalise the projected mass-weighted Hessian
!----------------------------------------------------------------------
    eigvec=phess
    e2=3*ncoo

    call dsyev('V','U',ncoo,eigvec,ncoo,eigval,work,e2,error)

    if (error.ne.0) then
       errmsg='Diagonalisation of the projected Hessian failed in &
            subroutine hess2nm'
       call error_control
    endif

!----------------------------------------------------------------------
! Sort the normal modes by absolute eigenvalues
!----------------------------------------------------------------------
    eigval=abs(eigval)
    call dsortindxa1('A',ncoo,eigval,indx)

!----------------------------------------------------------------------
! Frequencies in ev
!----------------------------------------------------------------------
    do i=7,ncoo
       freq(i-6)=sqrt(abs(eigval(indx(i))))*5140.8096d0*invcm2ev
    enddo

!-----------------------------------------------------------------------
! Normal mode vectors
!-----------------------------------------------------------------------
    do i=7,ncoo
       nmcoo(:,i-6)=eigvec(:,indx(i))
    enddo

    return

  end subroutine hess2nm

!######################################################################

  subroutine rtproj(proj)

    use constants
    use iomod
    use sysinfo

    implicit none

    integer                        :: i,j,k,l
    integer, dimension(6)          :: ipiv
    integer                        :: info
    real(dp), dimension(ncoo,ncoo) :: proj,rmat
    real(dp), dimension(6,ncoo)    :: rtvec
    real(dp), dimension(6,6)       :: smat,invsmat
    real(dp), dimension(ncoo,6)    :: bmat
    real(dp), dimension(6)         :: work
    logical(kind=4)                :: lcheck

!------------------------------------------------------------------
! Initialise arrays
!------------------------------------------------------------------
    rtvec=0.0d0

!------------------------------------------------------------------
! Vectors 1-3: translation along the three Cartesian axes
!------------------------------------------------------------------
    ! Loop over the translational DOFs
    do i=1,3
       ! Construct the vector for the current DOF
       do j=1,natm
          k=j*3-3+i
          rtvec(i,k)=sqrt(mass(j*3))
       enddo
    enddo

!------------------------------------------------------------------
! Vectors 4-6: infinitesimal displacements corresponding to
!              rotation about the three Cartesian axes
!------------------------------------------------------------------
    ! Rotation about the x-axis
    do i=1,natm
       j=i*3-1
       k=i*3
       rtvec(4,j)=sqrt(mass(i*3))*xcoo0(k)
       rtvec(4,k)=-sqrt(mass(i*3))*xcoo0(j)
    enddo

    ! Rotation about the y-axis
    do i=1,natm
       j=i*3-2
       k=i*3
       rtvec(5,j)=-sqrt(mass(i*3))*xcoo0(k)
       rtvec(5,k)=sqrt(mass(i*3))*xcoo0(j)
    enddo

    ! Rotation about the z-axis
    do i=1,natm
       j=i*3-2
       k=i*3-1
       rtvec(6,j)=sqrt(mass(i*3))*xcoo0(k)
       rtvec(6,k)=-sqrt(mass(i*3))*xcoo0(j)
    enddo

!------------------------------------------------------------------
! Calculate the projector R onto the translational and rotational
! DOFs using R=b*S^-1*b^T, where S=vec^T*vec.
!
! Here, R <-> rmat, S <-> smat, b <-> bmat (matrix of vectors)
!------------------------------------------------------------------
    ! Construct the b-matrix
    bmat=transpose(rtvec)

    ! Calculate the S-matrix
    smat=matmul(transpose(bmat),bmat)

    ! Invert the S-matrix
    invsmat=smat
    call dgetrf(6,6,invsmat,6,ipiv,info)
    if (info.ne.0) then
       errmsg='LU factorisation of the S-matrix failed'
       call error_control
    endif
    call dgetri(6,invsmat,6,ipiv,work,6,info)
    if (info.ne.0) then
       errmsg='Diagonalisation of the S-matrix failed'
       call error_control       
    endif
    
    ! Calculate the projection matrix R <-> rmat
    rmat=matmul(bmat,matmul(invsmat,transpose(bmat)))

!------------------------------------------------------------------
! Construct the projector onto the complement of the
! translation-rotation subspace: P=1-R
!------------------------------------------------------------------
    ! 1
    proj=0.0d0
    do i=1,ncoo
       proj(i,i)=1.0d0
    enddo

    ! 1-R
    proj=proj-rmat

    return

  end subroutine rtproj

!######################################################################

  subroutine nmortho

    use constants
    use sysinfo
    use channels

    implicit none

    integer  :: i,i1,j
    real(dp) :: x
      
    do i=1,nmodes
       
       do i1=1,i-1
          x=0.0d0
          do j=1,ncoo
             x=x+nmcoo(j,i)*nmcoo(j,i1)
          enddo
          do j=1,ncoo
             nmcoo(j,i)=nmcoo(j,i)-x*nmcoo(j,i1)
          enddo
       enddo
    
       x=0.0d0
       do j=1,ncoo
          x=x+nmcoo(j,i)*nmcoo(j,i)
       enddo
       
       x=1.0d0/sqrt(x)
       do j=1,ncoo
          nmcoo(j,i)=x*nmcoo(j,i)
       enddo

           enddo

    return

  end subroutine nmortho

!#######################################################################

  subroutine wrmodes

    use constants
    use sysinfo
    use iomod

    implicit none

    integer          :: unit,i,j,k
    character(len=2) :: amode

!----------------------------------------------------------------------
! Open the output file
!----------------------------------------------------------------------
    call freeunit(unit)
    open(unit,file='modes.xyz',form='formatted',status='unknown')

!----------------------------------------------------------------------
! Write the normal modes and frequencies to file
!----------------------------------------------------------------------
    do i=1,nmodes
       write(amode,'(i2)') i
       write(unit,'(i2)') natm
       write(unit,'(a,x,F10.4,x,a)') 'Q'//trim(adjustl(amode))&
            //', '//trim(adjustl(nmlab(i)))//',',freq(i)/invcm2ev,&
            'cm-1'
       do j=1,natm
          write(unit,'(a2,6(2x,F10.7))') atlbl(j),&
               (xcoo0(k)/ang2bohr,k=j*3-2,j*3),&
               (nmcoo(k,i)/sqrt(mass(k)),k=j*3-2,j*3)
       enddo
    enddo

!----------------------------------------------------------------------
! Close the output file
!----------------------------------------------------------------------
    close(unit)

    return

  end subroutine wrmodes
  
!######################################################################
! nmcoo transforms from nmodes to coo  x = nmcoo*Q
! coonm transforms from coo to nmodes  Q = coonm*(x-x0)
!
! freq in eV, mass in amu, x in non-mass-weighted Angstrom
!#######################################################################

  subroutine nm2xmat

    use constants
    use sysinfo

    implicit none

    integer :: i,j

!-----------------------------------------------------------------------
! Transposition of nmcoo
!-----------------------------------------------------------------------
    do j=1,nmodes         
       coonm(j,:)=nmcoo(:,j)
    enddo
      
!-----------------------------------------------------------------------
! Multiply matrix by mass and frequency factors
! "frequencies" in eV and masses in amu
!-----------------------------------------------------------------------
    do j=1,ncoo
       do i=1,nmodes
          nmcoo(j,i)=nmcoo(j,i)/(15.4644d0*sqrt(freq(i)*mass(j)))
          coonm(i,j)=coonm(i,j)*(15.4644d0*sqrt(freq(i)*mass(j)))
       enddo
    enddo
    
    return
  
  end subroutine nm2xmat

!######################################################################

  subroutine geten(v,nsta,filename)

    use constants
    use iomod
    
    implicit none

    integer, intent(in)          :: nsta
    integer                      :: entyp
    real(dp), dimension(nsta)    :: v
    character(len=*), intent(in) :: filename

!----------------------------------------------------------------------
! Determine the quantum chemistry program used to calculate the
! adiabatic energies
!----------------------------------------------------------------------
    call entype(entyp,filename)

!----------------------------------------------------------------------
! Read the energies
!----------------------------------------------------------------------
    if (entyp.eq.5) then
       ! DFT/MRCI
       call geten_dftmrci(v,nsta,filename)
    else if (entyp.eq.6) then
       ! Columbus MRCI
       call geten_colmrci(v,nsta,filename)
    endif
    
    return
    
  end subroutine geten

!######################################################################

  subroutine geten_dftmrci(v,nsta,filename)

    use constants
    use iomod

    implicit none

    integer, intent(in)          :: nsta
    integer                      :: unit,ista
    real(dp), dimension(nsta)    :: v
    character(len=*), intent(in) :: filename
    character(len=120)           :: string,searchstring
    
!----------------------------------------------------------------------
! Open the DFT/MRCI output file
!----------------------------------------------------------------------
    call freeunit(unit)
    open(unit,file=filename,form='formatted',status='old')

!----------------------------------------------------------------------
! Read the energies
!----------------------------------------------------------------------
    v=0.0d0
    ista=0

    searchstring='state         method'
    
5   read(unit,'(a)',end=100) string
    if (index(string,trim(searchstring)).ne.0) then
       ista=ista+1
       read(unit,*)
       read(unit,'(30x,F13.6)') v(ista)
    endif

    if (ista.eq.nsta) then
       goto 100
    else
       goto 5
    endif
       
100 continue

!----------------------------------------------------------------------
! Exit if not all state energies were found
!----------------------------------------------------------------------
    if (ista.ne.nsta) then
       errmsg='Not all state energies found in: '//trim(filename)
       call error_control
    endif
       
!----------------------------------------------------------------------
! Close the DFT/MRCI output file
!----------------------------------------------------------------------
    close(unit)

    return
    
  end subroutine geten_dftmrci

!######################################################################

  subroutine geten_colmrci(v,nsta,filename)

    use constants
    use iomod

    implicit none

    integer, intent(in)          :: nsta
    integer                      :: unit,ista
    real(dp), dimension(nsta)    :: v
    character(len=*), intent(in) :: filename
    character(len=120)           :: string,searchstring
    
!----------------------------------------------------------------------
! Open the Columbus MRCI output file
!----------------------------------------------------------------------
    call freeunit(unit)
    open(unit,file=filename,form='formatted',status='old')

!----------------------------------------------------------------------
! Read the energies
!----------------------------------------------------------------------
    v=0.0d0
    ista=0

    searchstring='eci       ='
    
5   read(unit,'(a)',end=100) string
    if (index(string,trim(searchstring)).ne.0) then
       ista=ista+1
       read(string,'(12x,F20.12)') v(ista)
    endif

    if (ista.eq.nsta) then
       goto 100
    else
       goto 5
    endif
       
100 continue

!----------------------------------------------------------------------
! Exit if not all state energies were found
!----------------------------------------------------------------------
    if (ista.ne.nsta) then
       errmsg='Not all state energies found in: '//trim(filename)
       call error_control
    endif

!----------------------------------------------------------------------
! Close the Columbus MRCI output file
!----------------------------------------------------------------------
    close(unit)

    return
    
  end subroutine geten_colmrci
  
!######################################################################
! entype: determines the quantum chemistry program used for the
!         energy calculation, and sets entyp accordingly:
!          
!         entyp = 5 <-> DFT/MRCI
!                 6 <-> Columbus MRCI
!######################################################################

  subroutine entype(entyp,filename)

    use constants
    use iomod
    
    implicit none

    integer                      :: entyp
    character(len=*), intent(in) :: filename

!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    entyp=-1

!----------------------------------------------------------------------
! Determine the quantum chemistry program used for the energy
! calculation
!---------------------------------------------------------------------
    if (isdftmrci(filename)) then
       ! DFT/MRCI
       entyp=5
    else if (iscolmrci(filename)) then
       ! Columbus MRCI
       entyp=6
    endif
    
!----------------------------------------------------------------------
! Check that the quantum chemistry program used is supported
!----------------------------------------------------------------------
    if (entyp.eq.-1) then
       errmsg='The quantum chemistry program used for the &
            energy calculation is not supported.'
       call error_control
    endif
    
    return
    
  end subroutine entype

!######################################################################

  function isdftmrci(filename) result(found)

    use constants
    use iomod

    implicit none

    integer            :: unit
    character(len=*)   :: filename
    character(len=120) :: string
    logical            :: found

!----------------------------------------------------------------------
! Open file
!----------------------------------------------------------------------
    call freeunit(unit)
    open(unit,file=filename,form='formatted',status='old')

!----------------------------------------------------------------------
! Check whether the calculation was performed using DFT/MRCI code
!----------------------------------------------------------------------
    found=.false.
    
5   read(unit,'(a)',end=10) string
    if (index(string,' *              M R - C I               *').ne.0) then
       found=.true.
    else
       goto 5
    endif
10  continue

!----------------------------------------------------------------------
! Close file
!----------------------------------------------------------------------
    close(unit)
    
    return
    
  end function isdftmrci

!######################################################################
  
 function iscolmrci(filename) result(found)

    use constants
    use iomod

    implicit none

    integer            :: unit
    character(len=*)   :: filename
    character(len=120) :: string
    logical            :: found

!----------------------------------------------------------------------
! Open file
!----------------------------------------------------------------------
    call freeunit(unit)
    open(unit,file=filename,form='formatted',status='old')

!----------------------------------------------------------------------
! Check whether the calculation was a Columbus MRCI calculation
!----------------------------------------------------------------------
    found=.false.
    
5   read(unit,'(a)',end=10) string
    if (index(string,'program ciudg').ne.0) then
       found=.true.
    else
       goto 5
    endif
10  continue
    
!----------------------------------------------------------------------
! Close file
!----------------------------------------------------------------------
    close(unit)
    
    return
    
  end function iscolmrci

!######################################################################
  
end module ioqc

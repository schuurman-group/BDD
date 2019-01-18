!######################################################################
! displace: A program for the creation of cuts along mass- and 
!           frequency-weighted normal modes.
!######################################################################

program displace

  use ioqc
  
  implicit none

!----------------------------------------------------------------------
! Open the input and log files
!----------------------------------------------------------------------
  call open_files

!----------------------------------------------------------------------
! Read the input file
!----------------------------------------------------------------------
  call rdinpfile

!----------------------------------------------------------------------
! Determine the normal mode file type
!----------------------------------------------------------------------
  call freqtype

!----------------------------------------------------------------------
! Determine the no. atoms and allocate xcoo0 and related arrays
!----------------------------------------------------------------------
  call getdim

!----------------------------------------------------------------------
! Read the Cartesian coordinates
!----------------------------------------------------------------------
  call getxcoo0
  
!----------------------------------------------------------------------
! Determine the number of normal modes from the moment of intertia
! tensor and allocate associated arrays
!----------------------------------------------------------------------
  call getnmodes

!----------------------------------------------------------------------
! Read the normal modes, frequencies and symmetry labels
!----------------------------------------------------------------------
  call getmodes

!----------------------------------------------------------------------
! Create the transformation matrices
!----------------------------------------------------------------------
  call nm2xmat

!----------------------------------------------------------------------
! Make the cuts
!----------------------------------------------------------------------
  call makecut
  
contains

!######################################################################

  subroutine open_files

    use constants
    use channels
    use iomod

    implicit none

    integer :: ilbl
    logical :: found
    
!----------------------------------------------------------------------
! Exit if no input file has been given
!----------------------------------------------------------------------
    if (iargc().eq.0) then
       write(6,'(/,2x,a,/)') 'No input file has been given'
       stop
    endif
    
!----------------------------------------------------------------------
! Read the name of the input file and set the name of the log file
!----------------------------------------------------------------------
    call getarg(1,ain)

    ilbl=index(ain,'.inp')
    if (ilbl.eq.0) then
       alog=trim(ain)//'.log'
       ain=trim(ain)//'.inp'
    else
       alog=ain(1:ilbl-1)//'.log'
    endif

!----------------------------------------------------------------------
! Exit if the input file does not exist
!----------------------------------------------------------------------
    inquire(file=trim(ain),exist=found)
    
    if (.not.found) then
       write(6,'(/,2x,a,/)') 'The file '//trim(ain)//' does not exist'
       stop
    endif

!----------------------------------------------------------------------
! Open the input and log files
!----------------------------------------------------------------------
    iin=1
    open(iin,file=ain,form='formatted',status='old')
    
    ilog=2
    open(ilog,file=alog,form='formatted',status='unknown')

    return
    
  end subroutine open_files

!######################################################################

  subroutine rdinpfile

    use constants
    use channels
    use iomod
    use parsemod
    use sysinfo
    use symmetry
    use dispglobal

    implicit none

    integer :: i,k,ilbl
    
!----------------------------------------------------------------------
! Set defaults
!----------------------------------------------------------------------
    freqfile=''
    icut=0
    npnts=0
    dq=0.0d0
    pntgrp=''
    nsta=0
    
!----------------------------------------------------------------------
! Read the input file
!----------------------------------------------------------------------
    rewind(iin)

15  continue
    call rdinp(iin)

    i=0
    if (.not.lend) then
    
20     continue
       i=i+1

       if (keyword(i).eq.'$freqfile') then
          if (keyword(i+1).eq.'=') then
             i=i+2
             freqfile=keyword(i)
          else
             goto 100
          endif

       else if (keyword(i).eq.'$cut') then
          if (keyword(i+1).eq.'=') then
             i=i+2
             if (keyword(i).eq.'all_1d') then
                ! One-dimensional cuts only
                icut=1
             else if (keyword(i).eq.'all_2d') then
                ! All cuts (1D and 2D)
                icut=2
             else
                goto 100
          endif
          ! Stepsize
          if (keyword(i+1).eq.',') then
             i=i+2
             read(keyword(i),*) dq
          else
             errmsg='The stepsize has not given with the &
                  all_1d keyword'
             call error_control
          endif
          ! Number of points in each direction
          if (keyword(i+1).eq.',') then
             i=i+2
             read(keyword(i),*) npnts
          else
             errmsg='The no. points has not given with the &
                  all_1d keyword'
             call error_control
          endif
       else
          errmsg='Unkown cut type: '//trim(keyword(i))
          call error_control
       endif
       
       else if (keyword(i).eq.'$point_group') then
          if (keyword(i+1).eq.'=') then
             i=i+2
             pntgrp=keyword(i)
          else
             goto 100
          endif

       else if (keyword(i).eq.'$state_sym') then
          ! First pass: determine the no. states and allocate arrays
          do
             call rdinp(iin)
             if (lend) then
                errmsg='End of file reached whilst reading the &
                     $state_sym section'
                call error_control
             endif
             if (keyword(1).eq.'$end') exit
             nsta=nsta+1
          enddo
          allocate(stalab(nsta))
          ! Second pass: read in the state symmetry labels
          do k=1,nsta+1
             backspace(iin)
          enddo
          do k=1,nsta
             call rdinp(iin)
             read(keyword(2),*) ilbl
             stalab(ilbl)=keyword(1)
          enddo
          ! Skip past the $end keyword
          call rdinp(iin)
          i=inkw
          
       else
          ! Exit if the keyword is not recognised
          errmsg='Unknown keyword: '//trim(keyword(i))
          call error_control
       endif
       
       ! If there are more keywords to be read on the current line,
       ! then read them, else read the next line
       if (i.lt.inkw) then
          goto 20
       else
          goto 15
       endif
       
       ! Exit if a required argument has not been given with a keyword
100    continue
       errmsg='No argument given with the keyword '//trim(keyword(i))
       call error_control
       
    endif

!----------------------------------------------------------------------
! Make sure that all the required information has been given
!----------------------------------------------------------------------
    ! Name of the frequency file
    if (freqfile.eq.'') then
       errmsg='The name of the frequency calculation file has not &
            been given'
       call error_control
    endif

    ! Name of the point group
    if (icut.eq.2.and.pntgrp.eq.'') then
       errmsg='For full cuts ($cut=all_2d), the point group must &
            be specified'
       call error_control
    endif
    
    return
    
  end subroutine rdinpfile
    
!######################################################################

  subroutine makecut

    use constants
    use sysinfo
    use symmetry
    use dispglobal
    
    implicit none

    integer :: n
    
!----------------------------------------------------------------------
! Reference geometry
!----------------------------------------------------------------------
    call write_1file(xcoo0/ang2bohr,'q0.xyz')
    
!----------------------------------------------------------------------
! 1D cuts: all included
!----------------------------------------------------------------------
    do n=1,nmodes
       call makecut_1d(n)
    enddo

!----------------------------------------------------------------------
! 2D cuts: only pairs of modes corresponding to symmetry-allowed
! coupling coefficients are included
!----------------------------------------------------------------------
    if (icut.eq.2) then

       ! Determine which coupling coefficients are zero by symmetry
       call create_mask

       ! Determine which 2D cuts to make
       call get_cutmask

       ! Make the cuts
       print*,'The 2D cut code still needs writing...'
       STOP
       
    endif
       
    return
    
  end subroutine makecut

!######################################################################

  subroutine makecut_1d(n)

    use constants
    use channels
    use iomod
    use sysinfo
    use dispglobal
    
    implicit none

    integer, intent(in)         :: n
    integer                     :: i
    real(dp), dimension(nmodes) :: q
    real(dp), dimension(ncoo)   :: x
    character(len=60)           :: filename
    character(len=3)            :: aq,ai
    character(len=1)            :: ad
    
    write(aq,'(i3)') n

    ! Loop over displacements
    do i=-npnts,npnts

       ! Skip the reference point
       if (i.eq.0) cycle
       
       ! Point in normal modes
       q=0.0d0
       q(n)=i*dq

       ! Cartesian coordinates
       x=xcoo0/ang2bohr+matmul(nmcoo,q)

       ! Filename
       write(ai,'(i3)') abs(i)
       if (i.lt.0) then
          ad='l'
       else if (i.gt.0) then
          ad='r'
       endif
       filename='q'//trim(adjustl(aq))//'_'//trim(adjustl(ai)) &
            //ad//'.xyz'
       
       ! Write the Cartesian coordinates to file
       call write_1file(x,filename)       
       
    enddo
       
    return
    
  end subroutine makecut_1d

!######################################################################

  subroutine write_1file(x,filename)

    use constants
    use iomod
    use sysinfo
    
    implicit none

    integer                   :: i,j,unit
    real(dp), dimension(ncoo) :: x
    character(len=*)          :: filename
    character(len=3)          :: an
    
!----------------------------------------------------------------------
! Open the xyz file
!----------------------------------------------------------------------
    call freeunit(unit)
    open(unit,file=filename,form='formatted',status='unknown')

!----------------------------------------------------------------------
! Write the xyz file
!----------------------------------------------------------------------
    write(an,'(i3)') natm
    write(unit,'(a)') trim(adjustl(an))
    write(unit,*)
    do i=1,natm
       write(unit,'(a2,3(2x,F10.7))') atlbl(i),(x(j),j=i*3-2,i*3)
    enddo
    
!----------------------------------------------------------------------
! Close the xyz file
!----------------------------------------------------------------------
    close(unit)
    
    return
    
  end subroutine write_1file
    
!######################################################################

  subroutine get_cutmask

    use constants
    use sysinfo
    use symmetry
    use dispglobal
    
    implicit none

    integer :: m1,m2,s,s1,s2
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(cut_mask(nmodes,nmodes))
    cut_mask=0

!----------------------------------------------------------------------
! Determine which 2D cuts to make based on the symmetry-allowed
! coupling coefficients
!----------------------------------------------------------------------
    ! gamma
    do s=1,nsta
       do m1=1,nmodes-1
          do m2=m1+1,nmodes
             if (gamma_mask(m1,m2,s).eq.1) then
                cut_mask(m1,m2)=1
                cut_mask(m2,m1)=1
             endif
          enddo
       enddo
    enddo
    
    ! mu
    do s1=1,nsta-1
       do s2=s1+1,nsta
          do m1=1,nmodes-1
             do m2=m1+1,nmodes
                if (mu_mask(m1,m2,s1,s2).eq.1) then
                   cut_mask(m1,m2)=1
                   cut_mask(m2,m1)=1
                endif
             enddo
          enddo
       enddo
    enddo
    
    return
    
  end subroutine get_cutmask
  
!######################################################################
  
end program displace

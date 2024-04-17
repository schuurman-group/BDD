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
  call open_files_displace

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

  subroutine open_files_displace

    use constants
    use channels
    use iomod

    implicit none

    integer :: ilbl
    logical :: found
    
!----------------------------------------------------------------------
! Exit if no input file has been given
!----------------------------------------------------------------------
    if (iargc() == 0) then
       write(6,'(/,2x,a,/)') 'No input file has been given'
       stop
    endif
    
!----------------------------------------------------------------------
! Read the name of the input file and set the name of the log file
!----------------------------------------------------------------------
    call getarg(1,ain)

    ilbl=index(ain,'.inp')
    if (ilbl == 0) then
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

!----------------------------------------------------------------------
! Create the geometries directory
!----------------------------------------------------------------------
    ! Works with ifort
    inquire(directory='geoms',exist=found)

    ! Works with gfortran
    !inquire(file='geoms/.',exist=found)

    if (found) then
       call system('rm -rf geoms/*')
    else
       call system('mkdir geoms')
    endif

    return
    
  end subroutine open_files_displace

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

       if (keyword(i) == '$freqfile') then
          if (keyword(i+1) == '=') then
             i=i+2
             freqfile=keyword(i)
          else
             goto 100
          endif

       else if (keyword(i) == '$cut') then
          if (keyword(i+1) == '=') then
             i=i+2
             if (keyword(i) == '1mode' &
                  .or. keyword(i) == 'all_1d') then
                ! One-mode cuts only
                icut=1
             else if (keyword(i) == '2mode') then
                ! Diagonal 2-mode cuts needed for the determination
                ! of all 2-mode terms
                icut=2
             else if (keyword(i) == '2mode_ondiag' &
                  .or. keyword(i) == 'gamma_diag_2d') then
                ! Diagonal 2-mode cuts needed for the determination
                ! of the on-diagonal two-mode terms
                icut=3
             else
                goto 100
          endif
          ! Stepsize
          if (keyword(i+1) == ',') then
             i=i+2
             read(keyword(i),*) dq
          else
             errmsg='The stepsize has not given with the &
                  all_1d keyword'
             call error_control
          endif
          ! Number of points in each direction
          if (keyword(i+1) == ',') then
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
       
       else if (keyword(i) == '$point_group') then
          if (keyword(i+1) == '=') then
             i=i+2
             pntgrp=keyword(i)
          else
             goto 100
          endif

       else if (keyword(i) == '$state_sym') then
          ! First pass: determine the no. states and allocate arrays
          do
             call rdinp(iin)
             if (lend) then
                errmsg='End of file reached whilst reading the &
                     $state_sym section'
                call error_control
             endif
             if (keyword(1) == '$end') exit
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
       if (i < inkw) then
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
    if (freqfile == '') then
       errmsg='The name of the frequency calculation file has not &
            been given'
       call error_control
    endif

    ! Name of the point group
    if (icut > 1 .and. pntgrp == '') then
       errmsg='For 2D cuts, the point group must be specified'
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

    integer :: n,n1,n2
    
!----------------------------------------------------------------------
! Reference geometry
!----------------------------------------------------------------------
    call write_1file(xcoo0/ang2bohr,'q0.xyz')
    
!----------------------------------------------------------------------
! 1-mode cuts
!----------------------------------------------------------------------
    if (icut == 1) then
       do n=1,nmodes
          call makecut_1d(n)
       enddo
    endif

!----------------------------------------------------------------------
! Diagonal 2-mode cuts needed for the determination of all 2-mode
! coupling coefficients
!----------------------------------------------------------------------
    if (icut == 2) then

       ! Determine which coupling coefficients are zero by symmetry
       call create_mask(2)

       ! Determine which 2D cuts to make
       call get_cutmask
       
       do n1=1,nmodes-1
          do n2=n1+1,nmodes
             if (cut_mask(n1,n2) == 1) call makecut_2d_diag(n1,n2)
          enddo
       enddo
       
    endif

    
!----------------------------------------------------------------------
! Diagonal 2-mode cuts needed for the determination of only the
! on-diagonal 2-mode coupling coefficients
!----------------------------------------------------------------------
    if (icut == 3) then

       ! Determine which coupling coefficients are zero by symmetry
       call create_mask(2)

       do n1=1,nmodes-1
          do n2=n1+1,nmodes
             !if (cut_mask(n1,n2) == 1) call makecut_2d_diag(n1,n2)

             if (coeff2_mask(n1,n2,1,1) == 1) &
                  call makecut_2d_diag(n1,n2)
             
          enddo
       enddo
       
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
    integer                     :: i,j,k,unit
    real(dp), dimension(nmodes) :: q
    real(dp), dimension(ncoo)   :: x
    character(len=60)           :: filename
    character(len=3)            :: aq,ai,an
    character(len=1)            :: ad

    write(aq,'(i3)') n

!----------------------------------------------------------------------
! Concatenated xyz files
!----------------------------------------------------------------------
    !
    ! Positive direction
    !
    ! Open the concatenated xyz file
    call freeunit(unit)
    filename='q'//trim(adjustl(aq))//'r.xyz'
    open(unit,file='geoms/'//trim(filename),form='formatted',&
         status='unknown')
    
    ! Loop over displacements
    do i=0,npnts

       ! Point in normal modes
       q=0.0d0
       q(n)=i*dq

       ! Cartesian coordinates
       x=xcoo0/ang2bohr+matmul(nmcoo,q)

       ! Write the Cartesian coordintates to file
       write(an,'(i3)') natm
       write(unit,'(a)') trim(adjustl(an))
       write(unit,*)
       do j=1,natm
          write(unit,'(a2,3(2x,F12.7))') atlbl(j),(x(k),k=j*3-2,j*3)
       enddo
       
    enddo
    
    ! Close the concatenated xyz file
    close(unit)

    !
    ! Negative direction
    !
    ! Open the concatenated xyz file
    call freeunit(unit)
    filename='q'//trim(adjustl(aq))//'l.xyz'
    open(unit,file='geoms/'//trim(filename),form='formatted',&
         status='unknown')

    ! Loop over displacements
    do i=0,-npnts,-1

       ! Point in normal modes
       q=0.0d0
       q(n)=i*dq
       
       ! Cartesian coordinates
       x=xcoo0/ang2bohr+matmul(nmcoo,q)
       
       ! Write the Cartesian coordintates to file
       write(an,'(i3)') natm
       write(unit,'(a)') trim(adjustl(an))
       write(unit,*)
       do j=1,natm
          write(unit,'(a2,3(2x,F12.7))') atlbl(j),(x(k),k=j*3-2,j*3)
       enddo
       
    enddo
    
    ! Close the concatenated xyz file
    close(unit)
    
    return
    
  end subroutine makecut_1d

!######################################################################

  subroutine makecut_2d_diag(n1,n2)

    use constants
    use channels
    use iomod
    use sysinfo
    use dispglobal
    
    implicit none

    integer, intent(in)         :: n1,n2
    integer                     :: i,j,k,unit
    real(dp), dimension(nmodes) :: q
    real(dp), dimension(ncoo)   :: x
    character(len=60)           :: filename
    character(len=3)            :: aq1,aq2,ai,an
    character(len=1)            :: ad1,ad2

    write(aq1,'(i3)') n1
    write(aq2,'(i3)') n2

!----------------------------------------------------------------------
! Concatenated xyz files
!----------------------------------------------------------------------
    !
    ! Positive direction
    !
    ! Open the concatenated xyz file
    call freeunit(unit)
    filename='q'//trim(adjustl(aq1))//'r_q'//trim(adjustl(aq2)) &
         //'r.xyz'
    open(unit,file='geoms/'//trim(filename),form='formatted',&
         status='unknown')

    ! Loop over displacements
    do i=0,npnts

       ! Point in normal modes
       q=0.0d0
       q(n1)=i*dq/sqrt(2.0d0)
       q(n2)=i*dq/sqrt(2.0d0)

       ! Cartesian coordinates
       x=xcoo0/ang2bohr+matmul(nmcoo,q)

       ! Write the Cartesian coordintates to file
       write(an,'(i3)') natm
       write(unit,'(a)') trim(adjustl(an))
       write(unit,*)
       do j=1,natm
          write(unit,'(a2,3(2x,F12.7))') atlbl(j),(x(k),k=j*3-2,j*3)
       enddo
       
    enddo
    
    ! Close the concatenated xyz file
    close(unit)

    !
    ! Negative direction
    !
    ! Open the concatenated xyz file
    call freeunit(unit)
    filename='q'//trim(adjustl(aq1))//'l_q'//trim(adjustl(aq2)) &
         //'l.xyz'
    open(unit,file='geoms/'//trim(filename),form='formatted',&
         status='unknown')

    ! Loop over displacements
    do i=0,-npnts,-1

       ! Point in normal modes
       q=0.0d0
       q(n1)=i*dq/sqrt(2.0d0)
       q(n2)=i*dq/sqrt(2.0d0)

       ! Cartesian coordinates
       x=xcoo0/ang2bohr+matmul(nmcoo,q)

       ! Write the Cartesian coordintates to file
       write(an,'(i3)') natm
       write(unit,'(a)') trim(adjustl(an))
       write(unit,*)
       do j=1,natm
          write(unit,'(a2,3(2x,F12.7))') atlbl(j),(x(k),k=j*3-2,j*3)
       enddo
       
    enddo
       
    ! Close the concatenated xyz file
    close(unit)
    
    return
    
  end subroutine makecut_2d_diag
    
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
    open(unit,file='geoms/'//trim(filename),form='formatted',&
         status='unknown')

!----------------------------------------------------------------------
! Write the xyz file
!----------------------------------------------------------------------
    write(an,'(i3)') natm
    write(unit,'(a)') trim(adjustl(an))
    write(unit,*)
    do i=1,natm
       write(unit,'(a2,3(2x,F12.7))') atlbl(i),(x(j),j=i*3-2,i*3)
    enddo
    
!----------------------------------------------------------------------
! Close the xyz file
!----------------------------------------------------------------------
    close(unit)
    
    return
    
  end subroutine write_1file
  
!######################################################################
  
end program displace

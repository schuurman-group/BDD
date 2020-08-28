!######################################################################
! gridfunc: A program for the calculation of various functions on
!           direct product grids:
!           (1) Projectors onto adiabatic states
!           (2) Rotational energies
!######################################################################
program gridfunc

  use constants
  use ioqc
  use dvr
  
  implicit none

!----------------------------------------------------------------------
! Open the input and log files
!----------------------------------------------------------------------
  call open_files_gridfunc

!----------------------------------------------------------------------
! Read the input file
!----------------------------------------------------------------------
  call rdinpfile_gridfunc

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
! Create the normal mode transformation matrices
!----------------------------------------------------------------------
  call nm2xmat

!----------------------------------------------------------------------
! Parse the primitive basis section in the input file
!----------------------------------------------------------------------
  call rdprimbasis

!----------------------------------------------------------------------
! Compute the 1D DVR grids
!----------------------------------------------------------------------
  call dvr_grids
  
contains

!######################################################################

  subroutine open_files_gridfunc

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
    
  end subroutine open_files_gridfunc
    
!######################################################################

  subroutine rdinpfile_gridfunc

    use constants
    use channels
    use iomod
    use parsemod
    use sysinfo
    use gridglobal
    
    implicit none

    integer :: i,i1,n
    logical :: primsection
    
!----------------------------------------------------------------------
! Set defaults
!----------------------------------------------------------------------
    ! Frequency file
    freqfile=''

    ! Function type
    ifunc=0

    ! Electronic state indices
    funcsta=0

    ! Number of modes entering into the function
    nfuncmode=0

    ! Potential parameter file
    abin=''

    ! Primitive section found
    primsection=.false.
    
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

       else if (keyword(i).eq.'$potfile') then
          if (keyword(i+1).eq.'=') then
             i=i+2
             abin=keyword(i)
          else
             goto 100
          endif

       else if (keyword(i).eq.'$function') then
          if (keyword(i+1).eq.'=') then
             i=i+2
             if (keyword(i).eq.'adiab_proj') then
                ifunc=1
             else
                errmsg='Unknown function type: '//trim(keyword(i))
                call error_control
             endif
          else
             goto 100
          endif

       else if (keyword(i).eq.'$states') then
          if (keyword(i+1).eq.'=') then
             i=i+2
             read(keyword(i),*) funcsta(1)
             if (keyword(i+1).eq.',') then
                i=i+2
                read(keyword(i),*) funcsta(2)
             endif
          else
             goto 100
          endif

       else if (keyword(i).eq.'$modes') then
          if (keyword(i+1).eq.'=') then
             i=i+2
             ! Determine the no. modes in the function
             i1=i
             nfuncmode=1
25           continue
             if (keyword(i1+1).eq.',') then
                nfuncmode=nfuncmode+1
                i1=i1+2
                goto 25
             endif
             ! Allocate the funcmode array
             allocate(funcmode(nfuncmode))
             ! Read the function modes
             n=1
             read(keyword(i),*) funcmode(n)
35           continue
             if (keyword(i+1).eq.',') then
                i=i+2
                n=n+1
                read(keyword(i),*) funcmode(n)
                goto 35
             endif
          else
             goto 100
          endif

       else if (index(keyword(i),'pbasis-section').ne.0) then
          ! Read past the primitive basis section: this will
          ! be parsed after we know the no. modes
          primsection=.true.
40        continue
          call rdinp(iin)
          if (index(keyword(1),'end-pbasis-section').eq.0) goto 40

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
    if (freqfile.eq.'') then
       errmsg='The name of the frequency calculation file has not &
            been given'
       call error_control
    endif

    if (abin.eq.'') then
       errmsg='The name of the potential parameter file has not been &
            given'
       call error_control
    endif

    if (ifunc.eq.0) then
        errmsg='The function type has not been given'
        call error_control
    endif

    if (funcsta(1).eq.0) then
       errmsg='The function states have not been given'
       call error_control
    endif

    if (nfuncmode.eq.0) then
       errmsg='The function modes have not been given'
       call error_control
    endif

    if (.not.primsection) then
       errmsg='The primitive basis section could not be found'
       call error_control
    endif
    
    return
    
  end subroutine rdinpfile_gridfunc
    
!######################################################################

  subroutine rdprimbasis

    use constants
    use channels
    use iomod
    use parsemod
    use sysinfo
    use gridglobal
    
    implicit none

    integer :: m,k,i

!----------------------------------------------------------------------
! Allocate and initialise arrays
!----------------------------------------------------------------------
    ! DVR types
    allocate(idvr(nmodes))
    idvr=0

    ! DVR parameters
    allocate(dvrpar(nmodes,maxdvrpar))
    dvrpar=0.0d0

    ! No. DVR basis functions
    allocate(ndvr(nmodes))
    ndvr=0
    
!----------------------------------------------------------------------
! Read to the primitive basis section
!----------------------------------------------------------------------
    rewind(iin)
5   call rdinp(iin)
    if (index(keyword(1),'pbasis-section').eq.0) goto 5

!----------------------------------------------------------------------
! Parse the primitive basis section
! Note that here we are assuming that the mode labels are the same as
! those written to the operator file by the KDC code, i.e., Q1, Q2,...
!----------------------------------------------------------------------
10  call rdinp(iin)
    if (index(keyword(1),'end-pbasis-section').eq.0) then

       ! Mode index (assuming a mode label Qm)
       read(keyword(1)(2:),*) m

       ! DVR type
       if (keyword(2).eq.'ho') then
          ! Harmonic oscillator DVR
          idvr(m)=1
       else
          errmsg='Unknown DVR type for mode '//trim(keyword(1))
          call error_control
       endif

       ! Number of DVR basis functions
       read(keyword(3),*) ndvr(m)

       ! DVR parameters (for now we will assume that all parameters
       ! are given in atomic units)
       k=0
       do i=4,inkw
          k=k+1
          read(keyword(i),*) dvrpar(m,k)
       enddo
       
       goto 10
    endif
    
    return
    
  end subroutine rdprimbasis
  
!######################################################################
  
end program gridfunc
  

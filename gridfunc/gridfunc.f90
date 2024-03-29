!######################################################################
! gridfunc: A program for the calculation of various functions on
!           direct product grids:
!           (1) Projectors onto adiabatic states
!           (2) Rotational energies
!######################################################################
program gridfunc

  use constants
  use gridglobal
  use ioqc
  use extfunc
  use dvr
  use func


  use sysinfo

  
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
  if (idiabfunc == 1) call freqtype

!----------------------------------------------------------------------
! Determine the no. atoms and allocate xcoo0 and related arrays
!----------------------------------------------------------------------
  if (idiabfunc == 1) then
     ! KDC potential or rotational energy functions
     call getdim
  else
     ! Model potential functions
     call getdim_extfunc(idiabfunc)
  endif
     
!----------------------------------------------------------------------
! Read the Cartesian coordinates
!----------------------------------------------------------------------
  if (idiabfunc == 1) call getxcoo0
  
!----------------------------------------------------------------------
! Determine the number of normal modes from the moment of intertia
! tensor and allocate associated arrays
!----------------------------------------------------------------------
  if (idiabfunc == 1) call getnmodes

!----------------------------------------------------------------------
! Read the normal modes, frequencies and symmetry labels
!----------------------------------------------------------------------
  if (idiabfunc == 1) call getmodes

!----------------------------------------------------------------------
! Create the normal mode transformation matrices
!----------------------------------------------------------------------
  if (idiabfunc == 1) call nm2xmat

!----------------------------------------------------------------------
! Parse the primitive basis section in the input file
!----------------------------------------------------------------------
  call rdprimbasis

!----------------------------------------------------------------------
! Read the potential file
!----------------------------------------------------------------------
  if (idiabfunc == 1) call rdpotfile
  
!----------------------------------------------------------------------
! Compute the 1D DVR grids
!----------------------------------------------------------------------
  call dvr_grids

!----------------------------------------------------------------------
! Compute and output the diabatic state representation of the
! requested function on the direct product (sub) grid
!----------------------------------------------------------------------
  call calc_func
  
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

    ! J
    Jval=0

    ! Diabatic function index
    ! Default: vibronic coupling Hamiltonian with parameters
    ! read from a KDC potential parameter file
    idiabfunc=1
    
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

       else if (keyword(i).eq.'$potfunc') then
          if (keyword(i+1).eq.'=') then
             i=i+2
             if (keyword(i).eq.'bma_lvc') then
                idiabfunc=2
             else if (keyword(i).eq.'c4h4_lvc') then
                idiabfunc=3
             else
                errmsg='Unknown potential function: '//trim(keyword(i))
                call error_control
             endif
          else
             goto 100
          endif

       else if (keyword(i).eq.'$function') then
          if (keyword(i+1).eq.'=') then
             i=i+2
             if (keyword(i).eq.'adiab_proj') then
                ! Projector onto an adiabatic state
                ifunc=1
             else if (keyword(i).eq.'adiab_exci') then
                ! Adiabatic state excitation operator
                ifunc=2
             else if (keyword(i).eq.'roten') then
                ! Rotational energy on the grid
                ifunc=3
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

       else if (keyword(i).eq.'$j') then
          if (keyword(i+1).eq.'=') then
             i=i+2
             read(keyword(i),*) Jval
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
    if (freqfile.eq.'' .and. idiabfunc.eq.1) then
       errmsg='The name of the frequency calculation file has not &
            been given'
       call error_control
    endif

    if (abin.eq.'' .and. idiabfunc.eq.1) then
       errmsg='The name of the potential parameter file has not been &
            given'
       call error_control
    endif

    if (ifunc.eq.0) then
        errmsg='The function type has not been given'
        call error_control
    endif

    if (ifunc.ne.3) then
       if (funcsta(1).eq.0) then
          errmsg='The function states have not been given'
          call error_control
       endif
    endif
    
    if (ifunc.eq.2.and.funcsta(2).eq.0) then
       errmsg='Adiabatic excitation operators require the &
            specification of two state labels'
       call error_control
    endif

    if (ifunc.eq.3.and.Jval.eq.0) then
       errmsg='The value of J has not been given'
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

  subroutine rdpotfile

    use constants
    use channels
    use iomod
    use sysinfo
    use symmetry
    use parameters
    
    implicit none

    integer  :: nmodes1,ndat1
    real(dp) :: freq1(nmodes)
    logical  :: ldip1
    
!----------------------------------------------------------------------
! Open the parameter file
!----------------------------------------------------------------------
    call freeunit(ibin)
    open(ibin,file=abin,form='unformatted',status='old')

!----------------------------------------------------------------------
! System dimensions
!----------------------------------------------------------------------
    ! Number of modes
    read(ibin) nmodes1

    ! Number of states
    read(ibin) nsta

    ! Exit if the number of modes does not match the number read from
    ! the freqency calculation output file
    if (nmodes1.ne.nmodes) then
       errmsg='Inconsistent normal mode numbers in rdpotfile'
       call error_control
    endif

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    ! Model diabatic potential arrays
    allocate(e0(nsta))
    allocate(kappa(nmodes,nsta))
    allocate(lambda(nmodes,nsta,nsta))
    allocate(gamma(nmodes,nmodes,nsta))
    allocate(mu(nmodes,nmodes,nsta,nsta))
    allocate(iota(nmodes,nsta))
    allocate(tau(nmodes,nsta,nsta))
    allocate(epsilon(nmodes,nsta))
    allocate(xi(nmodes,nsta,nsta))
    allocate(kappa_mask(nmodes,nsta))
    allocate(lambda_mask(nmodes,nsta,nsta))
    allocate(gamma_mask(nmodes,nmodes,nsta))
    allocate(mu_mask(nmodes,nmodes,nsta,nsta))
    allocate(iota_mask(nmodes,nsta))
    allocate(tau_mask(nmodes,nsta,nsta))
    allocate(epsilon_mask(nmodes,nsta))
    allocate(xi_mask(nmodes,nsta,nsta))
    
!----------------------------------------------------------------------
! Skip past variables that we do not need
!----------------------------------------------------------------------
    read(ibin) ndat1
    read(ibin) ldip1

!----------------------------------------------------------------------
! Vertical excitation energies
!----------------------------------------------------------------------
    read(ibin) e0
    
!----------------------------------------------------------------------
! Frequencies (read into a dummy array - we already have this
! information)
!----------------------------------------------------------------------
    read(ibin) freq1

!----------------------------------------------------------------------
! Coupling coefficients
!----------------------------------------------------------------------
    read(ibin) kappa
    read(ibin) lambda
    read(ibin) gamma
    read(ibin) mu
    read(ibin) iota
    read(ibin) tau
    read(ibin) epsilon
    read(ibin) xi
    
!----------------------------------------------------------------------
! Masks
!----------------------------------------------------------------------
    read(ibin) kappa_mask
    read(ibin) lambda_mask
    read(ibin) gamma_mask
    read(ibin) mu_mask
    read(ibin) iota_mask
    read(ibin) tau_mask
    read(ibin) epsilon_mask
    read(ibin) xi_mask

!----------------------------------------------------------------------
! Close the parameter file
!----------------------------------------------------------------------
    close(ibin)
    
    return
    
  end subroutine rdpotfile
    
!######################################################################
  
end program gridfunc
  

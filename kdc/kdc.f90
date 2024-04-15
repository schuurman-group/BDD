!######################################################################
! KDC: A program to compute coupling coefficients of the KDC vibronic
!      coupling Hamiltonian using the output of the blockdiag code.
!######################################################################

program kdc

  use ioqc
  use symmetry
  use parinfo
  use transform
  use opermod
  use parameters
  use kdcglobal
  
  implicit none

!----------------------------------------------------------------------
! Open the input and log files
!----------------------------------------------------------------------
  call open_files_kdc

!----------------------------------------------------------------------
! Determine the no. states from the input file
!----------------------------------------------------------------------
  call get_nsta_kdc
  
!----------------------------------------------------------------------
! Read the input file
!----------------------------------------------------------------------
  call rdinpfile_kdc

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
! Read the names of the BDD output files
!----------------------------------------------------------------------
  call rdfilenames

!----------------------------------------------------------------------
! Determine the BDD output file type
!----------------------------------------------------------------------
  call get_iprog
  
!----------------------------------------------------------------------
! Parse the BDD output files
!----------------------------------------------------------------------
  select case(iprog)
  case(1) ! blockdiag
     call parse_bdfiles
  case(2) ! GRaCI
     call parse_gracifiles
  end select
  
!----------------------------------------------------------------------
! Determine which coupling coefficients are zero by symmetry
!----------------------------------------------------------------------
  call create_mask(order1,ldipfit)

!----------------------------------------------------------------------
! Determine which pairs of modes give rise to non-zero coupling
! coefficients by symmetry
!----------------------------------------------------------------------
  call get_cutmask

!----------------------------------------------------------------------
! Optional shifts of the diabatic potential matrix elements
!----------------------------------------------------------------------
  if (lshift) call shift
  
!----------------------------------------------------------------------
! Optional block diagonalising transformation of the diabatic
! potential
!----------------------------------------------------------------------
  if (lblockdiag) call blockdiag
  
!----------------------------------------------------------------------
! Determine the coupling coefficients
!----------------------------------------------------------------------
  call get_coefficients
  
!----------------------------------------------------------------------
! Check to see if symmetry constraints are being satisfied
!----------------------------------------------------------------------
  call check_coefficients

!----------------------------------------------------------------------
! If the coupling coefficients were determined by fitting, then
! output the RMSDs to the log file
!----------------------------------------------------------------------
  if (ialgor.eq.2) call wrrmsd
  
!----------------------------------------------------------------------
! Output some useful information about the coupling coefficients to
! the log file
!----------------------------------------------------------------------
  call get_parinfo
  
!----------------------------------------------------------------------
! Write the MCTDH operator file
!----------------------------------------------------------------------
  call wroper

!----------------------------------------------------------------------
! Write the parameters to a binary file for reading by the pltkdc
! code
!----------------------------------------------------------------------
  call wrbinfile

!----------------------------------------------------------------------
! If a block diagonalisation of the diabatic potential was performed,
! then write the transformation matrices to disk for inspection
!----------------------------------------------------------------------
! IMPORTANT: write_Tmat must be called *after* the calculation of the
!            coupling coefficients due to the need for the findx1m
!            findx2m arrays to be filled
!----------------------------------------------------------------------
!  if (lblockdiag) call write_Tmat
  
contains

!######################################################################

  subroutine open_files_kdc
    
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
! Read the name of the input file and set the names of the log file
! operator file, and binary file
!----------------------------------------------------------------------
    call getarg(1,ain)

    ilbl=index(ain,'.inp')
    if (ilbl.eq.0) then
       alog=trim(ain)//'.log'
       aop=trim(ain)//'.op'
       abin=trim(ain)//'.dat'
       ain=trim(ain)//'.inp'
    else
       alog=ain(1:ilbl-1)//'.log'
       aop=ain(1:ilbl-1)//'.op'
       abin=ain(1:ilbl-1)//'.dat'
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
! Open the input, log, operator and binary files
!----------------------------------------------------------------------
    iin=1
    open(iin,file=ain,form='formatted',status='old')
    
    ilog=2
    open(ilog,file=alog,form='formatted',status='unknown')

    iop=3
    open(iop,file=aop,form='formatted',status='unknown')

    ibin=4
    open(ibin,file=abin,form='unformatted',status='unknown')
    
    return
    
  end subroutine open_files_kdc

!######################################################################
  
  subroutine get_nsta_kdc

    use constants
    use channels
    use iomod
    use parsemod
    use sysinfo
    use kdcglobal
    
    implicit none

!----------------------------------------------------------------------
! Determine the number of states from the input file
!----------------------------------------------------------------------
    rewind(iin)

    ! Determine the no. states from the $q0_ener section
    nsta=0
5   call rdinp(iin)
    if (keyword(1).eq.'$q0_ener') then
10     call rdinp(iin)
       if (keyword(1).ne.'$end') then
          nsta=nsta+1
          goto 10
       endif
    endif
    if (.not.lend) goto 5

    ! Quit if the $dets_ref section is missing
    if (nsta.eq.0) then
       errmsg='The $q0_ener section could not be found in '//trim(ain)
       call error_control
    endif

    return
    
  end subroutine get_nsta_kdc

!######################################################################

  subroutine rdinpfile_kdc

    use constants
    use channels
    use iomod
    use parsemod
    use sysinfo
    use parameters
    use kdcglobal

    implicit none

    integer :: i,ii,j,k,k1,k2
    
!----------------------------------------------------------------------
! Set defaults
!----------------------------------------------------------------------
    ! Frequency file
    freqfile=''

    ! Set file
    setfile=''
    lsetfile=.false.

    ! Q0 adiabatic energies
    allocate(q0pot(nsta))
    q0pot=0.0d0

    ! Zeroth-order parameter shifts
    lshift=.false.
    allocate(shift0(nsta,nsta))
    shift0=0.0d0
    
    ! Point group
    pntgrp=''

    ! Electronic state symmetry labels
    allocate(stalab(nsta))
    stalab=''

    ! Dipole symmetry labels
    diplab=''

    ! Q0 dipole matrix
    allocate(dip0(nsta,nsta,3))
    dip0=0.0d0

    ! Parameterisation algorithm: 1 <-> finite differences
    !                             2 <-> normal equations-based fitting
    ialgor=2

    ! Normal equation fitting weights
    iweight=0
    wfac=0.0d0

    ! GRaCI label stem to use
    lblstem=''

    ! Block diagonalisation
    lblockdiag=.false.

    ! Order of the 1-mode expansions
    order1=4
    
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

       else if (keyword(i).eq.'$bdfiles') then
          if (keyword(i+1).eq.'=') then
             ! Filenames are to be read from a set file
             lsetfile=.true.
             i=i+2
             setfile=keyword(i)
          else
             ! Filenames are given directly in the input file
             ! skip past for now and read these later
25           call rdinp(iin)
             if (keyword(1).ne.'$end') goto 25
          endif

       else if (keyword(i).eq.'$q0_ener') then
          do 
             call rdinp(iin)
             if (keyword(1).eq.'$end') exit
             if (lend) then
                errmsg='End of file reached whilst reading the &
                     $q0_ener section'
                call error_control
             endif
             read(keyword(2),*) k
             read(keyword(1),*) q0pot(k)
          enddo

       else if (keyword(i).eq.'$shifts') then
          lshift=.true.
          do
             call rdinp(iin)
             if (keyword(1).eq.'$end') exit
             if (lend) then
                errmsg='End of file reached whilst reading the &
                     $shifts section'
                call error_control
             endif
             read(keyword(2),*) k1
             read(keyword(3),*) k2
             read(keyword(1),*) shift0(k1,k2)
             shift0(k2,k1)=shift0(k1,k2)
          enddo
          
       else if (keyword(i).eq.'$point_group') then
          if (keyword(i+1).eq.'=') then
             i=i+2
             pntgrp=keyword(i)
          else
             goto 100
          endif

       else if (keyword(i).eq.'$state_sym') then
          do 
             call rdinp(iin)
             if (keyword(1).eq.'$end') exit
             if (lend) then
                errmsg='End of file reached whilst reading the &
                     $state_sym section'
                call error_control
             endif
             read(keyword(2),*) k
             stalab(k)=keyword(1)
          enddo

       else if (keyword(i).eq.'$algorithm') then
          if (keyword(i+1).eq.'=') then
             i=i+2
             if (keyword(i).eq.'fd') then
                ! Finite differences
                ialgor=1
             else if (keyword(i).eq.'normal') then
                ! Normal equations approach
                ialgor=2
             else
                errmsg='Unkown parameterisation algorithm: '&
                     //trim(keyword(i))
                call error_control
             endif                
          else
             goto 100
          endif

       else if (keyword(i).eq.'$dip_sym') then
          ldipfit=.true.
          do
             call rdinp(iin)
             if (keyword(1).eq.'$end') exit
             if (lend) then
                errmsg='End of file reached whilst reading the &
                     $dip_sym section'
                call error_control
             endif
             if (keyword(2).eq.'x') then
                k=1
             else if (keyword(2).eq.'y') then
                k=2
             else if (keyword(2).eq.'z') then
                k=3
             else
                errmsg='Unknown dipole component: '//trim(keyword(2))
                call error_control
             endif
             diplab(k)=keyword(1)
          enddo

       else if (keyword(i).eq.'$q0_dipole') then
          ldipfit=.true.
          do 
             call rdinp(iin)
             if (keyword(1).eq.'$end') exit
             if (lend) then
                errmsg='End of file reached whilst reading the &
                     $q0_dipole section'
                call error_control
             endif
             read(keyword(4),*) k1
             read(keyword(5),*) k2
             read(keyword(1),*) dip0(k1,k2,1)
             read(keyword(2),*) dip0(k1,k2,2)
             read(keyword(3),*) dip0(k1,k2,3)
             dip0(k2,k1,:)=dip0(k1,k2,:)
          enddo

       else if (keyword(i).eq.'$weight') then
          iweight=1
          wfac=1.0d0
          if (keyword(i+1).eq.'=') then
             i=i+2
             read(keyword(i),*) wfac
          endif

       else if (keyword(i).eq.'$label') then
          if (keyword(i+1).eq.'=') then
             i=i+2
             read(keyword(i),'(a)') lblstem
          else
             goto 100
          endif

       else if (keyword(i).eq.'$blockdiag') then
          if (keyword(i+1).eq.'=') then
             lblockdiag=.true.
             i=i+2
             ii=i
             if (keyword(i).eq.'{') then
                ! (1) Determine the no. states in each block
                nbd=0
                k=1
                do
                   i=i+1
                   if (keyword(i).eq.'}' .and. k.eq.2) exit
                   if (keyword(i).eq.'{') k=2
                   if (keyword(i).ne.',' &
                        .and. keyword(i).ne.'}' &
                        .and. keyword(i).ne.'{') then
                      nbd(k)=nbd(k)+1
                   endif
                enddo
                ! (2) Allocate arrays
                maxbd=maxval(nbd)
                allocate(ibd(maxbd,2))
                ibd=0
                ! (3) Read the indices of the states in each block
                k=1
                i=ii
                j=1
                do
                   i=i+1
                   if (keyword(i).eq.'}' .and. k.eq.2) exit
                   if (keyword(i).eq.'{') then
                      k=2
                      j=1
                   endif
                   if (keyword(i).ne.',' &
                        .and. keyword(i).ne.'}' &
                        .and. keyword(i).ne.'{') then
                      read(keyword(i),*) ibd(j,k)
                      j=j+1
                   endif
                enddo
             else
                errmsg='Error parsing the $blockdiag argument: '&
                     //'expected {i1,...,im} {j1,...,jn}'
                call error_control
             endif
          else
             goto 100
          endif

       else if (keyword(i).eq.'$order') then
          if (keyword(i+1).eq.'=') then
             i=i+2
             read(keyword(i),*) order1
          else
             goto 100
          endif
             
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

    if (pntgrp.eq.'') then
       errmsg='The point group has not been given'
       call error_control
    endif

    do i=1,nsta
       if (stalab(i).eq.'') then
          errmsg='Not all state symmetries have been given'
          call error_control
       endif
    enddo

    if (ldipfit) then
       do i=1,3
          if (diplab(i).eq.'') then
             errmsg='Not all dipole component symmetries have been &
                  given'
             call error_control
          endif
       enddo
       if (sum(abs(dip0)).eq.0.0d0) then
          errmsg='The Q0 dipole matrix elements have not been given'
          call error_control
       endif
    endif
    
    return
    
  end subroutine rdinpfile_kdc
    
!######################################################################

  subroutine rdfilenames

    use kdcglobal
    
    implicit none

    if (lsetfile) then
       ! Read the blockdiag filenames from a set files
       call rdfilenames_setfile
    else
       ! Read the blockdiag filenames from the input file
       call rdfilenames_inpfile
    endif
    
    return
    
  end subroutine rdfilenames

!######################################################################

  subroutine rdfilenames_setfile

    use constants
    use channels
    use parsemod
    use iomod
    use kdcglobal
    
    implicit none

    integer :: unit,ierr,i

!----------------------------------------------------------------------
! Open the set file
!----------------------------------------------------------------------
    call freeunit(unit)
    open(unit,file=setfile,form='formatted',status='old',iostat=ierr)

    if (ierr.ne.0) then
       errmsg='Error opening the set file: '//trim(setfile)
       call error_control
    endif
    
!----------------------------------------------------------------------
! First pass: determine the no. files and allocate arrays
!----------------------------------------------------------------------
    ! Determine the no. files
    nfiles=0
5   call rdinp(unit)
    if (.not.lend) then
       nfiles=nfiles+1
       goto 5
    endif

    ! Allocate arrays
    allocate(bdfiles(nfiles))
    bdfiles=''
    
!----------------------------------------------------------------------
! Second pass: read in the filenames
!----------------------------------------------------------------------
    rewind(unit)

    do i=1,nfiles
       call rdinp(unit)
       bdfiles(i)=keyword(1)
    enddo
    
!----------------------------------------------------------------------
! Close the set file
!----------------------------------------------------------------------
    close(unit)
    
    return
    
  end subroutine rdfilenames_setfile
    
!######################################################################

  subroutine rdfilenames_inpfile

    use constants
    use channels
    use parsemod
    use iomod
    use kdcglobal
    
    implicit none

    integer :: i
    
!----------------------------------------------------------------------
! Read to the bdfiles section
!----------------------------------------------------------------------
    rewind(iin)
5   call rdinp(iin)
    if (lend) goto 100
    if (keyword(1).ne.'$bdfiles') goto 5
       
!----------------------------------------------------------------------
! First pass: determine the no. files and allocate arrays
!----------------------------------------------------------------------
    ! Determine the no. files
    nfiles=0
10  call rdinp(iin)
    if (keyword(1).ne.'$end') then
       nfiles=nfiles+1
       goto 10
    endif

    ! Allocate arrays
    allocate(bdfiles(nfiles))
    bdfiles=''

!----------------------------------------------------------------------
! Second pass: read in the filenames
!----------------------------------------------------------------------
    do i=1,nfiles+1
       backspace(iin)
    enddo

    do i=1,nfiles
       call rdinp(iin)
       bdfiles(i)=keyword(1)
    enddo
       
    return

100 continue
    errmsg='The $bdfiles section could not be found'
    call error_control
    
  end subroutine rdfilenames_inpfile

!######################################################################

  subroutine get_iprog

    use constants
    use iomod
    use kdcglobal
    
    implicit none

    integer            :: unit,ierr
    character(len=120) :: string
    
!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    iprog=-1

!----------------------------------------------------------------------    
! Open the first output file
!----------------------------------------------------------------------
    call freeunit(unit)
    open(unit,file=bdfiles(1),form='formatted',status='old')
    
!----------------------------------------------------------------------
! Determine the file type
!----------------------------------------------------------------------
10  read(unit,'(a)',end=20) string

    if (index(string, 'B L O C K D I A G') /= 0) then
       ! blockdiag output
       iprog=1
       goto 20
    else if (index(string, 'GRaCI') /= 0) then
       ! GRaCI output
       iprog=2
       goto 20
    endif
    
    goto 10
    
20  continue

!----------------------------------------------------------------------    
! Close the first output file
!----------------------------------------------------------------------
    close(unit)

!----------------------------------------------------------------------
! Exit if the program type was not recognised
!----------------------------------------------------------------------
    if (iprog == -1) then
       errmsg='Error: BDD program type not recognised'
       call error_control
    endif
    
    return
    
  end subroutine get_iprog
  
!######################################################################

  subroutine parse_bdfiles

    use constants
    use iomod
    use sysinfo
    use kdcglobal
    
    implicit none

    integer :: n,unit,ierr
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    ! Number of displaced geometries: equal to the number of output
    ! files for blockdiag calculations
    ngeom=nfiles

    ! Diabatic potential matrix at the displaced geometries
    allocate(diabpot(nsta,nsta,ngeom))
    diabpot=0.0d0
    
    ! Normal mode coordinates at the displaced geometries
    allocate(qvec(nmodes,ngeom))
    qvec=0.0d0

    ! Diabatic dipole matrix elements at the displaced geometries
    if (ldipfit) then
       allocate(diabdip(nsta,nsta,3,ngeom))
       diabdip=0.0d0
    endif
    
!----------------------------------------------------------------------
! Parse the blockdiag log files
!----------------------------------------------------------------------
    do n=1,nfiles

       ! Open the blockdiag log file
       call freeunit(unit)
       open(unit,file=bdfiles(n),form='formatted',status='old',&
            iostat=ierr)
       if (ierr.ne.0) goto 999          

       ! Parse the blockdiag log file
       call parse_bdfiles_1file(unit,n)

       ! Close the blockdiag log file
       close(unit)
       
    enddo

    return

999 continue
    errmsg='Error opening the file: '//trim(bdfiles(n))
    call error_control

  end subroutine parse_bdfiles

!######################################################################

  subroutine parse_bdfiles_1file(unit,n)

    use constants
    use iomod
    use kdcglobal
    
    implicit none

    integer, intent(in) :: unit,n

!----------------------------------------------------------------------
! Point in normal modes
!----------------------------------------------------------------------
    call get_qvec_1file(unit,n)

!----------------------------------------------------------------------
! Diabatic potential matrix
!----------------------------------------------------------------------
    call get_diabpot_1file(unit,n)

!----------------------------------------------------------------------
! Diabatic dipole matrix
!----------------------------------------------------------------------
    if (ldipfit) call get_diabdip_1file(unit,n)
    
    return
    
  end subroutine parse_bdfiles_1file

!######################################################################

  subroutine parse_gracifiles

    use constants
    use iomod
    use sysinfo
    use kdcglobal
    
    implicit none

    integer :: n,unit,ierr
    integer :: igeom
    
!----------------------------------------------------------------------
! Determine the total number of displaced geometries
!----------------------------------------------------------------------
    call get_ngeom_graci

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    ! Diabatic potential matrix at the displaced geometries
    allocate(diabpot(nsta,nsta,ngeom))
    diabpot=0.0d0
    
    ! Normal mode coordinates at the displaced geometries
    allocate(qvec(nmodes,ngeom))
    qvec=0.0d0

    ! Diabatic dipole matrix elements at the displaced geometries
    if (ldipfit) then
       allocate(diabdip(nsta,nsta,3,ngeom))
       diabdip=0.0d0
    endif

!----------------------------------------------------------------------
! Parse the GRaCI output files
!----------------------------------------------------------------------
    call freeunit(unit)

    ! Displaced geometry counter
    igeom=0
    
    ! Loop over output files
    do n=1,nfiles

       ! Open the GRaCI output file
       open(unit,file=bdfiles(n),form='formatted',status='old',&
            iostat=ierr)
       if (ierr.ne.0) goto 999          

       ! Parse the GRaCI output file
       ! To do: add diabatic dipole matrices
       call parse_gracifiles_1file(unit,igeom)

       ! Close the blockdiag log file
       close(unit)
       
    enddo
    
    return

999 continue
    errmsg='Error opening the file: '//trim(bdfiles(n))
    call error_control
    
  end subroutine parse_gracifiles

!######################################################################

  subroutine get_ngeom_graci

    use constants
    use iomod
    use kdcglobal
    
    implicit none

    integer            :: n,unit,il,i1,i
    character(len=120) :: string

!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    ngeom=0
    
!----------------------------------------------------------------------
! No label stem given: count the number of diabatic potential matrix
! sections
!----------------------------------------------------------------------
    if (lblstem == '') then
    
       call freeunit(unit)

       ! Loop over files
       do n=1,nfiles
          
          ! Open the output file
          open(unit,file=bdfiles(n),form='formatted',status='old')

          ! Search for diabatic potential sections
10        read(unit,'(a)',end=20) string
          if (index(string,'Diabatic potential matrix elements') /= 0) &
            ngeom=ngeom+1
          goto 10

20        continue
       
          ! Close the output file
          close(unit)
       
       enddo

    endif

!----------------------------------------------------------------------
! Label stem given: count the number of diabatic potential matrix
! sections for DFT/MRCI(2) calculations with the correct label stem
!----------------------------------------------------------------------
    if (lblstem /= '') then

       ! Length of the label stem
       il=len_trim(lblstem)

       call freeunit(unit)
       
       ! Loop over files
       do n=1,nfiles
       
       ! Open the ouput file
          open(unit,file=bdfiles(1),form='formatted',status='old')
       
          ! Search for the DFT/MRCI(2) sections with the correct labels
30        read(unit,'(a)',end=40) string
          if (index(string,'DFT/MRCI(2) computation') /= 0) then
             i1=index(string,'=')
             if (string(i1+2:i1+1+il) == trim(lblstem)) then
                
                ! If we are here, then we are in a DFT/MRCI(2) section
                ! with the correct label stem. Now look for a diabatic
                ! potential matrix section
                
                do i=1,3
                   read(unit,*)
                enddo
             
50              read(unit,'(a)') string
                
                if (index(string,'Diabatic potential matrix elements') &
                     /= 0) ngeom=ngeom+1
                
                if (index(string,'**********') == 0) goto 50

             
             endif
          endif
       
          goto 30
          
40        continue
    
          ! Close the output file
          close(unit)

       enddo
       
    endif
       
    return
    
  end subroutine get_ngeom_graci

!######################################################################

  subroutine parse_gracifiles_1file(unit,igeom)

    use constants
    use iomod
    use sysinfo
    use kdcglobal
    
    implicit none

    integer, intent(in)         :: unit
    integer, intent(inout)      :: igeom
    integer                     :: i,j,itmp,jtmp,i1,il
    real(dp), dimension(ncoo)   :: xcoo
    real(dp), dimension(nmodes) :: qcoo
    real(dp), parameter         :: q0thrsh=1e-4_dp
    character(len=120)          :: string
    character(len=2)            :: atmp
    logical                     :: lstem
    
!----------------------------------------------------------------------
! Read in the points in normal mode coordinates and the corresponding
! diabatic potential values
!----------------------------------------------------------------------
    rewind(unit)

10  read(unit,'(a)',end=100) string

    ! Are we at a new DFT/MRCI(2) section?
    if (index(string,'DFT/MRCI(2) computation') /= 0) then

       ! Does this DFT/MRCI(2) section have the correct label stem?
       lstem=.false.
       if (lblstem == '') then
          lstem=.true.
       else
          il=len_trim(lblstem)
          i1=index(string,'=')
          if (string(i1+2:i1+1+il) == trim(lblstem)) lstem=.true.
       endif

       if (lstem) then

          ! Read to the Cartesian coordinates
15        read(unit,'(a)',end=100) string
          if (index(string,'Cartesian Coordinates') /= 0) then
             
             ! Parse the Cartesian coordinates
             read(unit,*)
             do i=1,natm
                read(unit,*) atmp,(xcoo(j),j=i*3-2,i*3)
             enddo
             
             ! Normal mode coordinates
             qcoo=matmul(coonm,(xcoo-xcoo0/ang2bohr))

             ! Skip if we are at Q0...
             if (sqrt(dot_product(qcoo,qcoo)) < q0thrsh) then
                goto 10
             endif
          
             ! ...otherwise save the normal mode coordinates and carry on
             igeom=igeom+1
             qvec(:,igeom)=qcoo
             
             ! Read to the diabatic potential section for this geometry
20           read(unit,'(a)',end=100) string
             if (index(string,'Diabatic potential matrix elements') == 0) &
                  goto 20
             
             ! Parse the diabatic potential (C1 symmetry assumed)
             do i=1,3
                read(unit,*)
             enddo
             do i=1,nsta
                do j=i,nsta
                   read(unit,*) itmp,jtmp,diabpot(i,j,igeom)
                   diabpot(j,i,igeom)=diabpot(i,j,igeom)
                enddo
             enddo
             
          else
             goto 15
          endif

       endif
          
    endif
    
    goto 10

100 continue
    
    return
    
  end subroutine parse_gracifiles_1file
    
!######################################################################

  subroutine get_qvec_1file(unit,n)

    use constants
    use iomod
    use ioqc
    use sysinfo
    use kdcglobal

    implicit none

    integer, intent(in)       :: unit,n
    integer                   :: i,j
    real(dp), dimension(ncoo) :: xcoo
    character(len=120)        :: string
    character(len=2)          :: atmp
    
!----------------------------------------------------------------------
! Read in the Cartesian coordinates (in Angstrom)
!----------------------------------------------------------------------
    rewind(unit)
5   read(unit,'(a)',end=999) string
    if (index(string,'Displaced geometry:').eq.0) goto 5

    do i=1,natm
       read(unit,'(1x,a2,3(2x,F10.7))') atmp,(xcoo(j),j=i*3-2,i*3)
    enddo
    
!----------------------------------------------------------------------
! Normal mode coordinates
! Note that xcoo0 is in Bohr, whilst we have read in xcoo in Angstrom
!----------------------------------------------------------------------
    qvec(:,n)=matmul(coonm,(xcoo-xcoo0/ang2bohr))

    return

999 continue
    errmsg='The displaced geometry section could not be found in: '&
         //trim(bdfiles(n))
    call error_control
    
  end subroutine get_qvec_1file

!######################################################################
  
  subroutine get_diabpot_1file(unit,n)

    use constants
    use iomod
    use ioqc
    use sysinfo
    use kdcglobal

    implicit none

    integer, intent(in) :: unit,n
    integer             :: i,j
    character(len=120)  :: string

!----------------------------------------------------------------------
! Read in the diabatic potential matrix (in a.u.)
!----------------------------------------------------------------------
    rewind(unit)
5   read(unit,'(a)',end=999) string
    if (index(string,'Quasi-Diabatic Potential Matrix').eq.0) goto 5
    read(unit,*)
    
    do i=1,nsta
       do j=i,nsta
          read(unit,'(10x,F15.10)') diabpot(i,j,n)
          diabpot(j,i,n)=diabpot(i,j,n)
       enddo
    enddo

    return

999 continue
    errmsg='The diabatic potential section could not be found in: '&
         //trim(bdfiles(n))
    call error_control

  end subroutine get_diabpot_1file

!######################################################################

  subroutine get_diabdip_1file(unit,n)

    use constants
    use iomod
    use ioqc
    use sysinfo
    use kdcglobal

    implicit none

    integer, intent(in) :: unit,n
    integer             :: i,j
    character(len=120)  :: string

!----------------------------------------------------------------------
! Read in the diabatic dipole matrix (in a.u.)
!----------------------------------------------------------------------
    rewind(unit)
5   read(unit,'(a)',end=999) string
    if (index(string,'Quasi-Diabatic Dipole Matrix').eq.0) goto 5
    read(unit,*)
    read(unit,*)

    do i=1,nsta
       do j=i,nsta
          read(unit,'(8x,3(2x,F15.10))') diabdip(i,j,:,n)
          diabdip(j,i,:,n)=diabdip(i,j,:,n)
       enddo
    enddo
    
    return

999 continue
    errmsg='The diabatic dipole section could not be found in: '&
         //trim(bdfiles(n))
    call error_control
    
  end subroutine get_diabdip_1file
    
!######################################################################

  subroutine get_coefficients

    use constants
    use sysinfo
    use parameters
    use kdcglobal
    use fdmod
    use nmeqmod
    
    implicit none

!----------------------------------------------------------------------
! Allocate and initialise arrays
!----------------------------------------------------------------------
    ! One-mode terms
    allocate(coeff1(nmodes,nsta,nsta,order1))
    coeff1=0.0d0

    ! Two-mode terms (2nd-order only)
    allocate(coeff2(nmodes,nmodes,nsta,nsta))
    coeff2=0.0d0
    
    ! First-order, intrastate
    allocate(kappa(nmodes,nsta))
    kappa=0.0d0

    ! First-order, interstate
    allocate(lambda(nmodes,nsta,nsta))
    lambda=0.0d0

    ! Second-order, intrastate
    allocate(gamma(nmodes,nmodes,nsta))
    gamma=0.0d0

    ! Second-order, interstate
    allocate(mu(nmodes,nmodes,nsta,nsta))
    mu=0.0d0

    ! Third-order, cubic, intrastate
    allocate(iota(nmodes,nsta))
    iota=0.0d0

    ! Third-order, cubic, interstate
    allocate(tau(nmodes,nsta,nsta))
    tau=0.0d0

    ! Fourth-order, quartic, intrastate
    allocate(epsilon(nmodes,nsta))
    epsilon=0.0d0

    ! Fourth-order, quartic, interstate
    allocate(xi(nmodes,nsta,nsta))
    xi=0.0d0

    ! Diabatic dipole matrix expansion coefficients
    if (ldipfit) then
       allocate(dip1(nmodes,nsta,nsta,3))
       dip1=0.0d0
       allocate(dip2(nmodes,nmodes,nsta,nsta,3))
       dip2=0.0d0
       allocate(dip3(nmodes,nsta,nsta,3))
       dip3=0.0d0
       allocate(dip4(nmodes,nsta,nsta,3))
       dip4=0.0d0
    endif
    
!----------------------------------------------------------------------
! Vertical excitation energies
!----------------------------------------------------------------------
    allocate(e0(nsta))
    e0(:)=(q0pot(:)-q0pot(1))*eh2ev
    
!----------------------------------------------------------------------
! Calculate the coupling coefficients
!----------------------------------------------------------------------
    if (ialgor.eq.1) then
       ! Finite differences
       call get_coefficients_fd
    else if (ialgor.eq.2) then
       ! Normal equations fitting
       call get_coefficients_nmeq
    endif
    
    return
    
  end subroutine get_coefficients
      
!######################################################################

  subroutine check_coefficients

    use constants
    use channels
    use iomod
    use sysinfo
    use parameters
    use kdcglobal
    use symmetry
    
    implicit none

    integer                        :: m,m1,m2,s,s1,s2,c,i
    real(dp), parameter            :: thrsh=1e-4_dp
    character(len=1), dimension(3) :: acomp

!----------------------------------------------------------------------
! Cartesian axis labels
!----------------------------------------------------------------------
    acomp=(/ 'x' , 'y' , 'z' /)
    
!----------------------------------------------------------------------
! First, output some information about the number of parameters of
! each class that are non-zero by symmetry
!----------------------------------------------------------------------
    ! Section header
    write(ilog,'(/,72a)') ('+',i=1,72)
    write(ilog,'(2x,a)') 'Coupling Coefficient Information'
    write(ilog,'(72a)') ('+',i=1,72)

    ! Get the total and symmetry-allowed number of coupling
    ! coefficients
    call getnpar(ldipfit)

    ! Output the total and symmetry-allowed number of coupling
    ! coefficients to the log file
    write(ilog,'(/,42a)') ('-',i=1,42)
    write(ilog,'(a)') ' Type   | Total number | Symmetry allowed'
    write(ilog,'(42a)') ('-',i=1,42)

    write(ilog,'(a,3x,a,4x,i5,6x,a1,4x,i5)') &
         'kappa','|',nkappa(1),'|',nkappa(2)

    write(ilog,'(a,2x,a,4x,i5,6x,a1,4x,i5)') &
         'lambda','|',nlambda(1),'|',nlambda(2)

    write(ilog,'(a,3x,a,4x,i5,6x,a1,4x,i5)') &
         'gamma','|',ngamma(1),'|',ngamma(2)

    write(ilog,'(a,6x,a,4x,i5,6x,a1,4x,i5)') &
         'mu','|',nmu(1),'|',nmu(2)

    write(ilog,'(a,4x,a,4x,i5,6x,a1,4x,i5)') &
         'iota','|',niota(1),'|',niota(2)

    write(ilog,'(a,5x,a,4x,i5,6x,a1,4x,i5)') &
         'tau','|',ntau(1),'|',ntau(2)

    write(ilog,'(a,1x,a,4x,i5,6x,a1,4x,i5)') &
         'epsilon','|',nepsilon(1),'|',nepsilon(2)

    write(ilog,'(a,6x,a,4x,i5,6x,a1,4x,i5)') &
         'xi','|',nxi(1),'|',nxi(2)

    if (ldipfit) then
       write(ilog,'(a,4x,a,4x,i5,6x,a1,4x,i5)') &
         'dip1','|',ndip1(1),'|',ndip1(2)

       write(ilog,'(a,4x,a,4x,i5,6x,a1,4x,i5)') &
         'dip2','|',ndip2(1),'|',ndip2(2)

       write(ilog,'(a,4x,a,4x,i5,6x,a1,4x,i5)') &
         'dip3','|',ndip3(1),'|',ndip3(2)

       write(ilog,'(a,4x,a,4x,i5,6x,a1,4x,i5)') &
         'dip4','|',ndip4(1),'|',ndip4(2)
    endif
    
    write(ilog,'(a,5x,a,4x,i5,6x,a1,4x,i5)') &
         'all','|',ntot(1),'|',ntot(2)

    write(ilog,'(42a,/)') ('-',i=1,42)
    
!----------------------------------------------------------------------
! Check for any non-zero coupling coefficients that should be zero
! by symmetry
!----------------------------------------------------------------------
    ! kappa
    do s=1,nsta
       do m=1,nmodes
          if (abs(kappa(m,s)).gt.thrsh.and.kappa_mask(m,s).eq.0) then
             write(ilog,'(a,2x,a,2(x,i2),2x,F10.7)') &
                  'WARNING non-zero parameter that should be zero &
                  by symmetry:','kappa',m,s,kappa(m,s)
          endif
       enddo
    enddo

    ! lambda
    do s1=1,nsta-1
       do s2=s1+1,nsta
          do m=1,nmodes
             if (abs(lambda(m,s1,s2)).gt.thrsh&
                  .and.lambda_mask(m,s1,s2).eq.0) then
                write(ilog,'(a,2x,a,3(x,i2),2x,F10.7)') &
                     'WARNING non-zero parameter that should be zero &
                     by symmetry:','lambda',m,s1,s2,lambda(m,s1,s2)
             endif
          enddo
       enddo
    enddo

    ! gamma
    do s=1,nsta
       do m1=1,nmodes
          do m2=m1,nmodes
             if (abs(gamma(m1,m2,s)).gt.thrsh&
                  .and.gamma_mask(m1,m2,s).eq.0) then
                write(ilog,'(a,2x,a,3(x,i2),2x,F10.7)') &
                     'WARNING non-zero parameter that should be zero &
                     by symmetry:','gamma',m1,m2,s,gamma(m1,m2,s)
             endif
          enddo
       enddo
    enddo

    ! mu
    do s1=1,nsta-1
       do s2=s1+1,nsta
          do m1=1,nmodes
             do m2=m1,nmodes
                if (abs(mu(m1,m2,s1,s2)).gt.thrsh&
                     .and.mu_mask(m1,m2,s1,s2).eq.0) then
                   write(ilog,'(a,2x,a,4(x,i2),2x,F10.7)') &
                        'WARNING non-zero parameter that should be &
                        zero by symmetry:','mu',m1,m2,s1,s2,&
                        mu(m1,m2,s1,s2)
                endif
             enddo
          enddo
       enddo
    enddo

    ! iota
    do s=1,nsta
       do m=1,nmodes
          if (abs(iota(m,s)).gt.thrsh&
               .and.iota_mask(m,s).eq.0) then
             write(ilog,'(a,2x,a,2(x,i2),2x,F10.7)') &
                  'WARNING non-zero parameter that should be &
                  zero by symmetry:','iota',m,s,&
                  iota(m,s)
          endif
       enddo
    enddo

    ! tau
    do s1=1,nsta-1
       do s2=s1+1,nsta
          do m=1,nmodes
             if (abs(tau(m,s1,s2)).gt.thrsh&
                  .and.tau_mask(m,s1,s2).eq.0) then
                write(ilog,'(a,2x,a,3(x,i2),2x,F10.7)') &
                     'WARNING non-zero parameter that should be &
                     zero by symmetry:','tau',m,s1,s2,&
                     tau(m,s1,s2)
             endif
          enddo
       enddo
    enddo

    ! espilon
    do s=1,nsta
       do m=1,nmodes
          if (abs(epsilon(m,s)).gt.thrsh&
               .and.epsilon_mask(m,s).eq.0) then
             write(ilog,'(a,2x,a,2(x,i2),2x,F10.7)') &
                  'WARNING non-zero parameter that should be &
                  zero by symmetry:','epsilon',m,s,&
                  epsilon(m,s)
          endif
       enddo
    enddo

    ! xi
    do s1=1,nsta-1
       do s2=s1+1,nsta
          do m=1,nmodes
             if (abs(xi(m,s1,s2)).gt.thrsh&
                  .and.xi_mask(m,s1,s2).eq.0) then
                write(ilog,'(a,2x,a,3(x,i2),2x,F10.7)') &
                     'WARNING non-zero parameter that should be &
                     zero by symmetry:','xi',m,s1,s2,&
                     xi(m,s1,s2)
             endif
          enddo
       enddo
    enddo

!----------------------------------------------------------------------
! Check for any non-zero diabatic dipole matrix expansion coefficients
! that should be zero by symmetry
!----------------------------------------------------------------------
    if (ldipfit) then

       ! 1st-order terms
       do c=1,3
          do s1=1,nsta
             do s2=s1,nsta
                do m=1,nmodes
                   if (abs(dip1(m,s1,s2,c)).gt.thrsh&
                        .and.dip1_mask(m,s1,s2,c).eq.0) then
                      write(ilog,'(a,2x,a,3(x,i2),x,a1,2x,F10.7)') &
                           'WARNING non-zero parameter that should be &
                           zero by symmetry:','dip1',m,s1,s2,acomp(c),&
                           dip1(m,s1,s2,c)
                   endif
                enddo
             enddo
          enddo
       enddo

       ! 2nd-order terms
       do c=1,3
          do s1=1,nsta
             do s2=s1,nsta
                do m1=1,nmodes
                   do m2=m1,nmodes
                      if (abs(dip2(m1,m2,s1,s2,c)).gt.thrsh &
                           .and.dip2_mask(m1,m2,s1,s2,c).eq.0) then
                         write(ilog,'(a,2x,a,4(x,i2),x,a1,2x,F10.7)') &
                              'WARNING non-zero parameter that should &
                              be zero by symmetry:','dip2',m1,m2,s1,s2,&
                              acomp(c),dip2(m1,m2,s1,s2,c)
                      endif
                   enddo
                enddo
             enddo
          enddo
       enddo

       ! 3rd-order terms
       do c=1,3
          do s1=1,nsta
             do s2=s1,nsta
                do m=1,nmodes
                   if (abs(dip3(m,s1,s2,c)).gt.thrsh&
                        .and.dip3_mask(m,s1,s2,c).eq.0) then
                      write(ilog,'(a,2x,a,3(x,i2),x,a1,2x,F10.7)') &
                           'WARNING non-zero parameter that should be &
                           zero by symmetry:','dip3',m,s1,s2,acomp(c),&
                           dip3(m,s1,s2,c)
                   endif
                enddo
             enddo
          enddo
       enddo

       ! 4th-order terms
       do c=1,3
          do s1=1,nsta
             do s2=s1,nsta
                do m=1,nmodes
                   if (abs(dip4(m,s1,s2,c)).gt.thrsh&
                        .and.dip4_mask(m,s1,s2,c).eq.0) then
                      write(ilog,'(a,2x,a,3(x,i2),x,a1,2x,F10.7)') &
                           'WARNING non-zero parameter that should be &
                           zero by symmetry:','dip4',m,s1,s2,acomp(c),&
                           dip4(m,s1,s2,c)
                   endif
                enddo
             enddo
          enddo
       enddo
       
    endif
    
    return
    
  end subroutine check_coefficients

!######################################################################

  subroutine wrrmsd

    use constants
    use channels
    use sysinfo
    use kdcglobal
    
    implicit none

    integer :: i,m

!----------------------------------------------------------------------
! Section header
!----------------------------------------------------------------------
    write(ilog,'(/,72a)') ('+',i=1,72)
    write(ilog,'(2x,a)') 'Fit RMSDs (Adiabatic Energies)'
    write(ilog,'(72a)') ('+',i=1,72)

!----------------------------------------------------------------------
! Total RMSD
!----------------------------------------------------------------------
    write(ilog,'(/,2x,a,2x,ES10.4,x,a)') 'Total RMSD:',rmsd,'eV'

!----------------------------------------------------------------------
! One-mode RMSDs
!----------------------------------------------------------------------
    write(ilog,'(/,2x,a)') 'One-Mode RMSDs:'
    write(ilog,'(20a)') ('-',i=1,20)
    write(ilog,'(a)') ' Mode | RMSD (eV)'
    write(ilog,'(20a)') ('-',i=1,20)
    do m=1,nmodes
       write(ilog,'(x,i3,2x,a,x,ES10.4)') m,'|',rmsd1m(m)
    enddo
    write(ilog,'(20a)') ('-',i=1,20)
    
    return
    
  end subroutine wrrmsd
  
!######################################################################

  subroutine wrbinfile

    use constants
    use channels
    use iomod
    use sysinfo
    use parameters
    use kdcglobal
    use symmetry
    
    implicit none

!----------------------------------------------------------------------
! System dimensions
!----------------------------------------------------------------------
    write(ibin) nmodes
    write(ibin) nsta

!----------------------------------------------------------------------
! No. ab initio diabatic potential values
!----------------------------------------------------------------------
    write(ibin) ngeom

!----------------------------------------------------------------------
! Diabatic dipole flag
!----------------------------------------------------------------------
    write(ibin) ldipfit
    
!----------------------------------------------------------------------
! Vertical excitation energies
!----------------------------------------------------------------------
    write(ibin) e0
    
!----------------------------------------------------------------------
! Frequencies
!----------------------------------------------------------------------
    write(ibin) freq

!----------------------------------------------------------------------
! Coupling coefficients
!----------------------------------------------------------------------
    write(ibin) kappa
    write(ibin) lambda
    write(ibin) gamma
    write(ibin) mu
    write(ibin) iota
    write(ibin) tau
    write(ibin) epsilon
    write(ibin) xi
    
!----------------------------------------------------------------------
! Masks
!----------------------------------------------------------------------
    write(ibin) kappa_mask
    write(ibin) lambda_mask
    write(ibin) gamma_mask
    write(ibin) mu_mask
    write(ibin) iota_mask
    write(ibin) tau_mask
    write(ibin) epsilon_mask
    write(ibin) xi_mask

!----------------------------------------------------------------------
! Ab initio diabatic potential values
!----------------------------------------------------------------------
    write(ibin) q0pot
    write(ibin) diabpot
    write(ibin) qvec
    
!----------------------------------------------------------------------
! Diabatic dipole coefficients and masks
!----------------------------------------------------------------------
    if (ldipfit) then

       ! Coefficients
       write(ibin) dip0
       write(ibin) dip1
       write(ibin) dip2
       write(ibin) dip3
       write(ibin) dip4

       ! Masks
       write(ibin) dip1_mask
       write(ibin) dip2_mask
       write(ibin) dip3_mask
       write(ibin) dip4_mask

       ! Ab initio diabatic dipole values
       write(ibin) diabdip
       
    endif
    
    return
    
  end subroutine wrbinfile
    
!######################################################################

end program kdc

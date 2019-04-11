!######################################################################
! KDC: A program to compute coupling coefficients of the KDC vibronic
!      coupling Hamiltonian using the output of the blockdiag code.
!######################################################################

program kdc

  use ioqc
  use symmetry
  use parinfo
  use opermod
  
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
! Read the names of the blockdiag files
!----------------------------------------------------------------------
  call rdbdfilenames
  
!----------------------------------------------------------------------
! Parse the blockdiag files
!----------------------------------------------------------------------
  call parse_bdfiles

!----------------------------------------------------------------------
! Determine which coupling coefficients are zero by symmetry
!----------------------------------------------------------------------
  call create_mask

!----------------------------------------------------------------------
! Determine which pairs of modes give rise to non-zero coupling
! coefficients by symmetry
!----------------------------------------------------------------------
  call get_cutmask

!----------------------------------------------------------------------
! Determine the coupling coefficients
!----------------------------------------------------------------------
  call get_coefficients
  


!----------------------------------------------------------------------
! Check to see if symmetry constraints are being satisfied
!----------------------------------------------------------------------
  call check_coefficients

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
    use kdcglobal

    implicit none

    integer :: i,k

!----------------------------------------------------------------------
! Set defaults
!----------------------------------------------------------------------
    freqfile=''
    setfile=''
    lsetfile=.false.

    allocate(q0pot(nsta))
    q0pot=0.0d0

    pntgrp=''
    
    allocate(stalab(nsta))
    stalab=''

    ialgor=1
    
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

    return
    
  end subroutine rdinpfile_kdc
    
!######################################################################

  subroutine rdbdfilenames

    use kdcglobal
    
    implicit none

    if (lsetfile) then
       ! Read the blockdiag filenames from a set files
       call rdbdfilenames_setfile
    else
       ! Read the blockdiag filenames from the input file
       call rdbdfilenames_inpfile
    endif
    
    return
    
  end subroutine rdbdfilenames

!######################################################################

  subroutine rdbdfilenames_setfile

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
    
  end subroutine rdbdfilenames_setfile
    
!######################################################################

  subroutine rdbdfilenames_inpfile

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
    
  end subroutine rdbdfilenames_inpfile

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
    ! Diabatic potential matrix at the displaced geometries
    allocate(diabpot(nsta,nsta,nfiles))
    diabpot=0.0d0

    ! Normal mode coordinates at the displaced geometries
    allocate(qvec(nmodes,nfiles))
    qvec=0.0d0
    
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
    
    return
    
  end subroutine parse_bdfiles_1file

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

  subroutine get_coefficients

    use constants
    use kdcglobal
    use fdmod
    
    implicit none

    if (ialgor.eq.1) then
       ! Finite differences
       call get_coefficients_fd
    else if (ialgor.eq.2) then
       ! Normal equations fitting
       print*,"HERE"
       stop
    endif
    
    return
    
  end subroutine get_coefficients
    
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
! Read in the Cartesian coordinates (in Angstrom)
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

  subroutine check_coefficients

    use constants
    use channels
    use iomod
    use sysinfo
    use kdcglobal    
    use symmetry
    
    implicit none

    integer             :: m,m1,m2,s,s1,s2,i
    real(dp), parameter :: thrsh=1e-4_dp

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
    call getnpar

    ! Output the total and symmetry-allowed number of coupling
    ! coefficients to the log file
    write(ilog,'(/,41a)') ('-',i=1,41)
    write(ilog,'(a)') ' Type  | Total number | Symmetry allowed'
    write(ilog,'(41a)') ('-',i=1,41)

    write(ilog,'(a,2x,a,4x,i4,6x,a1,4x,i4)') &
         'kappa','|',nkappa(1),'|',nkappa(2)

    write(ilog,'(a,1x,a,4x,i4,6x,a1,4x,i4)') &
         'lambda','|',nlambda(1),'|',nlambda(2)

    write(ilog,'(a,2x,a,4x,i4,6x,a1,4x,i4)') &
         'gamma','|',ngamma(1),'|',ngamma(2)

    write(ilog,'(a,5x,a,4x,i4,6x,a1,4x,i4)') &
         'mu','|',nmu(1),'|',nmu(2)

    write(ilog,'(a,4x,a,4x,i4,6x,a1,4x,i4)') &
         'all','|',ntot(1),'|',ntot(2)

    write(ilog,'(41a,/)') ('-',i=1,41)
    
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
       
    return
    
  end subroutine check_coefficients
    
!######################################################################

  subroutine wrbinfile

    use constants
    use channels
    use iomod
    use sysinfo
    use kdcglobal    
    use symmetry
    
    implicit none

    real(dp), dimension(nsta) :: e0
    
!----------------------------------------------------------------------
! System dimensions
!----------------------------------------------------------------------
    write(ibin) nmodes
    write(ibin) nsta

!----------------------------------------------------------------------
! Vertical excitation energies
!----------------------------------------------------------------------
    e0(:)=(q0pot(:)-q0pot(1))*eh2ev
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

!----------------------------------------------------------------------
! Masks
!----------------------------------------------------------------------
    write(ibin) kappa_mask
    write(ibin) lambda_mask
    write(ibin) gamma_mask
    write(ibin) mu_mask
    
    return
    
  end subroutine wrbinfile
    
!######################################################################

end program kdc

!######################################################################
! blockdiag: A program for the construction of quasidiabatic potential
!            matrices using the block diagonal diabatisation method
!            of Pacher, Cederbaum and Koppel.
!
!            See JCP, 89, 7367 (1988) for a decription of the method.
!
!            Uses wavefunctions expressed in terms of Slater
!            determinants, using the multigrid/superdyson format.
!
!            MOs are read from GAMESS checkpoint style files.
!######################################################################

program blockdiag

  use constants
  use channels
  use detparsing
  use mooverlaps
  use bdglobal
  use pacher
  use tamura
  use timingmod
  
  implicit none

  real(dp) :: tw1,tw2,tc1,tc2

!-----------------------------------------------------------------------
! Start timing
!-----------------------------------------------------------------------
    call times(tw1,tc1)
  
!----------------------------------------------------------------------
! Open the input and log files
!----------------------------------------------------------------------
  call open_files

!----------------------------------------------------------------------
! Write the log file header
!----------------------------------------------------------------------
  call wrheader
  
!----------------------------------------------------------------------
! Determine the number of electronic states
!----------------------------------------------------------------------
  call get_nsta
  
!----------------------------------------------------------------------
! Initialisation and allocation of arrays
!----------------------------------------------------------------------
  call initialise
  
!----------------------------------------------------------------------
! Read the input file
!----------------------------------------------------------------------
  call rdinpfile
  
!----------------------------------------------------------------------
! Read the displaced geometry adiabatic energies from a quantum
! chemistry ouput file if needed
!----------------------------------------------------------------------
  if (avmat.ne.'') call rdvmat
  
!----------------------------------------------------------------------
! Load the MOs for the reference and displaced geometries
!----------------------------------------------------------------------
  call loadmos

!----------------------------------------------------------------------
! Write the reference and displaced geometries to the log file
!----------------------------------------------------------------------
  call wrgeoms
  
!----------------------------------------------------------------------
! Read in the determinant files
!----------------------------------------------------------------------
  call rddetfiles

!----------------------------------------------------------------------
! Optional adjustment of the phases of the reference geometry
! wavefunctions. Note that this is only possible if a reference
! geometry blockdiag log file has been specified in the input file.
!----------------------------------------------------------------------
  if (lrdreftrans) call phase_refpsi
  
!----------------------------------------------------------------------
! Write the norms of the reference and displaced geometry
! wavefunctions to the log file
!----------------------------------------------------------------------
  call wrnorms
  
!----------------------------------------------------------------------
! Calculate the overlaps between the MOs at the reference and
! displaced geometries
!----------------------------------------------------------------------
  call mo_overlaps

!----------------------------------------------------------------------
! Perform the diabatisation
!----------------------------------------------------------------------
  if (algorithm.eq.'pacher') then
     call pacher_diabatisation
  else if (algorithm.eq.'tamura') then
     call tamura_diabatisation
  endif

!-----------------------------------------------------------------------
! Finalisation
!-----------------------------------------------------------------------
  call finalise
  
!-----------------------------------------------------------------------    
! Output timings
!-----------------------------------------------------------------------    
    call times(tw2,tc2)
    write(ilog,'(/,a,1x,F9.2,1x,a)') 'Wall Time:',tw2-tw1," s"
    write(ilog,'(a,2x,F9.2,1x,a)') 'CPU Time:',tc2-tc1," s"
    
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

  subroutine wrheader

    use channels
    use constants
    
    implicit none

    integer :: i

    write(ilog,'(82a)') ('-',i=1,82)
    write(ilog,'(32x,a)') 'B L O C K D I A G'
    write(ilog,'(82a)') ('-',i=1,82)
    write(ilog,'(20x,a)') 'A BLOCK DIAGONALISATION DIABATISATION CODE'
    write(ilog,'(82a,/)') ('-',i=1,82)
    
    return
    
  end subroutine wrheader
    
!######################################################################
  
  subroutine get_nsta

    use constants
    use channels
    use iomod
    use parsemod
    use bdglobal
    
    implicit none
    
!----------------------------------------------------------------------
! Determine the number of ref. geometry states from the input file
!----------------------------------------------------------------------
    rewind(iin)

    ! Determine the no. states from the $dets_ref section
    nsta_ref=0
5   call rdinp(iin)
    if (keyword(1).eq.'$dets_ref') then
10     call rdinp(iin)
       if (keyword(1).ne.'$end') then
          nsta_ref=nsta_ref+1
          goto 10
       endif
    endif
    if (.not.lend) goto 5

    ! Quit if the $dets_ref section is missing
    if (nsta_ref.eq.0) then
       errmsg='The $dets_ref section could not be found in '//trim(ain)
       call error_control
    endif
    
!----------------------------------------------------------------------
! Determine the number of disp. geometry states from the input file.
! Note that nsta_disp can be different to nsta_ref if we trying to
! follow states by character from geometry to geometry.
!----------------------------------------------------------------------
    rewind(iin)

    ! Determine the no. states from the $dets_disp section
    nsta_disp=0
15  call rdinp(iin)
    if (keyword(1).eq.'$dets_disp') then
20     call rdinp(iin)
       if (keyword(1).ne.'$end') then
          nsta_disp=nsta_disp+1
          goto 20
       endif
    endif
    if (.not.lend) goto 15

    ! Quit if the $dets_disp section is missing
    if (nsta_disp.eq.0) then
       errmsg='The $dets_disp section could not be found in '//trim(ain)
       call error_control
    endif

!----------------------------------------------------------------------
! Numbers of adiabatic and diabatic states
!----------------------------------------------------------------------
! For a Tamura-style diabatisation, the number of diabats can be
! greater than the number of adiabats
!----------------------------------------------------------------------
    rewind(iin)

    ! Check to see if the Tamura-style algorithm is being used and,
    ! if so, read the number of diabats
25  call rdinp(iin)
    if (keyword(1).eq.'$algorithm') then
       if (keyword(2).eq.'=') then
          read(keyword(3),'(a)') algorithm
          if (keyword(3).eq.'tamura') then
             if (keyword(4).eq.',') then
                read(keyword(5),*) nsta_diab
             else
                errmsg='The number of diabats has to be given with &
                     the tamura keyword'
                call error_control
             endif
          endif
       else
          errmsg='No argument given with the $algorithm keyword'
          call error_control
       endif
    endif

    ! Set the number of adiabats
    if (algorithm.eq.'tamura') then
       ! Tamura's projection diabatisation
       nsta_adiab=nsta_disp
       nsta=nsta_diab
    else
       ! Pacher, Cederbaum, Koppel propagative block diagonalisation
       ! diabatisation
       nsta_adiab=nsta_ref
       nsta_diab=nsta_adiab
       nsta=nsta_diab
    endif
    
    return
    
  end subroutine get_nsta
    
!######################################################################

  subroutine rdinpfile

    use constants
    use channels
    use iomod
    use parsemod
    use bdglobal
    
    implicit none

    integer :: i,k,l,k1,k2
    
!----------------------------------------------------------------------
! Set defaults
!----------------------------------------------------------------------
    amosref=''
    amosdisp=''
    adetref=''
    adetdisp=''
    lreftrans=.false.
    lrdreftrans=.false.
    ldiabpot=.false.
    ldipole=.false.
    ioverlap=1
    dthresh=1e-6_dp
    normcut=1.0d0
    ltruncate=.false.
    areftrans=''
    avmat=''
    ldmat=.false.
    algorithm='pacher'
    
!----------------------------------------------------------------------
! Second pass: read the input file
!----------------------------------------------------------------------
    rewind(iin)
    
15  continue
    call rdinp(iin,.false.)

    i=0
    if (.not.lend) then

20     continue
       i=i+1

       if (keyword(i).eq.'$mos_ref') then
          if (keyword(i+1).eq.'=') then
             i=i+2
             amosref=keyword(i)
          else
             goto 100
          endif

       else if (keyword(i).eq.'$mos_disp') then
          if (keyword(i+1).eq.'=') then
             i=i+2
             amosdisp=keyword(i)
          else
             goto 100
          endif

       else if (keyword(i).eq.'$dets_ref') then
          do 
             call rdinp(iin,.false.)
             if (keyword(1).eq.'$end') exit
             if (lend) then
                errmsg='End of file reached whilst reading the &
                     $dets_ref section'
                call error_control
             endif
             read(keyword(2),*) k
             adetref(k)=keyword(1)
          enddo
          
       else if (keyword(i).eq.'$dets_disp') then
          do 
             call rdinp(iin,.false.)
             if (keyword(1).eq.'$end') exit
             if (lend) then
                errmsg='End of file reached whilst reading the &
                     $dets_disp section'
                call error_control
             endif
             read(keyword(2),*) k
             adetdisp(k)=keyword(1)
          enddo

       else if (keyword(i).eq.'$energies') then
          ldiabpot=.true.

          if (keyword(i+1).eq.'=') then
             i=i+2
             avmat=keyword(i)
          else
             do
                call rdinp(iin,.false.)
                if (keyword(1).eq.'$end') exit
                if (lend) then
                   errmsg='End of file reached whilst reading the &
                        $energies section'
                   call error_control
                endif
                read(keyword(2),*) k
                read(keyword(1),*) Vmat1(k,k)
             enddo
          endif
             
       else if (keyword(i).eq.'$ref_trans') then
          lreftrans=.true.

          if (keyword(i+1).eq.'=') then
             lrdreftrans=.true.
             i=i+2
             areftrans=keyword(i)
          else
             do
                call rdinp(iin,.false.)
                if (keyword(1).eq.'$end') exit
                if (lend) then
                   errmsg='End of file reached whilst reading the &
                        $ref_trans section'
                   call error_control
                endif
                read(keyword(1),*) k
                read(keyword(2),*) l
                read(keyword(3),*) reftrans(k,l)
             enddo
          endif

       else if (keyword(i).eq.'$overlap') then
          if (keyword(i+1).eq.'=') then
             i=i+2
             if (keyword(i).eq.'fast') then
                ioverlap=1
             else if (keyword(i).eq.'slow') then
                ioverlap=2
             else
                errmsg='Unknown wavefunction overlap algorithm: '&
                     //trim(keyword(i))
                call error_control
             endif
          else
             goto 100 
          endif

       else if (keyword(i).eq.'$dthresh') then
          if (keyword(i+1).eq.'=') then
             i=i+2
             read(keyword(i),*) dthresh
          else
             goto 100 
          endif

       else if (keyword(i).eq.'$norm_cutoff') then
          if (keyword(i+1).eq.'=') then
             ltruncate=.true.
             i=i+2
             read(keyword(i),*) normcut
          else
             goto 100 
          endif
          
       else if (keyword(i).eq.'$dipole') then
          ldipole=.true.
          do
             call rdinp(iin,.false.)
             if (keyword(1).eq.'$end') exit
             if (lend) then
                errmsg='End of file reached whilst reading the &
                     $dipole section'
                call error_control
             endif
             read(keyword(4),*) k1
             read(keyword(5),*) k2
             read(keyword(1),*) adip1(k1,k2,1)
             read(keyword(2),*) adip1(k1,k2,2)
             read(keyword(3),*) adip1(k1,k2,3)
             adip1(k2,k1,:)=adip1(k1,k2,:)
          enddo

       else if (keyword(i).eq.'$dmat_trans') then
          ldmat=.true.

       else if (keyword(i).eq.'$algorithm') then
          if (keyword(i+1).eq.'=') then
             ltruncate=.true.
             i=i+2
             read(keyword(i),'(a)') algorithm
             ! Note that for algorithm='tamura', the no. diabats is
             ! also given here, but that this is now read in get_nsta
             ! but that we need to get to the end of the line here, so
             ! we will just read it in again
             if (algorithm.eq.'tamura') then
                if (keyword(i+1).eq.',') then
                   i=i+2
                   read(keyword(i),*) nsta_diab
                else
                   errmsg='The number of diabats has to be given with &
                        the tamura keyword'
                   call error_control
                endif
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
! Check that all the required information has been given
!----------------------------------------------------------------------
    ! Reference geometry MOs file
    if (amosref.eq.'') then
       errmsg='The name of the reference geometry MO file has not &
            been given'
       call error_control
    endif

    ! Displaced geometry MOs file
    if (amosdisp.eq.'') then
       errmsg='The name of the displaced geometry MO file has not &
            been given'
       call error_control
    endif

    ! Reference geometry determinant file
    do i=1,nsta_ref
       if (adetref(i).eq.'') then
          write(errmsg,'(a,x,i2)') 'The reference geometry &
               determinant file is missing for state',i
          call error_control
       endif
    enddo

    ! Displaced geometry determinant file
    do i=1,nsta_disp
       if (adetdisp(i).eq.'') then
          write(errmsg,'(a,x,i2)') 'The displaced geometry &
               determinant file is missing for state',i
          call error_control
       endif
    enddo

    ! Wavefunction truncation
    if (ltruncate) then
       if (normcut.lt.0.0d0.or.normcut.gt.1.0d0) then
          write(errmsg,'(a,x,ES11.4)') &
               'Nonsensical wavefunction cutoff value:',normcut
          call error_control
       endif
    endif

    ! Diabatisation algorithm: check that the given character string
    ! label is a valid choice
    if (algorithm.ne.'pacher'.and.algorithm.ne.'tamura') then
       errmsg='Unknown diabatisation algorithm: '//trim(algorithm)
       call error_control
    endif
    
    return
    
  end subroutine rdinpfile

!######################################################################

  subroutine initialise

    use constants
    use bdglobal
    use accuracy

    implicit none

!-----------------------------------------------------------------------
! Initialise multigrid data types
!-----------------------------------------------------------------------
    call accuracyInitialize

!-----------------------------------------------------------------------
! Allocate and initialise arrays
!-----------------------------------------------------------------------
    ! Names of the determinant files
    allocate(adetref(nsta_ref))
    allocate(adetdisp(nsta_disp))
    adetref=''
    adetdisp=''
    
    ! Number of determinants for each state at the reference geometry
    allocate(ndet_ref(nsta_ref))
    ndet_ref=0

    ! Number of determinants for each state at the displaced geometry
    allocate(ndet_disp(nsta_disp))
    ndet_disp=0

    ! Reference geometry transformation matrix
    allocate(reftrans(nsta_adiab,nsta_diab))
    reftrans=0.0d0
    
    ! Norms of the wavefunctions at the reference geometry
    allocate(norm_ref(nsta_ref))
    norm_ref=0.0d0

    ! Norms of the wavefunctions at the displaced geometry
    allocate(norm_disp(nsta_disp))
    norm_disp=0.0d0
    
    ! Overlaps between the electronic states at the two geometries
    if (algorithm.eq.'tamura') then
       allocate(spsi(nsta_disp,nsta_ref))
    else
       allocate(spsi(nsta_adiab,nsta_adiab))
    endif
    spsi=0.0d0

    ! Adiabatic potential matrix
    allocate(Vmat(nsta_adiab,nsta_adiab))
    Vmat=0.0d0
    
    ! Temporary adiabatic potential matrix
    allocate(Vmat1(nsta_disp,nsta_disp))
    Vmat1=0.0d0

    ! Diabatic potential matrix
    allocate(Wmat(nsta_diab,nsta_diab))
    Wmat=0.0d0

    ! ADT matrix
    allocate(adt(nsta_adiab,nsta_diab))
    adt=0.0d0

    ! Adiabatic dipole matrix
    allocate(adip(nsta_adiab,nsta_adiab,3))
    adip=0.0d0

    ! Temporary adiabatic dipole matrix
    allocate(adip1(nsta_disp,nsta_disp,3))
    adip1=0.0d0

    ! Diabatic dipole matrix
    allocate(ddip(nsta_diab,nsta_diab,3))
    ddip=0.0d0
    
    return
    
  end subroutine initialise

!######################################################################

  subroutine finalise

    use constants
    use bdglobal
    use accuracy

    implicit none

    deallocate(adetref)
    deallocate(adetdisp)
    deallocate(ndet_ref)
    deallocate(ndet_disp)
    deallocate(reftrans)
    deallocate(norm_ref)
    deallocate(norm_disp)
    deallocate(spsi)
    deallocate(Vmat)
    deallocate(Vmat1)
    deallocate(Wmat)
    deallocate(adt)
    deallocate(adip)
    deallocate(adip1)
    deallocate(ddip)
    deallocate(c_ref)
    deallocate(c_disp)
    if (allocated(det_ref)) deallocate(det_ref)
    if (allocated(det_disp)) deallocate(det_disp)
    deallocate(iocca_ref)
    deallocate(ioccb_ref)
    deallocate(iocca_disp)
    deallocate(ioccb_disp)
    if (allocated(isel)) deallocate(isel)
    
    return
    
  end subroutine finalise

!######################################################################

  subroutine rdvmat

    use constants
    use channels
    use ioqc
    use bdglobal
    
    implicit none

    integer                        :: i
    real(dp), dimension(nsta_disp) :: vtmp

    call geten(vtmp,nsta_disp,avmat)
    
    do i=1,nsta_disp
       Vmat1(i,i)=vtmp(i)
    enddo

    return
    
  end subroutine rdvmat
  
!######################################################################

  subroutine loadmos

    use constants
    use channels
    use bdglobal
    use import_gamess

    implicit none

    integer :: i
    
!-----------------------------------------------------------------------
! (1) Reference geometry MOs
!-----------------------------------------------------------------------
    write(ilog,'(/,82a)') ('+',i=1,82)
    write(ilog,'(2x,a)') 'Loading Reference Geometry MOs'
    write(ilog,'(82a)') ('+',i=1,82)
    call gamess_load_orbitals(file=amosref,structure=gam_ref)
    
!-----------------------------------------------------------------------
! (2) Displaced geometry MOs
!-----------------------------------------------------------------------
    write(ilog,'(/,82a)') ('+',i=1,82)
    write(ilog,'(2x,a)') 'Loading Displaced Geometry MOs'
    write(ilog,'(82a)') ('+',i=1,82)
    call gamess_load_orbitals(file=amosdisp,structure=gam_disp)
    
    return
    
  end subroutine loadmos

!######################################################################

  subroutine wrgeoms

    use constants
    use channels
    use bdglobal

    integer :: i,j
    
    write(ilog,'(/,82a)') ('+',i=1,82)
    write(ilog,'(2x,a)') 'Geometries'
    write(ilog,'(82a)') ('+',i=1,82)

    write(ilog,'(/,a)') 'Reference geometry:'
    do i=1,gam_ref%natoms
       write(ilog,'(1x,a2,3(2x,F10.7))') &
            gam_ref%atoms(i)%name,(gam_ref%atoms(i)%xyz(j),j=1,3)
    enddo
    
    write(ilog,'(/,a)') 'Displaced geometry:'
    do i=1,gam_disp%natoms
       write(ilog,'(1x,a2,3(2x,F10.7))') &
            gam_disp%atoms(i)%name,(gam_disp%atoms(i)%xyz(j),j=1,3)
    enddo
    
    return

  end subroutine wrgeoms
    
!######################################################################

  subroutine wrnorms

    use constants
    use channels
    use bdglobal
    
    implicit none

    integer :: i
    
!----------------------------------------------------------------------
! Section header
!----------------------------------------------------------------------
    write(ilog,'(/,82a)') ('+',i=1,82)
    write(ilog,'(2x,a)') 'Wavefunction Norms'
    write(ilog,'(82a)') ('+',i=1,82)

!----------------------------------------------------------------------
! Reference geometry wavefunctions
!----------------------------------------------------------------------
    write(ilog,'(/,47a)') ('-',i=1,47)
    write(ilog,'(a)') '   Ref. State  |   Norm'
    write(ilog,'(a)') '     |I>       | || |I> ||'
    write(ilog,'(47a)') ('-',i=1,47)
    do i=1,nsta_ref
       write(ilog,'(5x,i2,11x,F13.10)') i,norm_ref(i)
    enddo
    
!----------------------------------------------------------------------
! Displaced geometry wavefunctions
!----------------------------------------------------------------------
    write(ilog,'(/,47a)') ('-',i=1,47)
    write(ilog,'(a)') '  Disp. State  |   Norm'
    write(ilog,'(a)') '     |I>       | || |I> ||'
    write(ilog,'(47a)') ('-',i=1,47)
    do i=1,nsta_disp
       write(ilog,'(5x,i2,11x,F13.10)') i,norm_disp(i)
    enddo

    return
    
  end subroutine wrnorms  

!######################################################################

  subroutine phase_refpsi

    use constants
    use channels
    use iomod
    use bdglobal
    
    implicit none

    integer                      :: unit,i
    integer, dimension(nsta_ref) :: refphase
    character(len=120)           :: string
    logical                      :: found

!-----------------------------------------------------------------------
! Exit if the previous log file does not exist
!-----------------------------------------------------------------------
    inquire(file=trim(areftrans),exist=found)
    
    if (.not.found) then
       write(6,'(/,2x,a,/)') 'The file '//trim(areftrans)&
            //' does not exist'
       stop
    endif

!-----------------------------------------------------------------------
! Open the previous log file
!-----------------------------------------------------------------------
    call freeunit(unit)
    open(unit,file=areftrans,form='formatted',status='old')

!-----------------------------------------------------------------------
! Read the phase factors from the old log file
!-----------------------------------------------------------------------
    refphase=1
    
5   read(unit,'(a)',end=999) string
    if (index(string,'Phase Factors').eq.0) goto 5
    
    read(unit,*)

    do i=1,nsta_ref
       read(unit,'(a)') string
       if (string.eq.'') exit
       read(string,'(6x,i2)') refphase(i)
    enddo
    
!-----------------------------------------------------------------------
! Close the previous log file
!-----------------------------------------------------------------------
    close(unit)

!-----------------------------------------------------------------------
! Adjust the phases of the ref. geometry wavefunctions
!-----------------------------------------------------------------------
    do i=1,nsta_ref
       c_ref(:,i)=c_ref(:,i)*refphase(i)
    enddo
    
    return

999 continue
    errmsg='The Phase Factors section could not be found in: '&
         //trim(areftrans)
    call error_control
    
  end subroutine phase_refpsi
    
!######################################################################
  
end program blockdiag


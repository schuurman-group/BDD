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
  use mooverlaps
  use wfoverlaps
  use adtmod
  use bdglobal
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
! Determine the disp. geometry states to use by trying to follow
! electronic characters from geometry to geometry
!----------------------------------------------------------------------
  if (nsta_disp.ne.nsta_ref) then
     call trackwfs
  else
     Vmat=Vmat1
  endif
     
!----------------------------------------------------------------------
! Calculate the overlaps between the electronic states at the
! reference and displaced geometries
!----------------------------------------------------------------------
  call get_overlaps
  
!----------------------------------------------------------------------
! Write the wavefunction overlaps to the log file
!----------------------------------------------------------------------
  call wroverlaps
  
!----------------------------------------------------------------------
! New rephasing algorithm
!----------------------------------------------------------------------
  call rephase

!----------------------------------------------------------------------
! Optional transformation of the reference geometry wavefunctions
! (via the transformation of the wavefunction overlap matrix)
!
! Maybe this should be moved to psi_overlaps?
!----------------------------------------------------------------------
  if (lreftrans) call trans_refpsi
  
!----------------------------------------------------------------------
! Calculate the ADT matrix
!----------------------------------------------------------------------
  call get_adt

!----------------------------------------------------------------------
! Check on the swapping of diabats. Note that this can only be done if
! we have access to the ADT matrix of the previous geometry
!----------------------------------------------------------------------
  if (lreftrans) call switch_diabats
  
!----------------------------------------------------------------------
! Write the ADT matrix to the log file
!----------------------------------------------------------------------
  call write_adt
  
!----------------------------------------------------------------------
! Optional: calculate and output the quasi-diabatic potential matrix
!----------------------------------------------------------------------
  if (ldiabpot) call diabpotmat

!----------------------------------------------------------------------
! Optional: calculate and output the quasi-diabatic dipole matrix
!----------------------------------------------------------------------
  if (ldipole) call diabdipmat
  
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

    ! nsta (the dimension of the diabatic potential matrix) is going
    ! to be equal to nsta_ref no matter how we are operating
    nsta=nsta_ref

!----------------------------------------------------------------------
! Determine the number of disp. geometry states from the input file.
! Note that nsta_disp can be different to nsta_ref if we trying to
! follow states by character from geometry to geometry.
!----------------------------------------------------------------------
    rewind(iin)

    ! Determine the no. states from the $dets_disp section
    nsta_disp=0
15   call rdinp(iin)
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
    allocate(reftrans(nsta,nsta))
    reftrans=0.0d0
    
    ! Norms of the wavefunctions at the reference geometry
    allocate(norm_ref(nsta_ref))
    norm_ref=0.0d0

    ! Norms of the wavefunctions at the displaced geometry
    allocate(norm_disp(nsta_disp))
    norm_disp=0.0d0
    
    ! Overlaps between the electronic states at the two geometries
    allocate(spsi(nsta,nsta))
    spsi=0.0d0

    ! Adiabatic potential matrix
    allocate(Vmat(nsta,nsta))
    Vmat=0.0d0
    
    ! Temporary adiabatic potential matrix
    allocate(Vmat1(nsta_disp,nsta_disp))
    Vmat1=0.0d0

    ! Diabatic potential matrix
    allocate(Wmat(nsta,nsta))
    Wmat=0.0d0

    ! ADT matrix
    allocate(adt(nsta,nsta))
    adt=0.0d0

    ! Adiabatic dipole matrix
    allocate(adip(nsta,nsta,3))
    adip=0.0d0

    ! Temporary adiabatic dipole matrix
    allocate(adip1(nsta_disp,nsta_disp,3))
    adip1=0.0d0

    ! Diabatic dipole matrix
    allocate(ddip(nsta,nsta,3))
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
    
    return
    
  end subroutine finalise
    
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
  
  subroutine rddetfiles

    use constants
    use iomod
    use bdglobal
    use timingmod
    
    implicit none

    integer  :: i,k,n,idet,ilbl
    real(dp) :: tw1,tw2,tc1,tc2
    
!----------------------------------------------------------------------
! Output what we are doing
!----------------------------------------------------------------------
    write(ilog,'(/,82a)') ('+',i=1,82)
    write(ilog,'(2x,a)') 'Parsing the determinant files'
    write(ilog,'(82a)') ('+',i=1,82)
    
!-----------------------------------------------------------------------
! Start timing
!-----------------------------------------------------------------------
    call times(tw1,tc1)

!-----------------------------------------------------------------------
! Try to determine whether we have binary or ascii determinant files.
! We are here assuming that binary files have a .bin file extension, as
! is the case for binary determinant files written by the mrci code.
!-----------------------------------------------------------------------
    ilbl=len_trim(adetref(1))
    if (adetref(1)(ilbl-3:ilbl).eq.'.bin') then
       lbinary=.true.
    else
       lbinary=.false.
    endif

!-----------------------------------------------------------------------
! Call to the appropriate parsing routine
!-----------------------------------------------------------------------
    if (lbinary) then
       call rddetfiles_binary
    else
       call rddetfiles_ascii
    endif
    
!-----------------------------------------------------------------------
! Norms
!-----------------------------------------------------------------------
    ! Reference geometry
    do i=1,nsta_ref
       norm_ref(i)=sqrt(sum(c_ref(:,i)**2))
    enddo

    ! Displaced geometry
    do i=1,nsta_disp
       norm_disp(i)=sqrt(sum(c_disp(:,i)**2))
    enddo
    
!-----------------------------------------------------------------------
! Normalisation of the wavefunctions
!-----------------------------------------------------------------------
    do i=1,nsta_ref
       c_ref(:,i)=c_ref(:,i)/norm_ref(i)
    enddo

    do i=1,nsta_disp
       c_disp(:,i)=c_disp(:,i)/norm_disp(i)
    enddo
    
!-----------------------------------------------------------------------    
! Output timings
!-----------------------------------------------------------------------    
    call times(tw2,tc2)
    write(ilog,'(/,2x,a,1x,F9.2,1x,a)') &
         'Wall Time For Determinant Parsing:',tw2-tw1," s"
    write(ilog,'(2x,a,2x,F9.2,1x,a)') &
         'CPU Time For Determinant Parsing:',tc2-tc1," s"

    return
    
  end subroutine rddetfiles

!######################################################################

  subroutine rddetfiles_ascii

    use constants
    use iomod
    use parsemod
    use bdglobal
    use utils
    use timingmod
    
    implicit none

    integer  :: i,k,n,idet

!-----------------------------------------------------------------------
! First pass: determine the no. determinants for each file
!-----------------------------------------------------------------------
    do i=1,nsta_ref
       ndet_ref(i)=nlines(adetref(i))
    enddo

    do i=1,nsta_disp
       ndet_disp(i)=nlines(adetdisp(i))
    enddo
    
    ! Maximum number of determinants
    maxdet=max(maxval(ndet_ref),maxval(ndet_disp))

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
    ! Number of MOs for the reference and displaced geometries
    nmo_ref=gam_ref%nvectors
    nmo_disp=gam_disp%nvectors

    ! Coefficient vectors
    allocate(c_ref(maxdet,nsta_ref))
    allocate(c_disp(maxdet,nsta_disp))
    c_ref=0.0d0
    c_disp=0.0d0

    ! Determinant vectors
    allocate(det_ref(nmo_ref,maxdet,nsta_ref))
    allocate(det_disp(nmo_disp,maxdet,nsta_disp))
    det_ref=0
    det_disp=0

!-----------------------------------------------------------------------
! Read in the determinants and coefficients
!-----------------------------------------------------------------------
    call freeunit(idet)

    ! Reference geometry
    do i=1,nsta_ref
       open(idet,file=adetref(i),form='formatted',status='old')
       do k=1,ndet_ref(i)
          call rdinp(idet)
          read(keyword(1),*) c_ref(k,i)
          do n=2,inkw
             read(keyword(n),*) det_ref(n-1,k,i)
          enddo          
       enddo
       close(idet)
    enddo
    
    ! Displaced geometry
    do i=1,nsta_disp
       open(idet,file=adetdisp(i),form='formatted',status='old')
       do k=1,ndet_disp(i)
          call rdinp(idet)
          read(keyword(1),*) c_disp(k,i)          
          do n=2,inkw
             read(keyword(n),*) det_disp(n-1,k,i)
          enddo          
       enddo
       close(idet)
    enddo

!----------------------------------------------------------------------
! Get the alpha and beta spinorbital indices for every determinant
!----------------------------------------------------------------------
    call alpha_beta_indices
    
    return
    
  end subroutine rddetfiles_ascii

!######################################################################

  subroutine alpha_beta_indices

    use constants
    use channels
    use iomod
    use bdglobal
    
    implicit none

    integer :: imo,nar,nad,nbr,nbd,na,nb
    integer :: i,k
    
!-----------------------------------------------------------------------
! Determine the no. alpha and beta spinorbitals
!-----------------------------------------------------------------------
    ! Ref. States
    nar=0
    nbr=0
    do imo=1,nmo_ref
       if (det_ref(imo,1,1).eq.2) then
          nar=nar+1
          nbr=nbr+1
       else if (det_ref(imo,1,1).eq.+1) then
          nar=nar+1
       else if (det_ref(imo,1,1).eq.-1) then
          nbr=nbr+1
       endif
    enddo

    ! Disp. states
    nad=0
    nbd=0
    do imo=1,nmo_disp
       if (det_disp(imo,1,1).eq.2) then
          nad=nad+1
          nbd=nbd+1
       else if (det_disp(imo,1,1).eq.+1) then
          nad=nad+1
       else if (det_disp(imo,1,1).eq.-1) then
          nbd=nbd+1
       endif
    enddo

    ! Exit if the numbers of alpha and beta electrons in the ref. and
    ! disp. states is not consistent
    if (nar.ne.nad.or.nbr.ne.nbd) then
       errmsg='Inconsistent numbers of alpha and beta electrons in &
            the ref. and disp. states'
       call error_control
    endif

    ! Set the number of alpha and beta spinorbitals
    nalpha=nar
    nbeta=nbr

!-----------------------------------------------------------------------
! Allocate and initialise the spinorbital index arrays
!-----------------------------------------------------------------------
    allocate(iocca_ref(nalpha,maxdet,nsta_ref))
    allocate(ioccb_ref(nbeta,maxdet,nsta_ref))
    allocate(iocca_disp(nalpha,maxdet,nsta_disp))
    allocate(ioccb_disp(nbeta,maxdet,nsta_disp))
    iocca_ref=0
    ioccb_ref=0
    iocca_disp=0
    ioccb_disp=0

!-----------------------------------------------------------------------
! Fill in the spinorbital index arrays
!-----------------------------------------------------------------------
    ! Ref. states
    !
    ! Loop over states
    do i=1,nsta_ref
       ! Loop over determinants
       do k=1,ndet_ref(i)
          ! Fill in the spinorbital indicies for the current
          ! determinant
          na=0
          nb=0
          do imo=1,nmo_ref
             if (det_ref(imo,k,i).eq.2) then
                na=na+1
                nb=nb+1
                iocca_ref(na,k,i)=imo
                ioccb_ref(nb,k,i)=imo
             else if (det_ref(imo,k,i).eq.+1) then
                na=na+1
                iocca_ref(na,k,i)=imo
             else if (det_ref(imo,k,i).eq.-1) then
                nb=nb+1
                ioccb_ref(nb,k,i)=imo
             endif
          enddo
       enddo
    enddo

    ! Disp. states
    !
    ! Loop over states
    do i=1,nsta_disp
       ! Loop over determinants
       do k=1,ndet_disp(i)
          ! Fill in the spinorbital indicies for the current
          ! determinant
          na=0
          nb=0
          do imo=1,nmo_disp
             if (det_disp(imo,k,i).eq.2) then
                na=na+1
                nb=nb+1
                iocca_disp(na,k,i)=imo
                ioccb_disp(nb,k,i)=imo
             else if (det_disp(imo,k,i).eq.+1) then
                na=na+1
                iocca_disp(na,k,i)=imo
             else if (det_disp(imo,k,i).eq.-1) then
                nb=nb+1
                ioccb_disp(nb,k,i)=imo
             endif
          enddo
       enddo
    enddo
    
    return
    
  end subroutine alpha_beta_indices
  
!######################################################################

  subroutine rddetfiles_binary

    use constants
    use iomod
    use bdglobal
    use timingmod
    
    implicit none

    integer :: i,n,idet,nbeta_last,nalpha_last,itmp

!-----------------------------------------------------------------------
! First pass: read dimensions and allocate arrays
!-----------------------------------------------------------------------
    call freeunit(idet)

    ! Reference geometry states
    do i=1,nsta_ref

       ! Open the next determinant file
       open(idet,file=adetref(i),form='unformatted',status='old')

       ! Read in the no. determinants and the no. alpha and beta
       ! electrons
       read(idet) ndet_ref(i)
       read(idet) nalpha
       read(idet) nbeta

       ! Check that the no. alpha and beta electrons is consistent
       if (i.gt.1) then
          if (nalpha.ne.nalpha_last.or.nbeta.ne.nbeta_last) then
             errmsg='Inconsistent numbers of alpha and beta electrons'
             call error_control
          endif
       endif

       ! Reset nalpha_last and nbeta_last for the next iteration
       nalpha_last=nalpha
       nbeta_last=nbeta

       ! Close the determinant file
       close(idet)

    enddo

    ! Displaced geometry states
    do i=1,nsta_disp

       ! Open the next determinant file
       open(idet,file=adetdisp(i),form='unformatted',status='old')

       ! Read in the no. determinants and the no. alpha and beta
       ! electrons
       read(idet) ndet_disp(i)
       read(idet) nalpha
       read(idet) nbeta

       ! Check that the no. alpha and beta electrons is consistent
       if (nalpha.ne.nalpha_last.or.nbeta.ne.nbeta_last) then
          errmsg='Inconsistent numbers of alpha and beta electrons'
          call error_control
       endif
       
       ! Reset nalpha_last and nbeta_last for the next iteration
       nalpha_last=nalpha
       nbeta_last=nbeta

       ! Close the determinant file
       close(idet)

    enddo

    ! Maximum number of determinants
    maxdet=max(maxval(ndet_ref),maxval(ndet_disp))

    ! Number of MOs for the reference and displaced geometries (needed
    ! elsewhere)
    nmo_ref=gam_ref%nvectors
    nmo_disp=gam_disp%nvectors
    
!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
    ! Indices of the occupied alpha and beta orbitals
    allocate(iocca_ref(nalpha,maxdet,nsta_ref))
    allocate(ioccb_ref(nbeta,maxdet,nsta_ref))
    allocate(iocca_disp(nalpha,maxdet,nsta_disp))
    allocate(ioccb_disp(nbeta,maxdet,nsta_disp))
    iocca_ref=0
    ioccb_ref=0
    iocca_disp=0
    ioccb_disp=0
    
    ! Coefficient vectors
    allocate(c_ref(maxdet,nsta_ref))
    allocate(c_disp(maxdet,nsta_disp))
    c_ref=0.0d0
    c_disp=0.0d0
    
!-----------------------------------------------------------------------
! Read in the determinants (in terms of alpha and beta strings) and
! coefficients
!-----------------------------------------------------------------------
    ! Reference geometry
    do i=1,nsta_ref
       open(idet,file=adetref(i),form='unformatted',status='old')
       read(idet) itmp
       read(idet) itmp
       read(idet) itmp
       do n=1,ndet_ref(i)
          read(idet) iocca_ref(:,n,i)
       enddo
       do n=1,ndet_ref(i)
          read(idet) ioccb_ref(:,n,i)
       enddo
       read(idet) c_ref(1:ndet_ref(i),i)
       close(idet)
    enddo

    ! Displaced geometry
    do i=1,nsta_disp
       open(idet,file=adetdisp(i),form='unformatted',status='old')
       read(idet) itmp
       read(idet) itmp
       read(idet) itmp
       do n=1,ndet_disp(i)
          read(idet) iocca_disp(:,n,i)
       enddo
       do n=1,ndet_disp(i)
          read(idet) ioccb_disp(:,n,i)
       enddo
       read(idet) c_disp(1:ndet_disp(i),i)
       close(idet)
    enddo

    return
    
  end subroutine rddetfiles_binary
    
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

  subroutine trackwfs

    use constants
    use channels
    use wfoverlaps
    use bdglobal
    use utils
    use iomod
    
    implicit none

    integer                         :: i,j,k,maxdet1,nok,n
    integer, allocatable            :: indx(:),ndet1_disp(:),&
                                       ndet1_ref(:)
    integer, allocatable            :: iocca1_disp(:,:,:),&
                                       ioccb1_disp(:,:,:),&
                                       iocca1_ref(:,:,:),&
                                       ioccb1_ref(:,:,:)
    integer, allocatable            :: iswapvec1(:)
    integer, allocatable            :: iswapvec3(:,:,:)
    integer, parameter              :: nsmall=200
    real(dp), allocatable           :: fswapvec2(:,:)
    real(dp), allocatable           :: spsi1(:,:)
    real(dp), allocatable           :: c1_disp(:,:),c1_ref(:,:)
    real(dp), allocatable           :: cabs(:)
    real(dp), parameter             :: thrsh=0.8d0
    real(dp)                        :: norm
    character(len=120), allocatable :: aswapvec1(:)
    
!----------------------------------------------------------------------
! Output what we are doing
!----------------------------------------------------------------------
    write(ilog,'(/,82a)') ('+',i=1,82)
    write(ilog,'(2x,a)') 'Tracking wavefunctions by overlap'
    write(ilog,'(82a)') ('+',i=1,82)
    
!----------------------------------------------------------------------
! Allocate arrays.
!    
! To speed things up, we will only use a 'handful' of determinants for
! each state in the calculation of the overlaps. This should be
! sufficient for the following of electronic state character from
! geometry to geometry.
!----------------------------------------------------------------------
    ! Wavefunction overlaps
    allocate(spsi1(nsta_disp,nsta_ref))
    spsi1=0.0d0
    
    ! Reduced number of determinants for the disp. states
    allocate(ndet1_disp(nsta_disp))
    do i=1,nsta_disp
       ndet1_disp(i)=min(nsmall,ndet_disp(i))
    enddo

    ! Reduced number of determinants for the ref. states
    allocate(ndet1_ref(nsta_ref))
    do i=1,nsta_ref
       ndet1_ref(i)=min(nsmall,ndet_ref(i))
    enddo

    ! Maximum reduced number of determinants
    maxdet1=max(maxval(ndet1_ref),maxval(ndet1_disp))
    
    ! Determinant coefficient and spinorbital index arrays
    allocate(c1_disp(maxdet1,nsta_disp))
    allocate(c1_ref(maxdet1,nsta_ref))
    allocate(iocca1_disp(nalpha,maxdet1,nsta_disp))
    allocate(iocca1_ref(nalpha,maxdet1,nsta_ref))
    allocate(ioccb1_disp(nbeta,maxdet1,nsta_disp))
    allocate(ioccb1_ref(nbeta,maxdet1,nsta_ref))
    c1_disp=0.0d0
    c1_ref=0.0d0
    iocca1_disp=0
    iocca1_ref=0
    ioccb1_disp=0
    ioccb1_ref=0

!----------------------------------------------------------------------
! Fill in the reduced determinant coefficient and spinorbital index
! arrays
!----------------------------------------------------------------------
    ! Allocate arrays
    allocate(indx(maxdet))
    allocate(cabs(maxdet))

    ! Disp. geometry arrays
    do i=1,nsta_disp

       indx=0
       cabs(:)=abs(c_disp(:,i))

       call dsortindxa1('D',ndet_disp(i),cabs(1:ndet_disp(i)),&
            indx(1:ndet_disp(i)))

       do k=1,ndet1_disp(i)
          c1_disp(k,i)=c_disp(indx(k),i)
          iocca1_disp(:,k,i)=iocca_disp(:,indx(k),i)
          ioccb1_disp(:,k,i)=ioccb_disp(:,indx(k),i)
       enddo
       
    enddo
    
    ! Ref. geometry arrays
    do i=1,nsta_ref

       indx=0
       cabs(:)=abs(c_ref(:,i))

       call dsortindxa1('D',ndet_ref(i),cabs(1:ndet_ref(i)),&
            indx(1:ndet_ref(i)))

       do k=1,ndet1_ref(i)
          c1_ref(k,i)=c_ref(indx(k),i)
          iocca1_ref(:,k,i)=iocca_ref(:,indx(k),i)
          ioccb1_ref(:,k,i)=ioccb_ref(:,indx(k),i)
       enddo
       
    enddo
    
    ! Deallocate arrays
    deallocate(indx)
    deallocate(cabs)

!----------------------------------------------------------------------
! Normalisation of the truncated wavefunctions
!----------------------------------------------------------------------
    ! Reference geometry
    do i=1,nsta_ref
       norm=sqrt(sum(c1_ref(:,i)**2))
       c1_ref(:,i)=c1_ref(:,i)/norm
    enddo

    ! Displaced geometry
    do i=1,nsta_disp
       norm=sqrt(sum(c1_disp(:,i)**2))
       c1_disp(:,i)=c1_disp(:,i)/norm
    enddo
    
!----------------------------------------------------------------------
! Calculate the overlaps between the (larger) set of electronic states
! at the displaced geometry and the (smaller) set of electronic states
! at the reference geometry.
!----------------------------------------------------------------------
    call psi_overlaps(spsi1,nsta_disp,nsta_ref,nalpha,nbeta,&
         ndet1_disp,ndet1_ref,nmo_disp,nmo_ref,maxdet1,c1_disp,c1_ref,&
         iocca1_disp,iocca1_ref,ioccb1_disp,ioccb1_ref,ioverlap,smo,&
         dthresh)

!----------------------------------------------------------------------
! Write the overlaps to the log file
!----------------------------------------------------------------------
    write(ilog,'(/,47a)') ('-',i=1,47)
    write(ilog,'(a)') '  Disp. State  |  Ref. State  |  Overlap'
    write(ilog,'(a)') '     |I>       |    |J>       |   <I|J>'
    write(ilog,'(47a)') ('-',i=1,47)
    do i=1,nsta_disp
       do j=1,nsta_ref
          write(ilog,'(5x,i2,13x,i2,11x,F13.10)') i,j,spsi1(i,j)
       enddo
    enddo
    write(ilog,'(47a)') ('-',i=1,47)

!----------------------------------------------------------------------
! Determine which disp. states correspond to the ref. states.
!----------------------------------------------------------------------
    allocate(indx(nsta_ref))
    indx=0

    ! Good overlap counter
    nok=0

    ! Loop over ref. states
    do i=1,nsta_ref
       ! Loop over disp. states
       do j=1,nsta_disp
          if (abs(spsi1(j,i)).gt.thrsh) then
             nok=nok+1
             indx(i)=j
          endif
       enddo
    enddo

    ! Exit if any ref. states do not correspond to a disp. state
    do i=1,nsta_ref
       n=0
       do j=1,nsta_disp
          if (indx(i).eq.j) n=n+1
       enddo
       if (n.ne.1) then
          errmsg='Not all ref. states correspond to a single disp. &
               state. Quitting.'
          call error_control
       endif
    enddo

!----------------------------------------------------------------------
! Write the selected state indices to the log file
!----------------------------------------------------------------------
    write(ilog,'(a)') ''
    do i=1,nsta_ref
       write(ilog,'(x,2(x,a,x,i0))') 'Selected state',i,':',indx(i)
    enddo
    
!----------------------------------------------------------------------
! Re-fill the disp. geometry arrays
!----------------------------------------------------------------------
    ! Number of determinants for the disp. states
    allocate(iswapvec1(nsta_ref))
    iswapvec1=0
    do i=1,nsta_ref
       iswapvec1(i)=ndet_disp(indx(i))
    enddo
    ndet_disp=0
    ndet_disp(1:nsta_ref)=iswapvec1
    deallocate(iswapvec1)

    ! Disp. state coefficient vectors
    allocate(fswapvec2(maxdet,nsta_ref))
    do i=1,nsta_ref
       fswapvec2(:,i)=c_disp(:,indx(i))
    enddo
    c_disp=0.0d0
    c_disp(:,1:nsta_ref)=fswapvec2
    deallocate(fswapvec2)

    ! Disp. state alpha spinorbital index arrays
    allocate(iswapvec3(nalpha,maxdet,nsta_ref))
    do i=1,nsta_ref
       iswapvec3(:,:,i)=iocca_disp(:,:,indx(i))
    enddo
    iocca_disp=0
    iocca_disp(:,:,1:nsta_ref)=iswapvec3
    deallocate(iswapvec3)

    ! Disp. state beta spinorbital index arrays
    allocate(iswapvec3(nbeta,maxdet,nsta_ref))
    do i=1,nsta_ref
       iswapvec3(:,:,i)=ioccb_disp(:,:,indx(i))
    enddo
    ioccb_disp=0
    ioccb_disp(:,:,1:nsta_ref)=iswapvec3
    deallocate(iswapvec3)

    ! Adiabatic potential array at the disp. geometry
    Vmat=0.0d0
    do i=1,nsta_ref
       Vmat(i,i)=Vmat1(indx(i),indx(i))
    enddo

    ! Adiabatic dipole matrix
    adip=0.0d0
    do i=1,nsta_ref
       do j=1,nsta_ref
          adip(i,j,:)=adip1(indx(i),indx(j),:)
       enddo
    enddo
    
    ! Names of the disp. determinant files: this has to be reset
    ! otherwise the rephasing routine will over-write files
    ! incorrectly
    allocate(aswapvec1(nsta_ref))
    do i=1,nsta_ref
       aswapvec1(i)=adetdisp(indx(i))
    enddo
    adetdisp=''
    adetdisp(1:nsta_ref)=aswapvec1
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(spsi1)
    deallocate(ndet1_disp)
    deallocate(ndet1_ref)
    deallocate(c1_disp)
    deallocate(c1_ref)
    deallocate(iocca1_disp)
    deallocate(iocca1_ref)
    deallocate(ioccb1_disp)
    deallocate(ioccb1_ref)
    deallocate(indx)
    
    return
    
  end subroutine trackwfs

!######################################################################

  subroutine get_overlaps

    use constants
    use channels
    use wfoverlaps
    use bdglobal
    use iomod
    
    implicit none

    integer :: i
    
!----------------------------------------------------------------------
! Output what we are doing
!----------------------------------------------------------------------
    write(ilog,'(/,82a)') ('+',i=1,82)
    write(ilog,'(2x,a)') 'Calculating Adiabatic Wavefunction Overlaps'
    write(ilog,'(82a)') ('+',i=1,82)
    
!----------------------------------------------------------------------
! Calculate the wavefunction overlaps
!----------------------------------------------------------------------
    call psi_overlaps(spsi,nsta,nsta,nalpha,nbeta,ndet_disp(1:nsta),&
         ndet_ref(1:nsta),nmo_disp,nmo_ref,maxdet,c_disp(:,1:nsta),&
         c_ref(:,1:nsta),iocca_disp(:,:,1:nsta),iocca_ref(:,:,1:nsta),&
         ioccb_disp(:,:,1:nsta),ioccb_ref(:,:,1:nsta),ioverlap,smo,&
         dthresh)
    
    return
    
  end subroutine get_overlaps
    
!######################################################################

  subroutine wroverlaps

    use constants
    use channels
    use bdglobal
    
    implicit none

    integer              :: i,j
    real(dp), parameter  :: ovrthrsh=0.5d0
    logical              :: lovrlp

!----------------------------------------------------------------------
! Print out the wavefunction overlaps
!----------------------------------------------------------------------
    ! Table header
    write(ilog,'(/,47a)') ('-',i=1,47)
    write(ilog,'(a)') '  Disp. State  |  Ref. State  |  Overlap'
    write(ilog,'(a)') '     |I>       |    |J>       |   <I|J>'
    write(ilog,'(47a)') ('-',i=1,47)

    ! Table entries
    do i=1,nsta
       do j=1,nsta
          write(ilog,'(5x,i2,13x,i2,11x,F13.10)') i,j,spsi(i,j)
       enddo
    enddo
       
    ! End of the table
    write(ilog,'(47a)') ('-',i=1,47)

!----------------------------------------------------------------------
! Print out a warning if it looks like a state from outside the group
! of interest has crossed in
!----------------------------------------------------------------------
    ! Loop over disp. states
    do i=1,nsta

       lovrlp=.false.
       
       ! Loop over ref states
       do j=1,nsta
          if (abs(spsi(i,j)).ge.ovrthrsh) lovrlp=.true.
       enddo
       
       if (.not.lovrlp) write(ilog,'(/,2x,a)') &
            'WARNING: state crossing detected!'

    enddo
    
    return
    
  end subroutine wroverlaps
    
!######################################################################

  subroutine rephase

    use constants
    use iomod
    use bdglobal
    
    implicit none

    integer                   :: i,j,k,n,idet
    real(dp), dimension(nsta) :: phfac
    character(len=128)        :: fmt

!-----------------------------------------------------------------------
! Try to determine if the phase of a wavefunction has switched from the
! ref. geometry
!-----------------------------------------------------------------------    
    phfac=1.0d0
    do i=1,nsta
       do j=1,nsta
          if (abs(spsi(i,j)).gt.0.8d0) then
             if (spsi(i,j).lt.0.0d0) phfac(i)=-1.0d0
          endif
       enddo
    enddo

!-----------------------------------------------------------------------
! Re-phase the overlap matrix elements.
! Note that we are here assuming that the ref. states have the
! 'correct' phase.
!-----------------------------------------------------------------------
    do i=1,nsta
       spsi(i,:)=phfac(i)*spsi(i,:)
    enddo

!-----------------------------------------------------------------------
! Re-write the disp. state determinant files with the 'correct', i.e.,
! consistent, phases
!-----------------------------------------------------------------------
    write(fmt,'(a,i4,a)')'(ES15.8,',nmo_disp,'i5)'

    call freeunit(idet)

    if (lbinary) then
       ! Binary determinant files
       do i=1,nsta
          if (phfac(i).eq.1.0d0) cycle
          open(idet,file=adetdisp(i),form='unformatted',status='old')
          write(idet) ndet_disp(i)
          write(idet) nalpha
          write(idet) nbeta
          do n=1,ndet_disp(i)
             write(idet) iocca_disp(:,n,i)
          enddo
          do n=1,ndet_disp(i)
             write(idet) ioccb_disp(:,n,i)
          enddo
          write(idet) -c_disp(1:ndet_disp(i),i)
          close(idet)
       enddo
    else
       ! Ascii determinant files
       do i=1,nsta
          if (phfac(i).eq.1.0d0) cycle
          open(idet,file=adetdisp(i),form='formatted',status='old')
          do k=1,ndet_disp(i)
             write(idet,fmt) -c_disp(k,i),(det_disp(n,k,i),n=1,nmo_disp)
          enddo
          close(idet)
       enddo
    endif
       
    return
    
  end subroutine rephase
    
!######################################################################

  subroutine trans_refpsi

    use constants
    use iomod
    use bdglobal
    
    implicit none

    integer                        :: i,j,n
    real(dp), dimension(nsta,nsta) :: tmparr

!-----------------------------------------------------------------------
! If needed, read the transformation matrix from a previous log file
!-----------------------------------------------------------------------
    if (lrdreftrans) call rdreftrans
    
!-----------------------------------------------------------------------
! MGS orthonormalisation of the transformation to get rid of any
! issuses arising from the finite precision input of the transformation
! matrix
!-----------------------------------------------------------------------
    do i=1,nsta
       do j=1,i-1
          reftrans(:,i)=reftrans(:,i)&
               -dot_product(reftrans(:,i),reftrans(:,j))*reftrans(:,j)
       enddo
       reftrans(:,i)=reftrans(:,i)&
            /sqrt(dot_product(reftrans(:,i),reftrans(:,i)))
    enddo
    
!-----------------------------------------------------------------------
! Transformation of the wavefunction overlap matrix
!-----------------------------------------------------------------------
    tmparr=matmul(spsi,transpose(reftrans))
    spsi=tmparr
    
    return

  end subroutine trans_refpsi

!######################################################################

  subroutine rdreftrans

    use constants
    use channels
    use iomod
    use bdglobal
    
    implicit none

    integer            :: unit,i,j,itmp,jtmp
    character(len=120) :: string
    logical            :: found
    
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
! Read the transformation matrix from the old log file
!-----------------------------------------------------------------------
5   read(unit,'(a)',end=999) string
    if (index(string,'ADT Matrix').eq.0) goto 5

    read(unit,*)

    do i=1,nsta
       do j=1,nsta
          read(unit,*) itmp,jtmp,reftrans(i,j)
       enddo
    enddo
    
!-----------------------------------------------------------------------
! Close the previous log file
!-----------------------------------------------------------------------
    close(unit)
    
    return

999 continue
    errmsg='The ADT matrix section could not be found in: '&
         //trim(areftrans)
    call error_control
    
  end subroutine rdreftrans

!######################################################################

  subroutine switch_diabats

    use constants
    use channels
    use bdglobal
    
    implicit none

    integer                        :: i,j,ilbl
    real(dp), dimension(nsta,nsta) :: tau,tmpmat
    real(dp), dimension(nsta)      :: tmpvec
    real(dp)                       :: mxv
    
!----------------------------------------------------------------------
! Overlaps between the quasi-diabatic states at ref. and disp.
! geometries: tau_ij = < i_disp | j_ref >
!
! Note that if we are here, then the ref. wavefunctions in spsi
! have already been transformed using the reftrans transformation.
! i.e., we only need to transform the disp. wavefunctions
!----------------------------------------------------------------------
    tau=matmul(adt,spsi)
    
!----------------------------------------------------------------------
! Analysis of the ref. - disp. quasi-diabatic wavefunction overlap
! matrix
!----------------------------------------------------------------------
    write(ilog,'(/,82a)') ('+',i=1,82)
    write(ilog,'(2x,a)') 'Quasi-Diabatic Wavefunction Overlaps'
    write(ilog,'(82a)') ('+',i=1,82)

    write(ilog,'(/,47a)') ('-',i=1,47)
    write(ilog,'(a)') '  Disp. State  |  Ref. State  |  Overlap'
    write(ilog,'(a)') '     |I>       |    |J>       |   <I|J>'
    write(ilog,'(47a)') ('-',i=1,47)
    
    do i=1,nsta
       do j=1,nsta
          if (i.ne.j.and.abs(tau(i,j)).gt.0.8d0) then
             write(ilog,'(5x,i2,13x,i2,11x,F13.10,2x,a)') i,j,&
                  tau(i,j),'*'
          else
             write(ilog,'(5x,i2,13x,i2,11x,F13.10)') i,j,tau(i,j)
          endif
       enddo
    enddo

!----------------------------------------------------------------------
! Check: can we fix things by transposing states?
!----------------------------------------------------------------------
    do i=1,nsta
       mxv=0.0d0
       do j=1,nsta
          if (abs(tau(i,j)).gt.mxv) then
             mxv=abs(tau(i,j))
             ilbl=j
          endif
       enddo
       tmpmat(i,:)=adt(ilbl,:)
    enddo
    adt=tmpmat
    
    return
    
  end subroutine switch_diabats
    
!######################################################################
  
  subroutine write_adt

    use constants
    use channels
    use iomod
    use bdglobal
    
    implicit none

    integer :: i,j
    
!----------------------------------------------------------------------
! Write the ADT matrix to the log file
!----------------------------------------------------------------------
    write(ilog,'(/,82a)') ('+',i=1,82)
    write(ilog,'(2x,a)') 'ADT Matrix'
    write(ilog,'(82a)') ('+',i=1,82)

    do i=1,nsta
       do j=1,nsta
          write(ilog,'(2(2x,i2),2x,F15.10)') i,j,adt(i,j)
       enddo
    enddo
    
    return
    
  end subroutine write_adt
  
!######################################################################

  subroutine diabpotmat

    use constants
    use channels
    use iomod
    use bdglobal
    
    implicit none

    integer :: i,j
    
!----------------------------------------------------------------------
! Calculate the quasi-diabatic potential matrix
!----------------------------------------------------------------------
    Wmat=matmul(adt,matmul(Vmat,transpose(adt)))
    
!----------------------------------------------------------------------
! Write the quasi-diabatic potential matrix to the log file
!----------------------------------------------------------------------
    write(ilog,'(/,82a)') ('+',i=1,82)
    write(ilog,'(2x,a)') 'Quasi-Diabatic Potential Matrix'
    write(ilog,'(82a)') ('+',i=1,82)

    do i=1,nsta
       do j=i,nsta
          write(ilog,'(2(2x,i2),2x,F15.10)') i,j,Wmat(i,j)
       enddo
    enddo

    return
    
  end subroutine diabpotmat

!######################################################################

  subroutine diabdipmat

    use constants
    use channels
    use iomod
    use bdglobal
    
    implicit none

    integer :: i,j

!----------------------------------------------------------------------
! Calculate the quasi-diabatic dipole matrix
!----------------------------------------------------------------------
    do i=1,3
       ddip(:,:,i)=matmul(adt,matmul(adip(:,:,i),transpose(adt)))
    enddo

!----------------------------------------------------------------------
! Write the quasi-diabatic potential matrix to the log file
!----------------------------------------------------------------------
    write(ilog,'(/,82a)') ('+',i=1,82)
    write(ilog,'(2x,a)') 'Quasi-Diabatic Dipole Matrix'
    write(ilog,'(82a)') ('+',i=1,82)

    write(ilog,'(2x,3(16x,a1))') 'x','y','z'

    do i=1,nsta
       do j=i,nsta
          write(ilog,'(2(2x,i2),3(2x,F15.10))') i,j,ddip(i,j,:)
       enddo
    enddo
        
    return
    
  end subroutine diabdipmat
    
!######################################################################
    
end program blockdiag


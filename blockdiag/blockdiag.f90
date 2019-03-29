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
  use overlaps
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
! Adjust the phases of the wavefunctions to try and achieve
! consistency across geometries
!----------------------------------------------------------------------
  call rephase
  
!----------------------------------------------------------------------
! Calculate the overlaps between the MOs at the reference and
! displaced geometries
!----------------------------------------------------------------------
  call mo_overlaps

!----------------------------------------------------------------------
! Calculate the overlaps between the electronic states at the
! referenceand displaced geometries
!----------------------------------------------------------------------
  call psi_overlaps

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
! Determine the number of states from the input file
!----------------------------------------------------------------------
    rewind(iin)

    ! Determine the no. states from the $dets_ref section
    nsta=0
5   call rdinp(iin)
    if (keyword(1).eq.'$dets_ref') then
10     call rdinp(iin)
       if (keyword(1).ne.'$end') then
          nsta=nsta+1
          goto 10
       endif
    endif
    if (.not.lend) goto 5

    ! Quit if the $dets_ref section is missing
    if (nsta.eq.0) then
       errmsg='The $dets_ref section could not be found in '//trim(ain)
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

    integer :: i,k,l
    
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
    ioverlap=1
    
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
             read(keyword(1),*) Vmat(k,k)
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
    do i=1,nsta
       if (adetref(i).eq.'') then
          write(errmsg,'(a,x,i2)') 'The reference geometry &
               determinant file is missing for state',i
          call error_control
       endif
    enddo

    ! Displaced geometry determinant file
    do i=1,nsta
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
    allocate(adetref(nsta))
    allocate(adetdisp(nsta))
    adetref=''
    adetdisp=''
    
    ! Number of determinants for each state at the reference geometry
    allocate(ndet_ref(nsta))
    ndet_ref=0

    ! Number of determinants for each state at the displaced geometry
    allocate(ndet_disp(nsta))
    ndet_disp=0

    ! Reference geometry transformation matrix
    allocate(reftrans(nsta,nsta))
    reftrans=0.0d0
    
    ! Norms of the wavefunctions at the reference geometry
    allocate(norm_ref(nsta))
    norm_ref=0.0d0

    ! Norms of the wavefunctions at the displaced geometry
    allocate(norm_disp(nsta))
    norm_disp=0.0d0
    
    ! Overlaps between the electronic states at the two geometries
    allocate(spsi(nsta,nsta))
    spsi=0.0d0

    ! Adiabatic potential matrix
    allocate(Vmat(nsta,nsta))
    Vmat=0.0d0

    ! Diabatic potential matrix
    allocate(Wmat(nsta,nsta))
    Wmat=0.0d0

    ! ADT matrix
    allocate(adt(nsta,nsta))
    adt=0.0d0
    
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
    deallocate(Wmat)
    deallocate(adt)
    deallocate(c_ref)
    deallocate(c_disp)
    deallocate(det_ref)
    deallocate(det_disp)
    
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
    use parsemod
    use bdglobal
    use utils
    use timingmod
    
    implicit none

    integer  :: i,k,n,idet
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
! First pass: determine the no. determinants for each file
!-----------------------------------------------------------------------
    do i=1,nsta
       ndet_ref(i)=nlines(adetref(i))
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
    allocate(c_ref(maxdet,nsta))
    allocate(c_disp(maxdet,nsta))
    c_ref=0.0d0
    c_disp=0.0d0

    ! Determinant vectors
    allocate(det_ref(nmo_ref,maxdet,nsta))
    allocate(det_disp(nmo_disp,maxdet,nsta))
    det_ref=0
    det_disp=0

!-----------------------------------------------------------------------
! Read in the determinants and coefficients
!-----------------------------------------------------------------------
    call freeunit(idet)

    ! Reference geometry
    do i=1,nsta
       open(idet,file=adetref(i),form='formatted',status='old')
       do k=1,ndet_ref(i)
          call rdinp(idet)
          read(keyword(1),*) c_ref(k,i)
          do n=2,inkw
             read(keyword(n),*) det_ref(n-1,k,i)
          enddo          
       enddo
    enddo
    
    ! Displaced geometry
    do i=1,nsta
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

!-----------------------------------------------------------------------
! Norms
!-----------------------------------------------------------------------
    ! Reference geometry
    do i=1,nsta
       norm_ref(i)=sqrt(sum(c_ref(:,i)**2))
    enddo

    ! Displaced geometry
    do i=1,nsta
       norm_disp(i)=sqrt(sum(c_disp(:,i)**2))
    enddo
    
!-----------------------------------------------------------------------
! Normalisation of the wavefunctions
!-----------------------------------------------------------------------
    do i=1,nsta
       c_ref(:,i)=c_ref(:,i)/norm_ref(i)
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
    do i=1,nsta
       write(ilog,'(5x,i2,11x,F13.10)') i,norm_ref(i)
    enddo
    
!----------------------------------------------------------------------
! Displaced geometry wavefunctions
!----------------------------------------------------------------------
    write(ilog,'(/,47a)') ('-',i=1,47)
    write(ilog,'(a)') '  Disp. State  |   Norm'
    write(ilog,'(a)') '     |I>       | || |I> ||'
    write(ilog,'(47a)') ('-',i=1,47)
    do i=1,nsta
       write(ilog,'(5x,i2,11x,F13.10)') i,norm_disp(i)
    enddo

    return
    
  end subroutine wrnorms
    
!######################################################################

  subroutine rephase

    use constants
    use iomod
    use bdglobal
    use utils

    implicit none

    integer               :: i,n
    integer, allocatable  :: indx(:)
    real(dp), parameter   :: thrsh=1e-8_dp
    real(dp)              :: ftmp1,ftmp2
    real(dp), allocatable :: absval(:)

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
    allocate(indx(maxdet))
    indx=0

    allocate(absval(maxdet))
    absval=0.0d0
    
!-----------------------------------------------------------------------
! Adjustment of the phases of the wavefunctions at both the reference
! and displaced geometries
!-----------------------------------------------------------------------
    ! Reference geometry wavefunctions
    do i=1,nsta

       ! Number of determinants for the current wavefunction
       n=ndet_ref(i)

       ! Sort the absolute values of the coefficients for the current
       ! wavefunction in ascending order
       absval=abs(c_ref(:,i))
       call dsortindxa1('D',n,absval(1:n),indx(1:n))
       
       ! Adjust the phase of the wavefunction if the largest
       ! coefficient is not positive
       if (c_ref(indx(1),i).lt.0.0d0) c_ref(:,i)=-c_ref(:,i)
       
       ! Sanity check: exit here if the two leading coefficients are
       ! equal in magnitude but oposite in sign
       ftmp1=sign(1.0d0,c_ref(indx(1),i))*sign(1.0d0,c_ref(indx(2),i))
       ftmp2=abs(c_ref(indx(1),i))-abs(c_ref(indx(2),i))
       if (ftmp1.lt.0.0d0.and.ftmp2.le.thrsh) then
          errmsg='Something terrible has happened in subroutine &
               rephase...'
          call error_control
       endif
       
    enddo
    
    ! Displaced geometry wavefunctions
    do i=1,nsta

       ! Number of determinants for the current wavefunction
       n=ndet_disp(i)

       ! Sort the absolute values of the coefficients for the current
       ! wavefunction in ascending order
       absval=abs(c_disp(:,i))
       call dsortindxa1('D',n,absval(1:n),indx(1:n))

       ! Adjust the phase of the wavefunction if the largest
       ! coefficient is not positive
       if (c_disp(indx(1),i).lt.0.0d0) c_disp(:,i)=-c_disp(:,i)
       
       ! Sanity check: exit here if the two leading coefficients are
       ! equal in magnitude but oposite in sign
       ftmp1=sign(1.0d0,c_disp(indx(1),i))*sign(1.0d0,c_disp(indx(2),i))
       ftmp2=abs(c_disp(indx(1),i))-abs(c_disp(indx(2),i))
       if (ftmp1.lt.0.0d0.and.ftmp2.le.thrsh) then
          errmsg='Something terrible has happened in subroutine &
               rephase...'
          call error_control
       endif
       
    enddo

!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
    deallocate(indx)
    deallocate(absval)
    
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
    tmparr=matmul(spsi,reftrans)
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
    tau=matmul(transpose(adt),spsi)
    
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
    Wmat=matmul(transpose(adt),matmul(Vmat,adt))
    
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
    
end program blockdiag


!######################################################################
! overlaps: routines used in the calculation of the overlaps of
!           electronic states
!######################################################################

module overlaps

  implicit none
  
contains

!######################################################################

  subroutine mo_overlaps

    use constants
    use channels
    use bdglobal
    use import_gamess
    use timingmod
    
    implicit none

    integer               :: i,j,nao_ref,nao_disp
    real(dp), allocatable :: sao(:,:),ao2mo_ref(:,:),ao2mo_disp(:,:)
    real(dp)              :: tw1,tw2,tc1,tc2

!----------------------------------------------------------------------
! Output what we are doing
!----------------------------------------------------------------------
    write(ilog,'(/,82a)') ('+',i=1,82)
    write(ilog,'(2x,a)') 'Calculating MO Overlaps'
    write(ilog,'(82a)') ('+',i=1,82)
    
!-----------------------------------------------------------------------
! Start timing
!-----------------------------------------------------------------------
    call times(tw1,tc1)
    
!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
    ! Numbers of AOs
    nao_ref=gam_ref%nbasis
    nao_disp=gam_disp%nbasis

    ! AO overlap matrix
    allocate(sao(nao_ref,nao_disp))
    sao=0.0d0
    
    ! MO overlap matrix
    allocate(smo(nmo_ref,nmo_disp))
    smo=0.0d0
    
    ! AO-to-MO transformation matrices
    allocate(ao2mo_ref(nao_ref,nmo_ref))
    allocate(ao2mo_disp(nao_disp,nmo_disp))
    ao2mo_ref=0.0d0
    ao2mo_disp=0.0d0
    
!-----------------------------------------------------------------------
! AO-to-MO transformation matrices
!-----------------------------------------------------------------------
    ao2mo_ref=gam_ref%vectors(1:nao_ref,1:nmo_ref)
    ao2mo_disp=gam_disp%vectors(1:nao_disp,1:nmo_disp)

!-----------------------------------------------------------------------
! Calculation of the displaced (bra) - reference (ket) AO overlap matrix
!-----------------------------------------------------------------------
    call gamess_1e_integrals('AO OVERLAP',sao,gam_disp,gam_ref)

!-----------------------------------------------------------------------
! Calculation of the displaced (bra) - reference (ket) MO overlap matrix
!-----------------------------------------------------------------------
    smo=matmul(transpose(ao2mo_disp),matmul(sao,ao2mo_ref))

!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
    deallocate(sao)
    deallocate(ao2mo_ref)
    deallocate(ao2mo_disp)

!-----------------------------------------------------------------------    
! Output timings
!-----------------------------------------------------------------------    
    call times(tw2,tc2)
    write(ilog,'(/,a,1x,F9.2,1x,a)') 'Total Wall Time For MO Overlaps:'&
         ,tw2-tw1," s"
    write(ilog,'(a,2x,F9.2,1x,a)') 'Total CPU Time For MO Overlaps:',&
         tc2-tc1," s"
    
    return
    
  end subroutine mo_overlaps

!######################################################################

  subroutine psi_overlaps

    use constants
    use channels
    use bdglobal
    use timingmod
    
    implicit none

    integer             :: i,j
    real(dp), parameter :: ovrthrsh=0.5d0
    real(dp)            :: tw1,tw2,tc1,tc2
    logical             :: lovrlp

!----------------------------------------------------------------------
! Output what we are doing
!----------------------------------------------------------------------
    write(ilog,'(/,82a)') ('+',i=1,82)
    write(ilog,'(2x,a)') 'Calculating Adiabatic Wavefunction Overlaps'
    write(ilog,'(82a)') ('+',i=1,82)

!-----------------------------------------------------------------------
! Start timing
!-----------------------------------------------------------------------
    call times(tw1,tc1)
    
!----------------------------------------------------------------------
! Call to the requested wavefunction overlap code:
! 1 <-> fast, but memory intensive
! 2 <-> slow, but memory efficient
!----------------------------------------------------------------------
    if (ioverlap.eq.1) then
       call psi_overlaps_fast
    else if (ioverlap.eq.2) then
       call psi_overlaps_slow
    endif

!----------------------------------------------------------------------
! Print out a warning if it looks like a state from outside the group
! of interest has crossed in
!----------------------------------------------------------------------
    ! Loop over disp. states
    do i=1,nsta

       lovrlp=.false.
       
       ! Loop over ref. states
       do j=1,nsta
          if (abs(spsi(i,j)).ge.ovrthrsh) lovrlp=.true.
       enddo
       
       if (.not.lovrlp) then
          write(ilog,'(/,2x,a)') 'WARNING: state crossing detected!'
          !spsi(i,:)=0.0d0
          !spsi(:,i)=0.0d0
          !spsi(i,i)=1.0d0
       endif
       
    enddo

!-----------------------------------------------------------------------
! Lowdin orthogonalisation of the overlap matrix
!-----------------------------------------------------------------------
!    call lowdin_ortho
    
!-----------------------------------------------------------------------
! Output timings
!-----------------------------------------------------------------------
    call times(tw2,tc2)
    write(ilog,'(/,a,1x,F9.2,1x,a)') &
         'Total Wall Time For Wavefunction Overlaps:',tw2-tw1," s"
    write(ilog,'(a,2x,F9.2,1x,a)') &
         'Total CPU Time For Wavefunction Overlaps:',tc2-tc1," s"

    return
    
  end subroutine psi_overlaps

!#####################################################################

  subroutine psi_overlaps_fast

    use constants
    use channels
    use bdglobal
    use timingmod
    use omp_lib
    
    implicit none
    
    integer               :: i,j,m,k,nthreads,tid,iad,ibd,iar,ibr
    real(dp), allocatable :: spsi_1thread(:,:,:)
    real(dp)              :: tw1,tw2,tc1,tc2

!-----------------------------------------------------------------------
! Number of threads
!-----------------------------------------------------------------------  
    !$omp parallel
    nthreads=omp_get_num_threads()
    !$omp end parallel

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(spsi_1thread(nsta,nsta,nthreads))
    spsi_1thread=0.0d0
    
!----------------------------------------------------------------------
! Get the alpha and beta spinorbital indices for every determinant
!----------------------------------------------------------------------
    if (.not.allocated(iocca_ref)) call alpha_beta_indices

!----------------------------------------------------------------------
! Generate an integer label for every unique alpha and beta string
!----------------------------------------------------------------------
    call alpha_beta_labels

!----------------------------------------------------------------------
! Sort the alpha and beta strings
!----------------------------------------------------------------------
    call alpha_beta_sort

!----------------------------------------------------------------------
! Calculate the unique alpha and beta factors
!----------------------------------------------------------------------
    call get_unique_factors

!-----------------------------------------------------------------------
! Start timing for the contraction step
!-----------------------------------------------------------------------
    call times(tw1,tc1)
    
!----------------------------------------------------------------------
! Table header
!----------------------------------------------------------------------
    write(ilog,'(/,47a)') ('-',i=1,47)
    write(ilog,'(a)') '  Disp. State  |  Ref. State  |  Overlap'
    write(ilog,'(a)') '     |I>       |    |J>       |   <I|J>'
    write(ilog,'(47a)') ('-',i=1,47)
    
!----------------------------------------------------------------------
! Calculate overlaps
!----------------------------------------------------------------------
    spsi=0.0d0

    ! Loop disp. states
    do i=1,nsta
       ! Loop over ref. states
       do j=1,nsta

          !$omp parallel do &
          !$omp& firstprivate(i,j) &
          !$omp& private(m,k,iad,ibd,iar,ibr,tid) &
          !$omp& shared(ndet_ref,ndet_disp,afac,bfac,ia_disp,ib_disp,&
          !$omp&        ia_ref,ib_ref,c_ref,c_disp,spsi_1thread)
          
          ! Loop over determinants for the displaced geometry
          do m=1,ndet_disp(i)

             tid=1+omp_get_thread_num()
             
             ! Indices of the unique alpha and beta strings
             ! for the current disp. state/determinat pair
             iad=ia_disp(m,i)
             ibd=ib_disp(m,i)
             
             ! Loop over determinants for the reference geometry
             do k=1,ndet_ref(j)
                
                ! Indices of the unique alpha and beta strings
                ! for the current ref. state/determinat pair
                iar=ia_ref(k,j)
                ibr=ib_ref(k,j)
                
                ! Contibution to < i | j > from the current determinant
                ! pair
!                spsi(i,j)=spsi(i,j)+c_disp(m,i)*c_ref(k,j)&
!                     *afac(iar,iad)*bfac(ibr,ibd)
                spsi_1thread(i,j,tid)=spsi_1thread(i,j,tid)&
                     +c_disp(m,i)*c_ref(k,j)*afac(iar,iad)*bfac(ibr,ibd)
                     
             enddo
             
          enddo
          !$omp end parallel do

          ! Accumulate the contributions from each thread
          spsi(i,j)=sum(spsi_1thread(i,j,:))
          
          ! Table entry
          write(ilog,'(5x,i2,13x,i2,11x,F13.10)') i,j,spsi(i,j)
          
       enddo
    enddo
    
    ! End of the table
    write(ilog,'(47a)') ('-',i=1,47)

!-----------------------------------------------------------------------    
! Output timings for the contraction step
!-----------------------------------------------------------------------    
    call times(tw2,tc2)
    write(ilog,'(/,2x,a,1x,F9.2,1x,a)') 'Wall Time For Contraction:'&
         ,tw2-tw1," s"
    write(ilog,'(2x,a,2x,F9.2,1x,a)') 'CPU Time For Contraction:',&
         tc2-tc1," s"
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(spsi_1thread)

    deallocate(afac)
    deallocate(bfac)
    deallocate(stringa_ref)
    deallocate(stringb_ref)
    deallocate(stringa_disp)
    deallocate(stringb_disp)
    deallocate(ia_ref)
    deallocate(ib_ref)
    deallocate(ia_disp)
    deallocate(ib_disp)
        
    return
  
  end subroutine psi_overlaps_fast

!#####################################################################

  subroutine psi_overlaps_slow

    use constants
    use channels
    use bdglobal
    use omp_lib
    
    implicit none

    integer               :: i,j,m,k,nthreads,tid
    real(dp), allocatable :: smk(:,:)
    real(dp), allocatable :: spsi_1thread(:,:,:)
    real(dp)              :: detsmk

!-----------------------------------------------------------------------
! Number of threads
!-----------------------------------------------------------------------  
    !$omp parallel
    nthreads=omp_get_num_threads()
    !$omp end parallel

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    nel=sum(abs(det_ref(:,5,1)))
    allocate(smk(nel,nel))
    smk=0.0d0

    allocate(spsi_1thread(nsta,nsta,nthreads))
    spsi_1thread=0.0d0

!----------------------------------------------------------------------
! Table header
!----------------------------------------------------------------------
    write(ilog,'(/,47a)') ('-',i=1,47)
    write(ilog,'(a)') '  Disp. State  |  Ref. State  |  Overlap'
    write(ilog,'(a)') '     |I>       |    |J>       |   <I|J>'
    write(ilog,'(47a)') ('-',i=1,47)
    
!----------------------------------------------------------------------
! Calculate overlaps
!----------------------------------------------------------------------
    spsi=0.0d0
    
    ! Loop disp. states
    do i=1,nsta
       ! Loop over ref. states
       do j=1,nsta
          
          !$omp parallel do &
          !$omp& private(m,k,tid,smk,detsmk) &
          !$omp& shared(ndet_ref,ndet_disp,det_ref,det_disp,&
          !$omp&        smo,c_ref,c_disp,spsi_1thread)

          ! Loop over determinants for the displaced geometry
          do m=1,ndet_disp(i)
             
             ! Loop over determinants for the reference geometry
             do k=1,ndet_ref(j)

                tid=1+omp_get_thread_num()
                
                ! Construct the matrix S^mk of overlaps between
                ! the MOs occupied in the displaced bra <m| and the
                ! reference ket |k>
                call fill_spinorbital_integrals(det_disp(:,m,i),&
                    det_ref(:,k,j),smk,smo)

                ! Calculate det S^mk
                detsmk=determinant_overlap(smk)

                ! Add the contribution to < i | j >
                spsi_1thread(i,j,tid)=spsi_1thread(i,j,tid)&
                     +detsmk*c_disp(m,i)*c_ref(k,j)

!                spsi(i,j)=spsi(i,j)+detsmk*c_disp(m,i)*c_ref(k,j)
                
             enddo
                
          enddo
          !$omp end parallel do

          ! Accumulate the contributions from each thread
          spsi(i,j)=sum(spsi_1thread(i,j,:))

          ! Table entry
          write(ilog,'(5x,i2,13x,i2,11x,F13.10)') i,j,spsi(i,j)
          
       enddo
    enddo

    ! End of the table
    write(ilog,'(47a)') ('-',i=1,47)

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(smk)
    deallocate(spsi_1thread)
    
    return
    
  end subroutine psi_overlaps_slow
    
!#####################################################################

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
    allocate(iocca_ref(nalpha,maxdet,nsta))
    allocate(ioccb_ref(nbeta,maxdet,nsta))
    allocate(iocca_disp(nalpha,maxdet,nsta))
    allocate(ioccb_disp(nbeta,maxdet,nsta))
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
    do i=1,nsta
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
    do i=1,nsta
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

!#####################################################################

  subroutine alpha_beta_labels

    use constants
    use channels
    use bdglobal
    
    implicit none

    integer                 :: i,k
    character(len=nalpha*3) :: string_alpha    
    character(len=nbeta*3)  :: string_beta
    character(len=10)       :: fmat_alpha,fmat_beta

!----------------------------------------------------------------------
! Allocate and initialise the label arrays
!----------------------------------------------------------------------
    allocate(ilbla_ref(maxdet,nsta))
    allocate(ilblb_ref(maxdet,nsta))
    allocate(ilbla_disp(maxdet,nsta))
    allocate(ilblb_disp(maxdet,nsta))
    ilbla_ref=0
    ilblb_ref=0
    ilbla_disp=0
    ilblb_disp=0
    
!----------------------------------------------------------------------
! Generate a hash of every alpha and beta spinorbital string
!----------------------------------------------------------------------
    ! Format statements
    fmat_alpha=''
    fmat_beta=''
    write(fmat_alpha,'(a,i0,a)') '(',nalpha,'i0)'
    write(fmat_beta,'(a,i0,a)') '(',nbeta,'i0)'
    
    ! Ref. states
    !
    do i=1,nsta
       do k=1,ndet_ref(i)

          ! Alpha spinorbital character string
          string_alpha=''
          write(string_alpha,fmat_alpha) iocca_ref(:,k,i)

          ! Beta spinorbital character string
          string_beta=''
          write(string_beta,fmat_beta) ioccb_ref(:,k,i)

          ! Hashes of the alpha and beta spinorbital character string
          ilbla_ref(k,i)=djb_hash(trim(string_alpha))
          ilblb_ref(k,i)=djb_hash(trim(string_beta))
          
       enddo
    enddo

    ! Disp. states
    !
    do i=1,nsta
       do k=1,ndet_disp(i)

          ! Alpha spinorbital character string
          string_alpha=''
          write(string_alpha,fmat_alpha) iocca_disp(:,k,i)

          ! Beta spinorbital character string
          string_beta=''
          write(string_beta,fmat_beta) ioccb_disp(:,k,i)
          
          ! Hashes of the alpha and beta spinorbital character string
          ilbla_disp(k,i)=djb_hash(trim(string_alpha))
          ilblb_disp(k,i)=djb_hash(trim(string_beta))

       enddo
    enddo
    
    return
    
  end subroutine alpha_beta_labels

!#####################################################################
  
  function djb_hash(str) result(hash)

    implicit none
    
    character(len=*),intent(in) :: str
    integer*8                   :: hash
    integer                     :: i
    
    hash = 5381

    do i=1,len(str)
        hash=(ishft(hash,5)+hash)+ichar(str(i:i))
    enddo

    return
    
  end function djb_hash

!#####################################################################

  subroutine alpha_beta_sort

    use constants
    use channels
    use utils
    use bdglobal

    implicit none

    real(dp) :: mem
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    ! Indices of the unique alpha and beta strings for every
    ! determinant in the ref. and disp. states
    allocate(ia_ref(maxdet,nsta))
    allocate(ib_ref(maxdet,nsta))
    allocate(ia_disp(maxdet,nsta))
    allocate(ib_disp(maxdet,nsta))
    ia_ref=0
    ib_ref=0
    ia_disp=0
    ib_disp=0
    
!----------------------------------------------------------------------
! Determine the unique alpha and beta strings and the
! determinant/state pairs that they correspond to for the ref. states
! and the disp. states
!----------------------------------------------------------------------
    ! alpha, ref.
    call get_unique_strings(ndet_ref,ilbla_ref,na_ref,stringa_ref,&
         iocca_ref,ia_ref,nalpha)

    ! beta, ref.
    call get_unique_strings(ndet_ref,ilblb_ref,nb_ref,stringb_ref,&
         ioccb_ref,ib_ref,nbeta)
    
    ! alpha, disp.
    call get_unique_strings(ndet_disp,ilbla_disp,na_disp,stringa_disp,&
         iocca_disp,ia_disp,nalpha)
    
    ! beta, disp.
    call get_unique_strings(ndet_disp,ilblb_disp,nb_disp,stringb_disp,&
         ioccb_disp,ib_disp,nbeta)

!----------------------------------------------------------------------
! Output the no. unique alpha and beta strings
!----------------------------------------------------------------------
    write(ilog,'(/,2x,a,2x,i0)') &
         'Unique disp. state alpha strings:',na_disp
    write(ilog,'(2x,a,3x,i0)') &
         'Unique disp. state beta strings:',nb_disp
    write(ilog,'(2x,a,3x,i0)') &
         'Unique ref. state alpha strings:',na_ref
    write(ilog,'(2x,a,4x,i0)') &
         'Unique ref. state beta strings:',na_ref

!----------------------------------------------------------------------
! Output the ammount of memory required to store the unique alpha
! and beta factors
!----------------------------------------------------------------------
    mem=8.0d0*(na_ref*na_disp+nb_ref*nb_disp)/1024.0d0**2
    write(ilog,'(/,2x,a,2x,F8.1,x,a)') 'Memory required to store the &
         unique factors:',mem,'MB'

!----------------------------------------------------------------------
! Deallocate arrays that we no longer need
!----------------------------------------------------------------------
    deallocate(ilbla_ref)
    deallocate(ilblb_ref)
    deallocate(ilbla_disp)
    deallocate(ilblb_disp)

    return
    
  end subroutine alpha_beta_sort

!#####################################################################

  subroutine get_unique_strings(ndet,ilbl,nunique,string,iocc,&
       iunique,nspin)

    use constants
    use channels
    use utils
    use bdglobal
    
    implicit none

    integer, dimension(nsta)              :: ndet
    integer*8, dimension(maxdet,nsta)     :: ilbl
    integer                               :: nunique
    integer, allocatable                  :: string(:,:)
    integer, dimension(nspin,maxdet,nsta) :: iocc
    integer, dimension(maxdet,nsta)       :: iunique
    integer                               :: nspin

    integer                               :: i,k,k1,i1,i2,n,&
                                             istate,idet
    integer*8                             :: sumdet,last
    integer*8, allocatable                :: indx(:),ilbl_all(:)
    integer, allocatable                  :: info(:,:)

    ! CHECK
    integer, allocatable :: ipos(:)
    integer              :: m,istate1,idet1,istate2,idet2
    ! CHECK
    
!----------------------------------------------------------------------
! Allocate and initialise arrays
!----------------------------------------------------------------------
    sumdet=sum(ndet(:))
    allocate(indx(sumdet))
    indx=0
    allocate(info(sumdet,2))
    info=0
    allocate(ilbl_all(sumdet))
    ilbl_all=0

!----------------------------------------------------------------------
! Put together the combined list of alpha/beta string labels for the
! current set of states
!----------------------------------------------------------------------
    i1=0
    do i=1,nsta
       ! Next lot of alpha/beta strings
       i2=sum(ndet(1:i))
       i1=i2-ndet(i)+1
       ilbl_all(i1:i2)=ilbl(1:ndet(i),i)
       ! State number of this lot of alpha/beta strings
       info(i1:i2,2)=i
       ! Determinant numbers for this lot of alpha/beta strings
       do k=1,ndet(i)
          info(i1-1+k,1)=k
       enddo
    enddo

!----------------------------------------------------------------------
! Sort the combined list of alpha/beta string labels for the current
! set of states
!----------------------------------------------------------------------
    call i8sortindxa1('A',sumdet,ilbl_all,indx)
    
!----------------------------------------------------------------------  
! No. unique alpha/beta strings for the current set of states
!----------------------------------------------------------------------  
    last=0
    nunique=0
    do k=1,sumdet
       k1=indx(k)
       if (ilbl_all(k1).ne.last) then
          last=ilbl_all(k1)
          nunique=nunique+1
       endif
    enddo

!----------------------------------------------------------------------  
! Assign working alpha/beta string indices to all of the determinants
! for all ref. states
!----------------------------------------------------------------------  
    allocate(string(nspin,nunique))

    last=0
    n=0
    do k=1,sumdet

       k1=indx(k)
       
       ! Determinant index for the current state/determinant pair
       idet=info(k1,1)
       
       ! State index for the current state/determinant pair
       istate=info(k1,2)

       ! Are we at the start of a new unique alpha/beta string?
       if (ilbl_all(k1).ne.last) then
          ! If so, increment the unique alpha/beta string counter and
          ! fill in the nth unique alpha/beta string
          last=ilbl_all(k1)
          n=n+1
          ! nth unique alpha/beta string
          string(:,n)=iocc(:,idet,istate)
       endif
       
       ! Index of the unique alpha/beta string for the current
       ! state/determinant pair
       iunique(idet,istate)=n

    enddo

!----------------------------------------------------------------------
! Check for duplicates
!----------------------------------------------------------------------
!    allocate(ipos(nunique))
!
!    ipos=0
!    last=0
!    n=0
!    do k=1,sumdet
!       k1=indx(k)
!       if (ilbl_all(k1).ne.last) then
!          last=ilbl_all(k1)
!          n=n+1
!          ipos(n)=k
!       endif
!    enddo
!
!    do n=1,nunique
!
!       i1=ipos(n)
!       if (n.lt.nunique) then
!          i2=ipos(n+1)-1
!       else
!          i2=sumdet
!       endif
!
!       do i=i1+1,i2
!          idet1=info(indx(i-1),1)
!          istate1=info(indx(i-1),2)
!          idet2=info(indx(i),1)
!          istate2=info(indx(i),2)
!          do m=1,nspin
!             if (iocc(m,idet1,istate1).ne.iocc(m,idet2,istate2)) then
!                print*,'DUPLICATE FOUND!'
!                stop
!             endif
!          enddo
!       enddo
!       
!    enddo
!    
!    deallocate(ipos)
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(indx)
    deallocate(info)
    deallocate(ilbl_all)
    
    return
    
  end subroutine get_unique_strings

!#####################################################################

  subroutine get_unique_factors

    use constants
    use channels
    use utils
    use bdglobal
    use timingmod
    use omp_lib
    
    implicit none

    integer               :: i,j,k,l,mobra,moket,nthreads,tid
    real(dp), allocatable :: Sa(:,:),Sb(:,:)
    real(dp)              :: ubound
    real(dp)              :: tw1,tw2,tc1,tc2
    real(dp), parameter   :: dthrsh=1e-6_dp
    logical               :: lequiv
    
!-----------------------------------------------------------------------
! Start timing
!-----------------------------------------------------------------------
    call times(tw1,tc1)
    
!-----------------------------------------------------------------------
! Number of threads
!-----------------------------------------------------------------------  
    !$omp parallel
    nthreads=omp_get_num_threads()
    !$omp end parallel
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(afac(na_ref,na_disp))
    afac=0.0d0

    allocate(bfac(nb_ref,nb_disp))
    bfac=0.0d0

    allocate(Sa(nalpha,nalpha))
    Sa=0.0d0

    allocate(Sb(nbeta,nbeta))
    Sb=0.0d0
    
!----------------------------------------------------------------------
! Calculate the unique alpha factors
!----------------------------------------------------------------------
    !$omp parallel do &
    !$omp& private(i,j,k,l,moket,mobra,Sa,ubound,tid) &
    !$omp& shared(na_ref,na_disp,nalpha,stringa_ref,stringa_disp,&
    !$omp&        smo,afac)
    
    ! Loop over unique ref. state alpha strings
    do i=1,na_ref
       
       ! Loop over unique disp. state alpha strings
       do j=1,na_disp

          ! Matrix of orbital overlaps
          do k=1,nalpha
             moket=stringa_ref(k,i)
             do l=1,nalpha
                mobra=stringa_disp(l,j)
                Sa(k,l)=smo(mobra,moket)
             enddo
          enddo

          ! Hadamard-inequality screening
          ubound=hadamard_bound(Sa,nalpha)
          if (ubound.gt.dthrsh) then
             ! Determinant of the matrix of orbital overlaps
             afac(i,j)=determinant_new(Sa,nalpha)
          endif

       enddo

    enddo
    
    !$omp end parallel do

!----------------------------------------------------------------------
! Calculate the unique beta factors
!----------------------------------------------------------------------
    ! Check to see if the sets of alpha and beta strings are
    ! equivalent
    lequiv=isequivab()

    if (lequiv) then
       ! The sets of unique alpha and beta strings are equivalent.
       ! We can skip the calculation of the unique beta factors.
       bfac=afac
    else
       ! The sets of unique alpha and beta strings are not equivalent.
       ! We have to explicitly calculate the unique beta factors.
       !
       !$omp parallel do &
       !$omp& private(i,j,k,l,moket,mobra,Sb,tid) &
       !$omp& shared(nb_ref,nb_disp,nbeta,stringb_ref,stringb_disp,&
       !$omp&        smo,bfac)
       !
       ! Loop over unique ref. state beta strings
       do i=1,nb_ref
          
          ! Loop over unique disp. state beta strings
          do j=1,nb_disp
             
             ! Matrix of orbital overlaps
             do k=1,nbeta
                moket=stringb_ref(k,i)
                do l=1,nbeta
                   mobra=stringb_disp(l,j)
                   Sb(k,l)=smo(mobra,moket)
                enddo
             enddo
             
             ! Hadamard-inequality screening
             ubound=hadamard_bound(Sb,nbeta)
             if (ubound.gt.dthrsh) then
                ! Determinant of the matrix of orbital overlaps
                bfac(i,j)=determinant_new(Sb,nbeta)
             endif
             
          enddo
          
       enddo
       !
       !$omp end parallel do
    endif
       
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(Sa)
    deallocate(Sb)

!-----------------------------------------------------------------------    
! Output timings
!-----------------------------------------------------------------------    
    call times(tw2,tc2)
    write(ilog,'(/,2x,a,1x,F9.2,1x,a)') 'Wall Time For Unique Factors:'&
         ,tw2-tw1," s"
    write(ilog,'(2x,a,2x,F9.2,1x,a)') 'CPU Time For Unique Factors:',&
         tc2-tc1," s"
    
    return
    
  end subroutine get_unique_factors

!#####################################################################

  function hadamard_bound(mat,dim) result(func)

    use constants
    
    implicit none

    integer, intent(in)                      :: dim
    integer                                  :: i
    real(dp), dimension(dim,dim), intent(in) :: mat
    real(dp), dimension(dim)                 :: norm
    real(dp)                                 :: func

    ! Norms of the columns of the matrix
    do i=1,dim
       norm(i)=dot_product(mat(:,i),mat(:,i))
    enddo
    norm=sqrt(norm)

    ! Hadamard upper bound for the determinant of the matrix
    func=1.0d0
    do i=1,dim
       func=func*norm(i)
    enddo
    
    return
    
  end function hadamard_bound

!#####################################################################

  function isequivab()

    use constants
    use channels
    use bdglobal
    
    implicit none

    integer :: i,j
    logical :: isequivab

    ! Shortcut 1: The sets of unique alpha and beta strings are not
    ! equivalent if the number of elements in each is not equal
    if (na_disp.ne.nb_disp) then
       isequivab=.false.
       return
    endif

    ! Shortcut2: The sets of unique alpha and beta strings are not
    ! equivalent if the numbers of alpha and beta electronc are
    ! not equal
    if (nalpha.ne.nbeta) then
       isequivab=.false.
       return
    endif

    ! Check on the equivalency of the sets of unique alpha and beta
    ! strings
    isequivab=.true.
    do i=1,na_ref
       do j=1,nalpha
          if (stringa_disp(i,j).ne.stringb_disp(i,j)) then
             isequivab=.false.
             exit
          endif
       enddo
    enddo
    
    return
    
  end function isequivab
  
!#####################################################################
! Calculate 1-particle overlap for all spin-orbitals in a given pair 
! of determinants.
! Taken from the superdyson code (see sd_core.f90).
!######################################################################
  
  subroutine fill_spinorbital_integrals(occbra,occket,adet,amo)

    use constants
    use bdglobal
      
    implicit none

    integer, intent(in)   :: occbra(:) ! Parent determinant
    integer, intent(in)   :: occket(:) ! Ion determinant
    real(dp), intent(out) :: adet(:,:) ! Integral matrix
                                       ! for a given determinant
    real(dp), intent(in)  :: amo (:,:) ! Integral matrix for all MOs
    
    integer :: mobra, moket     ! Spatial MO indices
    integer :: spinbra, spinket ! Spin+space MO indices
    integer :: nalphabra        ! Number of spin-alpha electrons in the bra
    
    integer :: nmoket,nmobra

    nmobra=nmo_disp
    nmoket=nmo_ref
    !
    adet = 0.0d0
    !
    !  Alpha-alpha overlaps - upper left corner of sdet
    !
    spinket = 0
    spinbra = 0
    ket_alpha: do moket=1,nmoket
       if (occket(moket)/=1 .and. occket(moket)/=2) cycle ket_alpha
       !
       spinket = spinket + 1
       spinbra = 0
       bra_alpha: do mobra=1,nmobra
          if (occbra(mobra)/=1 .and. occbra(mobra)/=2) cycle bra_alpha
          !
          spinbra = spinbra + 1
          adet(spinbra,spinket) = amo(mobra,moket)
       end do bra_alpha
    end do ket_alpha
    nalphabra = spinbra
    !
    !  Beta-beta overlaps - bottom right corner of sdet
    !
    ket_beta: do moket=1,nmoket
       if (occket(moket)/=-1 .and. occket(moket)/=2) cycle ket_beta
       !
       spinket = spinket + 1
       spinbra = nalphabra
       bra_beta: do mobra=1,nmobra
          if (occbra(mobra)/=-1 .and. occbra(mobra)/=2) cycle bra_beta
          !
          spinbra = spinbra + 1
          adet(spinbra,spinket) = amo(mobra,moket)
       end do bra_beta
    end do ket_beta
    
    return
      
  end subroutine fill_spinorbital_integrals
  
!######################################################################
! determinant_overlap: Calculates the determinant of the matrix sdet of
!                      overlaps between the MOs in two determinants.
!                      Taken from superdyson.
!######################################################################

  function determinant_overlap(sdet) result(overdet)
      
    use constants
    use bdglobal, only: nel
    
    implicit none
    
    real(dp), intent(in) :: sdet(:,:)
    real(dp)             :: overdet
    real(dp)             :: scr(nel,nel) ! Buffer for the 1-particle 
                                         !  overlap matrix 
    scr     = sdet
    overdet = determinant(scr)
    
    return

  end function determinant_overlap
    
!######################################################################
! determinant: Calculates the determinant of a sqaure matrix using
!              Gaussian elimination.
!######################################################################

  real(dp) function determinant(mat)

    use constants
    
    real(dp), intent(inout) :: mat(:,:) ! Matrix destroyed on exit
    integer                 :: order,info
    integer                 :: ipvt(size(mat,dim=1))
    integer                 :: i
    real(dp)                :: work(size(mat,dim=1))
    real(dp)                :: detx(2)
    
    order = size(mat,dim=1)

    call dgefa(mat,order,order,ipvt,info)    

    call dgedi(mat,order,order,ipvt,detx,work,10)

    determinant = detx(1) * 10.0d0**detx(2)
    
    return

  end function determinant

!######################################################################

  function determinant_new(mat,dim) result(detval)

    use constants

    implicit none

    integer, intent(in)     :: dim
    integer                 :: order,info
    integer                 :: ipvt(dim)
    real(dp)                :: detval
    real(dp), intent(inout) :: mat(dim,dim)
    real(dp)                :: tmp(dim,dim)
    real(dp)                :: work(dim)
    real(dp)                :: detx(2)

!----------------------------------------------------------------------
! The input matrix is overwritten in dgefa/dgedi, so use a copy
!----------------------------------------------------------------------
    tmp=mat

!----------------------------------------------------------------------
! Calculate the determinant of the input matrix using Gaussian
! elimination
!----------------------------------------------------------------------
    call dgefa(tmp,dim,dim,ipvt,info)
    call dgedi(tmp,dim,dim,ipvt,detx,work,10)
    detval=detx(1)*10.0d0**detx(2)
    
    return
    
  end function determinant_new

!######################################################################

  subroutine lowdin_ortho

    use constants
    use channels
    use bdglobal
    use iomod
    
    implicit none

    integer                        :: info,i
    real(dp), dimension(nsta,nsta) :: tmp,umat,vtmat
    real(dp), dimension(nsta)      :: sigma
    real(dp), dimension(5*nsta)    :: work
    
!-----------------------------------------------------------------------
! SVD of the wavefunction overlap matrix
!-----------------------------------------------------------------------
    tmp=spsi
    call dgesvd('A','A',nsta,nsta,tmp,nsta,sigma,umat,nsta,vtmat,nsta,&
         work,5*nsta,info)

    ! Exit if the SVD failed
    if (info.ne.0) then
       errmsg='SVD failed in subroutine lowdin_ortho'
       call error_control
    endif

!-----------------------------------------------------------------------
! Compute the orthogonalised waverfunction overlap matrix
!-----------------------------------------------------------------------
    spsi=matmul(vtmat,umat)
    
    return
    
  end subroutine lowdin_ortho
  
!######################################################################
  
end module overlaps

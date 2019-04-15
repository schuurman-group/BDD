!######################################################################
! wfoverlaps: routines used in the calculation of the overlaps of
!             electronic states
!
!             Note that this module is now completely detached from
!             the global module. This means that we have to pass lots
!             off arrays around, but that we can, e.g., calculate
!             overlaps between different subsets of states.
!######################################################################

module wfoverlaps

  use constants
  
  implicit none

  integer, allocatable   :: ia_bra(:,:),ib_bra(:,:),&
                            ia_ket(:,:),ib_ket(:,:)
  integer, allocatable   :: stringa_bra(:,:),stringb_bra(:,:),&
                            stringa_ket(:,:),stringb_ket(:,:)
  integer*8, allocatable :: ilbla_bra(:,:),ilbla_ket(:,:),&
                            ilblb_bra(:,:),ilblb_ket(:,:)
  integer                :: na_bra,na_ket,nb_bra,nb_ket
  real(dp), allocatable  :: afac(:,:),bfac(:,:)
  
contains

!######################################################################

  subroutine psi_overlaps(spsi,nsta,nalpha,nbeta,ndet_bra,ndet_ket,&
       nmo_bra,nmo_ket,maxdet,c_bra,c_ket,det_bra,det_ket,iocca_bra,&
       iocca_ket,ioccb_bra,ioccb_ket,ioverlap,smo)

    use constants
    use channels
    use iomod
    use timingmod
    
    implicit none

    integer, intent(in)  :: nsta,nalpha,nbeta,nmo_bra,nmo_ket,maxdet,&
                            ioverlap
    integer, intent(in)  :: ndet_bra(nsta),ndet_ket(nsta)
    integer              :: i,j
    real(dp)             :: spsi(nsta,nsta)
    integer, intent(in)  :: det_bra(nmo_bra,maxdet,nsta)
    integer, intent(in)  :: det_ket(nmo_ket,maxdet,nsta)
    integer, intent(in)  :: iocca_bra(nalpha,maxdet,nsta)
    integer, intent(in)  :: iocca_ket(nalpha,maxdet,nsta)
    integer, intent(in)  :: ioccb_bra(nbeta,maxdet,nsta)
    integer, intent(in)  :: ioccb_ket(nbeta,maxdet,nsta)
    real(dp), intent(in) :: c_bra(maxdet,nsta),c_ket(maxdet,nsta)
    real(dp), intent(in) :: smo(nmo_bra,nmo_ket)
    real(dp), parameter  :: ovrthrsh=0.5d0
    real(dp)             :: tw1,tw2,tc1,tc2
    logical              :: lovrlp
    
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
       call psi_overlaps_fast(spsi,nsta,nalpha,nbeta,ndet_bra,&
            ndet_ket,nmo_bra,nmo_ket,maxdet,c_bra,c_ket,iocca_bra,&
            iocca_ket,ioccb_bra,ioccb_ket,smo)
    else if (ioverlap.eq.2) then
       errmsg='The slow overlaps code needs re-writing to detach &
            it from the bdglobal module...'
       call error_control
       !call psi_overlaps_slow
    endif

!----------------------------------------------------------------------
! Print out a warning if it looks like a state from outside the group
! of interest has crossed in
!----------------------------------------------------------------------
    ! Loop over bra states
    do i=1,nsta

       lovrlp=.false.
       
       ! Loop over ket states
       do j=1,nsta
          if (abs(spsi(i,j)).ge.ovrthrsh) lovrlp=.true.
       enddo
       
       if (.not.lovrlp) write(ilog,'(/,2x,a)') &
            'WARNING: state crossing detected!'

    enddo

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

  subroutine psi_overlaps_fast(spsi,nsta,nalpha,nbeta,ndet_bra,&
       ndet_ket,nmo_bra,nmo_ket,maxdet,c_bra,c_ket,iocca_bra,&
       iocca_ket,ioccb_bra,ioccb_ket,smo)

    use constants
    use channels
    use timingmod
    use omp_lib
    
    implicit none

   
    
    integer, intent(in)   :: nsta,nalpha,nbeta,nmo_bra,nmo_ket,maxdet
    integer, intent(in)   :: ndet_bra(nsta),ndet_ket(nsta)
    integer, intent(in)   :: iocca_bra(nalpha,maxdet,nsta)
    integer, intent(in)   :: iocca_ket(nalpha,maxdet,nsta)
    integer, intent(in)   :: ioccb_bra(nbeta,maxdet,nsta)
    integer, intent(in)   :: ioccb_ket(nbeta,maxdet,nsta)
    real(dp)              :: spsi(nsta,nsta)
    real(dp), intent(in)  :: c_bra(maxdet,nsta),c_ket(maxdet,nsta)
    real(dp), intent(in)  :: smo(nmo_bra,nmo_ket)
    
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
! Generate an integer label for every unique alpha and beta string
!----------------------------------------------------------------------
    call alpha_beta_labels(nsta,maxdet,nalpha,nbeta,ndet_bra,ndet_ket,&
         iocca_bra,iocca_ket,ioccb_bra,ioccb_ket)

!----------------------------------------------------------------------
! Sort the alpha and beta strings
!----------------------------------------------------------------------
    call alpha_beta_sort(nsta,maxdet,nalpha,nbeta,ndet_bra,ndet_ket,&
         iocca_bra,iocca_ket,ioccb_bra,ioccb_ket)

!----------------------------------------------------------------------
! Calculate the unique alpha and beta factors
!----------------------------------------------------------------------
    call get_unique_factors(nalpha,nbeta,nmo_bra,nmo_ket,smo)

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

    ! Loop bra states
    do i=1,nsta
       ! Loop over ket states
       do j=1,nsta

          !$omp parallel do &
          !$omp& firstprivate(i,j) &
          !$omp& private(m,k,iad,ibd,iar,ibr,tid) &
          !$omp& shared(ndet_ket,ndet_bra,afac,bfac,ia_bra,ib_bra,&
          !$omp&        ia_ket,ib_ket,c_ket,c_bra,spsi_1thread)
          
          ! Loop over determinants for the bra state
          do m=1,ndet_bra(i)

             tid=1+omp_get_thread_num()
             
             ! Indices of the unique alpha and beta strings
             ! for the current bra state/determinat pair
             iad=ia_bra(m,i)
             ibd=ib_bra(m,i)
             
             ! Loop over determinants for the ket state
             do k=1,ndet_ket(j)
                
                ! Indices of the unique alpha and beta strings
                ! for the current ket state/determinat pair
                iar=ia_ket(k,j)
                ibr=ib_ket(k,j)
                
                ! Contibution to < i | j > from the current determinant
                ! pair
!                spsi(i,j)=spsi(i,j)+c_bra(m,i)*c_ket(k,j)&
!                     *afac(iar,iad)*bfac(ibr,ibd)
                spsi_1thread(i,j,tid)=spsi_1thread(i,j,tid)&
                     +c_bra(m,i)*c_ket(k,j)*afac(iar,iad)*bfac(ibr,ibd)
                     
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
    deallocate(stringa_ket)
    deallocate(stringb_ket)
    deallocate(stringa_bra)
    deallocate(stringb_bra)
    deallocate(ia_ket)
    deallocate(ib_ket)
    deallocate(ia_bra)
    deallocate(ib_bra)
        
    return
  
  end subroutine psi_overlaps_fast

!#####################################################################

!  subroutine psi_overlaps_slow
!
!    use constants
!    use channels
!    use bdglobal
!    use omp_lib
!    
!    implicit none
!
!    integer               :: i,j,m,k,nthreads,tid
!    real(dp), allocatable :: smk(:,:)
!    real(dp), allocatable :: spsi_1thread(:,:,:)
!    real(dp)              :: detsmk
!
!!-----------------------------------------------------------------------
!! Number of threads
!!-----------------------------------------------------------------------  
!    !$omp parallel
!    nthreads=omp_get_num_threads()
!    !$omp end parallel
!
!!----------------------------------------------------------------------
!! Allocate arrays
!!----------------------------------------------------------------------
!    nel=sum(abs(det_ket(:,5,1)))
!    allocate(smk(nel,nel))
!    smk=0.0d0
!
!    allocate(spsi_1thread(nsta,nsta,nthreads))
!    spsi_1thread=0.0d0
!
!!----------------------------------------------------------------------
!! Table header
!!----------------------------------------------------------------------
!    write(ilog,'(/,47a)') ('-',i=1,47)
!    write(ilog,'(a)') '  Disp. State  |  Ref. State  |  Overlap'
!    write(ilog,'(a)') '     |I>       |    |J>       |   <I|J>'
!    write(ilog,'(47a)') ('-',i=1,47)
!    
!!----------------------------------------------------------------------
!! Calculate overlaps
!!----------------------------------------------------------------------
!    spsi=0.0d0
!    
!    ! Loop bra states
!    do i=1,nsta
!       ! Loop over ket states
!       do j=1,nsta
!          
!          !$omp parallel do &
!          !$omp& private(m,k,tid,smk,detsmk) &
!          !$omp& shared(ndet_ket,ndet_bra,det_ket,det_bra,&
!          !$omp&        smo,c_ket,c_bra,spsi_1thread)
!
!          ! Loop over determinants for the bra state
!          do m=1,ndet_bra(i)
!             
!             ! Loop over determinants for the ket state
!             do k=1,ndet_ket(j)
!
!                tid=1+omp_get_thread_num()
!                
!                ! Construct the matrix S^mk of overlaps between
!                ! the MOs occupied in the bra <m| and the ket |k>
!                call fill_spinorbital_integrals(det_bra(:,m,i),&
!                    det_ket(:,k,j),smk,smo)
!
!                ! Calculate det S^mk
!                detsmk=determinant_overlap(smk)
!
!                ! Add the contribution to < i | j >
!                spsi_1thread(i,j,tid)=spsi_1thread(i,j,tid)&
!                     +detsmk*c_bra(m,i)*c_ket(k,j)
!
!!                spsi(i,j)=spsi(i,j)+detsmk*c_bra(m,i)*c_ket(k,j)
!                
!             enddo
!                
!          enddo
!          !$omp end parallel do
!
!          ! Accumulate the contributions from each thread
!          spsi(i,j)=sum(spsi_1thread(i,j,:))
!
!          ! Table entry
!          write(ilog,'(5x,i2,13x,i2,11x,F13.10)') i,j,spsi(i,j)
!          
!       enddo
!    enddo
!
!    ! End of the table
!    write(ilog,'(47a)') ('-',i=1,47)
!
!!----------------------------------------------------------------------
!! Deallocate arrays
!!----------------------------------------------------------------------
!    deallocate(smk)
!    deallocate(spsi_1thread)
!    
!    return
!    
!  end subroutine psi_overlaps_slow
    
!#####################################################################

  subroutine alpha_beta_labels(nsta,maxdet,nalpha,nbeta,ndet_bra,&
       ndet_ket,iocca_bra,iocca_ket,ioccb_bra,ioccb_ket)

    use constants
    use channels
!    use bdglobal
    
    implicit none

    integer, intent(in)     :: nsta,nalpha,nbeta,maxdet
    integer, intent(in)     :: ndet_bra(nsta),ndet_ket(nsta)
    integer                 :: i,k
    integer, intent(in)     :: iocca_bra(nalpha,maxdet,nsta)
    integer, intent(in)     :: iocca_ket(nalpha,maxdet,nsta)
    integer, intent(in)     :: ioccb_bra(nbeta,maxdet,nsta)
    integer, intent(in)     :: ioccb_ket(nbeta,maxdet,nsta)
    character(len=nalpha*3) :: string_alpha    
    character(len=nbeta*3)  :: string_beta
    character(len=10)       :: fmat_alpha,fmat_beta

!----------------------------------------------------------------------
! Allocate and initialise the label arrays
!----------------------------------------------------------------------
    allocate(ilbla_ket(maxdet,nsta))
    allocate(ilblb_ket(maxdet,nsta))
    allocate(ilbla_bra(maxdet,nsta))
    allocate(ilblb_bra(maxdet,nsta))
    ilbla_ket=0
    ilblb_ket=0
    ilbla_bra=0
    ilblb_bra=0
    
!----------------------------------------------------------------------
! Generate a hash of every alpha and beta spinorbital string
!----------------------------------------------------------------------
    ! Format statements
    fmat_alpha=''
    fmat_beta=''
    write(fmat_alpha,'(a,i0,a)') '(',nalpha,'i0)'
    write(fmat_beta,'(a,i0,a)') '(',nbeta,'i0)'
    
    ! Ket states
    !
    do i=1,nsta
       do k=1,ndet_ket(i)

          ! Alpha spinorbital character string
          string_alpha=''
          write(string_alpha,fmat_alpha) iocca_ket(:,k,i)

          ! Beta spinorbital character string
          string_beta=''
          write(string_beta,fmat_beta) ioccb_ket(:,k,i)

          ! Hashes of the alpha and beta spinorbital character string
          ilbla_ket(k,i)=djb_hash(trim(string_alpha))
          ilblb_ket(k,i)=djb_hash(trim(string_beta))
          
       enddo
    enddo

    ! Bra states
    !
    do i=1,nsta
       do k=1,ndet_bra(i)

          ! Alpha spinorbital character string
          string_alpha=''
          write(string_alpha,fmat_alpha) iocca_bra(:,k,i)

          ! Beta spinorbital character string
          string_beta=''
          write(string_beta,fmat_beta) ioccb_bra(:,k,i)
          
          ! Hashes of the alpha and beta spinorbital character string
          ilbla_bra(k,i)=djb_hash(trim(string_alpha))
          ilblb_bra(k,i)=djb_hash(trim(string_beta))

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

  subroutine alpha_beta_sort(nsta,maxdet,nalpha,nbeta,ndet_bra,&
       ndet_ket,iocca_bra,iocca_ket,ioccb_bra,ioccb_ket)

    use constants
    use channels
    use utils
!    use bdglobal

    implicit none
    
    integer, intent(in) :: nsta,nalpha,nbeta,maxdet
    integer, intent(in) :: ndet_bra(nsta),ndet_ket(nsta)
    integer, intent(in) :: iocca_bra(nalpha,maxdet,nsta)
    integer, intent(in) :: iocca_ket(nalpha,maxdet,nsta)
    integer, intent(in) :: ioccb_bra(nbeta,maxdet,nsta)
    integer, intent(in) :: ioccb_ket(nbeta,maxdet,nsta)    
    real(dp)            :: mem
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    ! Indices of the unique alpha and beta strings for every
    ! determinant in the bra and ket states
    allocate(ia_ket(maxdet,nsta))
    allocate(ib_ket(maxdet,nsta))
    allocate(ia_bra(maxdet,nsta))
    allocate(ib_bra(maxdet,nsta))
    ia_ket=0
    ib_ket=0
    ia_bra=0
    ib_bra=0
    
!----------------------------------------------------------------------
! Determine the unique alpha and beta strings and the
! determinant/state pairs that they correspond to for the ket states
! and the bra states
!----------------------------------------------------------------------
    ! alpha, ket
    call get_unique_strings(ndet_ket,ilbla_ket,na_ket,stringa_ket,&
         iocca_ket,ia_ket,nalpha,nsta,maxdet)

    ! beta, ket
    call get_unique_strings(ndet_ket,ilblb_ket,nb_ket,stringb_ket,&
         ioccb_ket,ib_ket,nbeta,nsta,maxdet)
    
    ! alpha, bra
    call get_unique_strings(ndet_bra,ilbla_bra,na_bra,stringa_bra,&
         iocca_bra,ia_bra,nalpha,nsta,maxdet)
    
    ! beta, bra
    call get_unique_strings(ndet_bra,ilblb_bra,nb_bra,stringb_bra,&
         ioccb_bra,ib_bra,nbeta,nsta,maxdet)

!----------------------------------------------------------------------
! Output the no. unique alpha and beta strings
!----------------------------------------------------------------------
    write(ilog,'(/,2x,a,2x,i0)') &
         'Unique bra state alpha strings:',na_bra
    write(ilog,'(2x,a,3x,i0)') &
         'Unique bra state beta strings:',nb_bra
    write(ilog,'(2x,a,3x,i0)') &
         'Unique ket state alpha strings:',na_ket
    write(ilog,'(2x,a,4x,i0)') &
         'Unique ket state beta strings:',na_ket

!----------------------------------------------------------------------
! Output the ammount of memory required to store the unique alpha
! and beta factors
!----------------------------------------------------------------------
    mem=8.0d0*(na_ket*na_bra+nb_ket*nb_bra)/1024.0d0**2
    write(ilog,'(/,2x,a,2x,F8.1,x,a)') 'Memory required to store the &
         unique factors:',mem,'MB'

!----------------------------------------------------------------------
! Deallocate arrays that we no longer need
!----------------------------------------------------------------------
    deallocate(ilbla_ket)
    deallocate(ilblb_ket)
    deallocate(ilbla_bra)
    deallocate(ilblb_bra)

    return
    
  end subroutine alpha_beta_sort

!#####################################################################

  subroutine get_unique_strings(ndet,ilbl,nunique,string,iocc,&
       iunique,nspin,nsta,maxdet)

    use constants
    use channels
    use utils
!    use bdglobal
    
    implicit none

    integer, intent(in)                   :: nsta,maxdet
    integer, intent(in), dimension(nsta)  :: ndet
    integer*8, dimension(maxdet,nsta)     :: ilbl
    integer                               :: nunique
    integer, allocatable                  :: string(:,:)
    integer, intent(in), dimension(nspin,maxdet,nsta) :: iocc
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
! for all ket states
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

  subroutine get_unique_factors(nalpha,nbeta,nmo_bra,nmo_ket,smo)

    use constants
    use channels
    use utils
!    use bdglobal
    use timingmod
    use omp_lib
    
    implicit none

    integer, intent(in)   :: nalpha,nbeta,nmo_bra,nmo_ket
    integer               :: i,j,k,l,mobra,moket,nthreads,tid
    real(dp), intent(in)  :: smo(nmo_bra,nmo_ket)
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
    allocate(afac(na_ket,na_bra))
    afac=0.0d0

    allocate(bfac(nb_ket,nb_bra))
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
    !$omp& shared(na_ket,na_bra,nalpha,stringa_ket,stringa_bra,&
    !$omp&        smo,afac)
    
    ! Loop over unique ket state alpha strings
    do i=1,na_ket
       
       ! Loop over unique bra state alpha strings
       do j=1,na_bra

          ! Matrix of orbital overlaps
          do k=1,nalpha
             moket=stringa_ket(k,i)
             do l=1,nalpha
                mobra=stringa_bra(l,j)
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
    lequiv=isequivab(nalpha,nbeta)

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
       !$omp& shared(nb_ket,nb_bra,nbeta,stringb_ket,stringb_bra,&
       !$omp&        smo,bfac)
       !
       ! Loop over unique ket state beta strings
       do i=1,nb_ket
          
          ! Loop over unique bra state beta strings
          do j=1,nb_bra
             
             ! Matrix of orbital overlaps
             do k=1,nbeta
                moket=stringb_ket(k,i)
                do l=1,nbeta
                   mobra=stringb_bra(l,j)
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

  function isequivab(nalpha,nbeta)

    use constants
    use channels
!    use bdglobal
    
    implicit none

    integer, intent(in) :: nalpha,nbeta
    integer             :: i,j
    logical             :: isequivab

    ! Shortcut 1: The sets of unique alpha and beta strings are not
    ! equivalent if the number of elements in each is not equal
    if (na_bra.ne.nb_bra) then
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
    do i=1,na_ket
       do j=1,nalpha
          if (stringa_bra(i,j).ne.stringb_bra(i,j)) then
             isequivab=.false.
             exit
          endif
       enddo
    enddo
    
    return
    
  end function isequivab
  
!!#####################################################################
!! Calculate 1-particle overlap for all spin-orbitals in a given pair 
!! of determinants.
!! Taken from the superdyson code (see sd_core.f90).
!!######################################################################
!
!  subroutine fill_spinorbital_integrals(occbra,occket,adet,amo)
!
!    use constants
!    use bdglobal
!      
!    implicit none
!
!    integer, intent(in)   :: occbra(:) ! Parent determinant
!    integer, intent(in)   :: occket(:) ! Ion determinant
!    real(dp), intent(out) :: adet(:,:) ! Integral matrix
!                                       ! for a given determinant
!    real(dp), intent(in)  :: amo (:,:) ! Integral matrix for all MOs
!    
!    integer :: mobra, moket     ! Spatial MO indices
!    integer :: spinbra, spinket ! Spin+space MO indices
!    integer :: nalphabra        ! Number of spin-alpha electrons in the bra
!    
!    integer :: nmoket,nmobra
!
!    nmobra=nmo_bra
!    nmoket=nmo_ket
!    !
!    adet = 0.0d0
!    !
!    !  Alpha-alpha overlaps - upper left corner of sdet
!    !
!    spinket = 0
!    spinbra = 0
!    ket_alpha: do moket=1,nmoket
!       if (occket(moket)/=1 .and. occket(moket)/=2) cycle ket_alpha
!       !
!       spinket = spinket + 1
!       spinbra = 0
!       bra_alpha: do mobra=1,nmobra
!          if (occbra(mobra)/=1 .and. occbra(mobra)/=2) cycle bra_alpha
!          !
!          spinbra = spinbra + 1
!          adet(spinbra,spinket) = amo(mobra,moket)
!       end do bra_alpha
!    end do ket_alpha
!    nalphabra = spinbra
!    !
!    !  Beta-beta overlaps - bottom right corner of sdet
!    !
!    ket_beta: do moket=1,nmoket
!       if (occket(moket)/=-1 .and. occket(moket)/=2) cycle ket_beta
!       !
!       spinket = spinket + 1
!       spinbra = nalphabra
!       bra_beta: do mobra=1,nmobra
!          if (occbra(mobra)/=-1 .and. occbra(mobra)/=2) cycle bra_beta
!          !
!          spinbra = spinbra + 1
!          adet(spinbra,spinket) = amo(mobra,moket)
!       end do bra_beta
!    end do ket_beta
!    
!    return
!      
!  end subroutine fill_spinorbital_integrals
!  
!!######################################################################
!! determinant_overlap: Calculates the determinant of the matrix sdet of
!!                      overlaps between the MOs in two determinants.
!!                      Taken from superdyson.
!!######################################################################
!
!  function determinant_overlap(sdet) result(overdet)
!      
!    use constants
!    use bdglobal, only: nel
!    
!    implicit none
!    
!    real(dp), intent(in) :: sdet(:,:)
!    real(dp)             :: overdet
!    real(dp)             :: scr(nel,nel) ! Buffer for the 1-particle 
!                                         !  overlap matrix 
!    scr     = sdet
!    overdet = determinant(scr)
!    
!    return
!
!  end function determinant_overlap
!    
!!######################################################################
!! determinant: Calculates the determinant of a sqaure matrix using
!!              Gaussian elimination.
!!######################################################################
!
!  real(dp) function determinant(mat)
!
!    use constants
!    
!    real(dp), intent(inout) :: mat(:,:) ! Matrix destroyed on exit
!    integer                 :: order,info
!    integer                 :: ipvt(size(mat,dim=1))
!    integer                 :: i
!    real(dp)                :: work(size(mat,dim=1))
!    real(dp)                :: detx(2)
!    
!    order = size(mat,dim=1)
!
!    call dgefa(mat,order,order,ipvt,info)    
!
!    call dgedi(mat,order,order,ipvt,detx,work,10)
!
!    determinant = detx(1) * 10.0d0**detx(2)
!    
!    return
!
!  end function determinant

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

end module wfoverlaps

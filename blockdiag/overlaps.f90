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
    use bdglobal
    use import_gamess    

    implicit none

    integer               :: i,j,nao_ref,nao_disp
    real(dp), allocatable :: sao(:,:),ao2mo_ref(:,:),ao2mo_disp(:,:)

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
    
    return
    
  end subroutine mo_overlaps

!######################################################################

  subroutine psi_overlaps

    use constants
    use channels
    use bdglobal
    use omp_lib
    
    implicit none

    integer               :: i,j,m,k,nthreads,tid
    real(dp), allocatable :: smk(:,:)
    real(dp), allocatable :: spsi_1thread(:,:,:)
    real(dp)              :: detsmk
    real(dp), parameter   :: ovrthrsh=0.5d0
    logical               :: lovrlp
    
!----------------------------------------------------------------------
! Output what we are doing
!----------------------------------------------------------------------
    write(ilog,'(/,82a)') ('+',i=1,82)
    write(ilog,'(2x,a)') 'Calculating Adiabatic Wavefunction Overlaps'
    write(ilog,'(82a)') ('+',i=1,82)
    
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
! Get the alpha and beta spinorbital indices for every determinant
!----------------------------------------------------------------------
  call alpha_beta_indices

!----------------------------------------------------------------------
! Generate an integer label for every unique alpha and beta string
!----------------------------------------------------------------------
  call alpha_beta_labels

!----------------------------------------------------------------------
! Sort the alpha and beta strings
!----------------------------------------------------------------------
  call alpha_beta_sort
  
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
        
    return
    
  end subroutine psi_overlaps

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
    allocate(ia_ref(nalpha,maxdet,nsta))
    allocate(ib_ref(nbeta,maxdet,nsta))
    allocate(ia_disp(nalpha,maxdet,nsta))
    allocate(ib_disp(nbeta,maxdet,nsta))
    ia_ref=0
    ib_ref=0
    ia_disp=0
    ib_disp=0

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
                ia_ref(na,k,i)=imo
                ib_ref(nb,k,i)=imo
             else if (det_ref(imo,k,i).eq.+1) then
                na=na+1
                ia_ref(na,k,i)=imo
             else if (det_ref(imo,k,i).eq.-1) then
                nb=nb+1
                ib_ref(nb,k,i)=imo
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
                ia_disp(na,k,i)=imo
                ib_disp(nb,k,i)=imo
             else if (det_disp(imo,k,i).eq.+1) then
                na=na+1
                ia_disp(na,k,i)=imo
             else if (det_disp(imo,k,i).eq.-1) then
                nb=nb+1
                ib_disp(nb,k,i)=imo
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
    write(fmat_alpha,'(a,i0,a)') '(',nalpha,'i3)'
    write(fmat_beta,'(a,i0,a)') '(',nbeta,'i3)'
    
    ! Ref. states
    !
    do i=1,nsta
       do k=1,ndet_ref(i)

          ! Alpha spinorbital character string
          string_alpha=''
          write(string_alpha,fmat_alpha) ia_ref(:,k,i)

          ! Beta spinorbital character string
          string_beta=''
          write(string_beta,fmat_beta) ib_ref(:,k,i)
          
          ! Hashes of the alpha and beta spinorbital character string
          ilbla_ref(k,i)=djb_hash(string_alpha)
          ilblb_ref(k,i)=djb_hash(string_beta)
          
       enddo
    enddo

    ! Disp. states
    !
    do i=1,nsta
       do k=1,ndet_disp(i)

          ! Alpha spinorbital character string
          string_alpha=''
          write(string_alpha,fmat_alpha) ia_disp(:,k,i)

          ! Beta spinorbital character string
          string_beta=''
          write(string_beta,fmat_beta) ib_disp(:,k,i)
          
          ! Hashes of the alpha and beta spinorbital character string
          ilbla_disp(k,i)=djb_hash(string_alpha)
          ilblb_disp(k,i)=djb_hash(string_beta)
          
       enddo
    enddo
    
    return
    
  end subroutine alpha_beta_labels

!#####################################################################
  
  function djb_hash(str) result(hash)

    implicit none
    
    character(len=*),intent(in) :: str
    integer                     :: hash
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

    integer              :: i,k,sumdet,i1,i2,last,n,istate,idet
    integer, allocatable :: indx(:),ilbl(:)
    integer, allocatable :: info(:,:)

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    ! Indices of the unique alpha and beta strings for every
    ! determinant in the ref. and disp. states
    allocate(ias_ref(maxdet,nsta))
    allocate(ibs_ref(maxdet,nsta))
    allocate(ias_disp(maxdet,nsta))
    allocate(ibs_disp(maxdet,nsta))
    ias_ref=0
    ibs_ref=0
    ias_disp=0
    ibs_disp=0
    
!----------------------------------------------------------------------
! Determine the unique alpha and beta strings and the
! determinant/state pairs that they correspond to for the ref. states
! and the disp. states
!----------------------------------------------------------------------
    ! alpha, ref.
    call get_unique_strings(ndet_ref,ilbla_ref,na_ref,as_ref,ia_ref,&
         ias_ref,nalpha)

    ! beta, ref.
    call get_unique_strings(ndet_ref,ilblb_ref,nb_ref,bs_ref,ib_ref,&
         ibs_ref,nbeta)
    
    ! alpha, disp.
    call get_unique_strings(ndet_disp,ilbla_disp,na_disp,as_disp,&
         ia_disp,ias_disp,nalpha)
    
    ! beta, disp.
    call get_unique_strings(ndet_disp,ilblb_disp,nb_disp,bs_disp,&
         ib_disp,ibs_disp,nbeta)

    return

  end subroutine alpha_beta_sort

!#####################################################################

  subroutine get_unique_strings(ndet,ilblab,nab,abs,iab,iabs,nspin)

    use constants
    use channels
    use utils
    use bdglobal
    
    implicit none

    integer, dimension(nsta)            :: ndet
    integer, dimension(maxdet,nsta)     :: ilblab
    integer                             :: nab
    integer, allocatable                :: abs(:,:)
    integer, dimension(nab,maxdet,nsta) :: iab
    integer, dimension(maxdet,nsta)     :: iabs
    integer                             :: nspin

    integer              :: i,k,sumdet,i1,i2,last,n,istate,idet
    integer, allocatable :: indx(:),ilbl(:)
    integer, allocatable :: info(:,:)

!----------------------------------------------------------------------
! Allocate and initialise arrays
!----------------------------------------------------------------------
    sumdet=sum(ndet(:))
    allocate(indx(sumdet))
    indx=0
    allocate(info(sumdet,2))
    info=0
    allocate(ilbl(sumdet))
    ilbl=0

!----------------------------------------------------------------------
! Put together the combined list of alpha/beta string labels for the
! current set of states
!----------------------------------------------------------------------
    i1=0
    do i=1,nsta
       ! Next lot of alpha/beta strings
       i2=sum(ndet(1:i))
       i1=i2-ndet(i)+1
       ilbl(i1:i2)=ilblab(1:ndet(i),i)
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
    call isortindxa1('A',sumdet,ilbl,indx)
    
!----------------------------------------------------------------------  
! No. unique alpha/beta strings for the current set of states
!----------------------------------------------------------------------  
    last=0
    nab=0
    do k=1,sumdet
       if (ilbl(indx(k)).ne.last) then
          last=ilbl(indx(k))
          nab=nab+1
       endif
    enddo

!----------------------------------------------------------------------  
! Assign working alpha string indices to all of the determinants
! for all ref. states
!----------------------------------------------------------------------  
    allocate(abs(nspin,na_ref))

    last=0
    n=0
    do k=1,sumdet

       ! Determinant index for the current state/determinant pair
       idet=info(k,1)
       
       ! State index for the current state/determinant pair
       istate=info(k,2)

       ! Are we at the start of a new unique alpha string?
       if (ilbl(indx(k)).ne.last) then
          ! If so, increment the unique alpha string counter and
          ! fill in the nth unique alpha string
          last=ilbl(indx(k))
          n=n+1
          ! nth unique alpha string
          abs(:,n)=iab(:,idet,istate)
       endif
       
       ! Index of the unique alpha string for the current
       ! state/determinant pair
       iabs(idet,istate)=n

    enddo
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(indx)
    deallocate(info)
    deallocate(ilbl)
    
    return
    
  end subroutine get_unique_strings
    
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
! determinant: Calculates the determinant of a sqaure matrix.
!              Taken from superdyson.
!######################################################################

  real(dp) function determinant(mat)

    use constants
    
    real(dp), intent(inout) :: mat(:,:) ! Matrix destroyed on exit
    integer                 :: order, info
    integer                 :: ipvt(size(mat,dim=1))
    real(dp)                :: work(size(mat,dim=1))
    real(dp)                :: detx(2)
      
    order = size(mat,dim=1)
    if (order==0) stop 'null matrix passed to determinant'
    
    call dgefa(mat,order,order,ipvt,info)
    
    call dgedi(mat,order,order,ipvt,detx,work,10)
    
    determinant = detx(1) * 10.0d0**detx(2)
    
    return

  end function determinant
  
!######################################################################
  
end module overlaps

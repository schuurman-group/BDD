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

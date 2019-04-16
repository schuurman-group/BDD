!######################################################################
! mooverlaps: routines used in the calculation of the overlaps of
!             MOs at two different geometries
!######################################################################

module mooverlaps

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
    allocate(sao(nao_disp,nao_ref))
    sao=0.0d0
    
    ! MO overlap matrix
    allocate(smo(nmo_disp,nmo_ref))
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
  
end module mooverlaps

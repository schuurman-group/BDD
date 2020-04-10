!######################################################################
! Projection diabatisation using the algorithm of Tamura. In the limit
! of equal numbers of adiabats and diabats, this method is identical
! to propagative block diagonalisation diabatisation. However, we
! may also use more adiabats than diabats, which might be advantageous
! in some situations
!######################################################################
module tamura

contains

!######################################################################

  subroutine tamura_diabatisation
    
    use constants
    use channels
    use bdglobal
    use adtmod
    
    implicit none

    real(dp), allocatable :: protocoeff(:,:)
    
!----------------------------------------------------------------------
! Set the adiabatic potential and dipole matrices
!----------------------------------------------------------------------
    Vmat=Vmat1
    adip=adip1

!----------------------------------------------------------------------
! Calculate the overlaps between the adiabatic electronic states at the
! reference and displaced geometries
!----------------------------------------------------------------------
    call get_overlaps_tamura

!----------------------------------------------------------------------
! Write the wavefunction overlaps to the log file
!----------------------------------------------------------------------
    call wroverlaps_tamura

!----------------------------------------------------------------------
! Rephasing of the disp. states
!----------------------------------------------------------------------
    call rephase_tamura

!----------------------------------------------------------------------
! Calculation of the prototype functions as the overlaps between the
!  reference geometry diabats and the displaced geometry adiabats
!----------------------------------------------------------------------
    allocate(protocoeff(nsta_adiab,nsta_diab))
    call prototype_coefficients(protocoeff)
    
!----------------------------------------------------------------------
! Calculation of the projection diabatisation ADT matrix
!----------------------------------------------------------------------
    call adt_tamura(protocoeff)
    
    return

  end subroutine tamura_diabatisation

!######################################################################

  subroutine get_overlaps_tamura

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
    call psi_overlaps(spsi,nsta_disp,nsta_ref,nalpha,nbeta,&
         ndet_disp,ndet_ref,nmo_disp,nmo_ref,maxdet,c_disp,c_ref,&
         iocca_disp,iocca_ref,ioccb_disp,ioccb_ref,ioverlap,smo,&
         dthresh)

    return
    
  end subroutine get_overlaps_tamura
    
!######################################################################

    subroutine wroverlaps_tamura
    
    use constants
    use channels
    use bdglobal
    
    implicit none

    integer :: i,j

!----------------------------------------------------------------------
! Print out the wavefunction overlaps
!----------------------------------------------------------------------
    ! Table header
    write(ilog,'(/,47a)') ('-',i=1,47)
    write(ilog,'(a)') '  Disp. State  |  Ref. State  |  Overlap'
    write(ilog,'(a)') '     |I>       |    |J>       |   <I|J>'
    write(ilog,'(47a)') ('-',i=1,47)

    ! Table entries
    do i=1,nsta_disp
       do j=1,nsta_ref
          write(ilog,'(5x,i2,13x,i2,11x,F13.10)') i,j,spsi(i,j)
       enddo
    enddo
       
    ! End of the table
    write(ilog,'(47a)') ('-',i=1,47)
    
    return
    
  end subroutine wroverlaps_tamura

!######################################################################

    subroutine rephase_tamura

    use constants
    use channels
    use iomod
    use bdglobal
    
    implicit none

    integer                        :: i,j,k,n,idet
    real(dp), dimension(nsta_disp) :: phfac
    character(len=128)             :: fmt

!-----------------------------------------------------------------------
! Try to determine if the phase of a wavefunction has switched from the
! ref. geometry
!-----------------------------------------------------------------------    
    phfac=1.0d0
    do i=1,nsta_disp
       do j=1,nsta_ref
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
    do i=1,nsta_disp
       spsi(i,:)=phfac(i)*spsi(i,:)
    enddo

!-----------------------------------------------------------------------
! Re-phase the adiabatic dipole matrix elements
!-----------------------------------------------------------------------
    if (ldipole) then
       do i=1,nsta_adiab
          do j=1,nsta_adiab
             adip(i,j,:)=adip(i,j,:)*phfac(i)*phfac(j)
          enddo
       enddo
    endif
    
!-----------------------------------------------------------------------
! Write the phase factors to the log file to be used at the next
! geometry
!-----------------------------------------------------------------------
    write(ilog,'(/,82a)') ('+',i=1,82)
    write(ilog,'(2x,a)') 'Phase Factors'
    write(ilog,'(82a)') ('+',i=1,82)
    do i=1,nsta
       write(ilog,'(2x,i2,2x,i2)') i,int(phfac(i))
    enddo

    return

  end subroutine rephase_tamura
    
!######################################################################

  subroutine prototype_coefficients(protocoeff)

    use constants
    use iomod
    use bdglobal
    
    implicit none

    integer               :: i,j    
    real(dp), intent(out) :: protocoeff(nsta_adiab,nsta_diab)

!----------------------------------------------------------------------
! If we are not reading in the ADT matrix form the reference geometry
! then we will assume that it is the unit matrix
!----------------------------------------------------------------------
    if (.not.lreftrans) then
       protocoeff=spsi
       return
    endif
       
!----------------------------------------------------------------------
! If needed, read the transformation matrix from a previous log file
!----------------------------------------------------------------------
    if (lrdreftrans) call rdreftrans_tamura

!-----------------------------------------------------------------------
! MGS orthonormalisation of the transformation to get rid of any
! issuses arising from the finite precision input of the transformation
! matrix
!-----------------------------------------------------------------------
    do i=1,nsta_diab
       do j=1,i-1
          reftrans(:,i)=reftrans(:,i)&
               -dot_product(reftrans(:,i),reftrans(:,j))*reftrans(:,j)
       enddo
       reftrans(:,i)=reftrans(:,i)&
            /sqrt(dot_product(reftrans(:,i),reftrans(:,i)))
    enddo

!-----------------------------------------------------------------------
! Get the prototype function expansion coefficients
!-----------------------------------------------------------------------
    protocoeff=matmul(spsi,reftrans)

    return
    
  end subroutine prototype_coefficients

!######################################################################

  subroutine rdreftrans_tamura

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

    do i=1,nsta_adiab
       do j=1,nsta_diab
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
    
  end subroutine rdreftrans_tamura
  
!######################################################################

  subroutine adt_tamura(protocoeff)

    use constants
    use channels
    use iomod
    use utils
    use bdglobal
    
    implicit none

    real(dp), intent(in) :: protocoeff(nsta_adiab,nsta_diab)
    real(dp)             :: coecoeT(nsta_adiab,nsta_adiab)
    real(dp)             :: inv_coecoeT(nsta_adiab,nsta_adiab)
    real(dp)             :: invsqrt_coecoeT(nsta_adiab,nsta_adiab)
    logical              :: pseudo
    
    coecoeT=matmul(protocoeff,transpose(protocoeff))

    call invert_matrix(coecoeT,inv_coecoeT,nsta_adiab,pseudo)

    if (pseudo) write(ilog,'(/,a)') 'WARNING: Pseudoinverse used in &
         the construction of the ADT matrix!'

    call sqrt_matrix(inv_coecoeT,invsqrt_coecoeT,nsta_adiab)

    adt=matmul(invsqrt_coecoeT,protocoeff)
    
    return
    
  end subroutine adt_tamura
  
!######################################################################
  
end module tamura

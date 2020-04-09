!######################################################################
! Propagative block diagonalisation using the original algorithm of
! Pacher, Koppel and Cederbaum. That is, using equal numbers of
! adiabatic and diabatic states.
!######################################################################
module pacher

contains

!######################################################################

  subroutine pacher_diabatisation
    
    use constants
    use channels
    use bdglobal
    use adtmod
    
    implicit none

!----------------------------------------------------------------------
! Determine the disp. geometry states to use by trying to follow
! electronic characters from geometry to geometry
!----------------------------------------------------------------------
    if (nsta_disp.ne.nsta_ref) then
       ! More disp. states than ref. states: select the subset of disp.
       ! states that map onto the ref. states
       call trackwfs
    else
       ! Equal numbers of disp. and ref. states: set the adiabatic
       ! potential and dipole matrices here
       Vmat=Vmat1
       adip=adip1
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
! Calculate the ADT matrix
!----------------------------------------------------------------------
    call get_adt

!----------------------------------------------------------------------
! Check on the swapping of diabats. Note that this can only be done if
! we have access to the ADT matrix of the previous geometry
!----------------------------------------------------------------------
    if (lreftrans) call switch_diabats

!----------------------------------------------------------------------
! Write the adiabatic potentials to the log file
!----------------------------------------------------------------------
    call write_adiabpot
  
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

!----------------------------------------------------------------------
! Optional: write the ADT matrix to file in a format compatible with
!           the DFT/MRCI code
!----------------------------------------------------------------------
    if (ldmat) call write_dmat_trans
      
    return
  
  end subroutine pacher_diabatisation
    
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
    real(dp), allocatable           :: tmp(:,:),sumsq(:)
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

    ! Overlap analysis arrays
    allocate(tmp(nsta_disp,nsta_ref))
    allocate(sumsq(nsta_ref))
    tmp=0.0d0
    sumsq=0.0d0
    
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
! Determine which disp. states are needed to represent the ref. states
!----------------------------------------------------------------------
    allocate(indx(nsta_ref))
    indx=0
    
    tmp=spsi1
    do i=1,nsta_ref
       indx(i)=maxloc(abs(tmp(:,i)),dim=1)
       tmp(indx(i),:)=0.0d0
    enddo
    
!----------------------------------------------------------------------
! Write the selected state indices to the log file
!----------------------------------------------------------------------
    write(ilog,'(a)') ''
    do i=1,nsta_ref
       write(ilog,'(x,2(x,a,x,i0))') 'Selected state',i,':',indx(i)
    enddo

!----------------------------------------------------------------------
! Write the sums of the squares of the overlaps of the selected disp.
! states with the ref. states to the log file
!----------------------------------------------------------------------
    sumsq=0.0d0
    do i=1,nsta_ref
       do j=1,nsta_ref
          sumsq(i)=sumsq(i)+spsi1(indx(j),i)**2
       enddo
    enddo
       
    write(ilog,'(/,47a)') ('-',i=1,47)
    write(ilog,'(a)') '  Ref. State  |  Sum of Squared Overlaps'
    write(ilog,'(a)') '    |J>       |      Sum_n <I_n|J>**2'
    write(ilog,'(47a)') ('-',i=1,47)
    do i=1,nsta_ref
       write(ilog,'(5x,i2,17x,F13.10)'),i,sumsq(i)
    enddo
    write(ilog,'(47a)') ('-',i=1,47)
    
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
! Make a copy of the indices of the selected disp. states
!----------------------------------------------------------------------
    allocate(isel(nsta_ref))
    isel=indx
    
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
    deallocate(tmp)
    deallocate(sumsq)
    
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

!----------------------------------------------------------------------
! Rephasing of the disp. states
!----------------------------------------------------------------------
    call rephase

!----------------------------------------------------------------------
! Optional transformation of the reference geometry wavefunctions
! (via the transformation of the wavefunction overlap matrix)
!----------------------------------------------------------------------
  if (lreftrans) call trans_refpsi
    
    return
    
  end subroutine wroverlaps

!######################################################################

  subroutine rephase

    use constants
    use channels
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
! Re-phase the adiabatic dipole matrix elements
!-----------------------------------------------------------------------
    if (ldipole) then
       do i=1,nsta
          do j=1,nsta
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

  subroutine write_adiabpot

    use constants
    use channels
    use iomod
    use bdglobal
    
    implicit none

    integer :: i

!----------------------------------------------------------------------
! Write the energies of the selected adiabatic states to the log file
!----------------------------------------------------------------------
    write(ilog,'(/,82a)') ('+',i=1,82)
    write(ilog,'(2x,a)') 'Adiabatic Potentials (Selected States Only)'
    write(ilog,'(82a)') ('+',i=1,82)

    do i=1,nsta
       write(ilog,'(2x,i2,2x,F15.10)') i,vmat(i,i)
    enddo
    
    return
    
  end subroutine write_adiabpot

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

  subroutine write_dmat_trans

    use constants
    use channels
    use iomod
    use bdglobal
    
    implicit none

    integer            :: ilbl,unit,i,j
    character(len=120) :: admat
    
!----------------------------------------------------------------------
! Open the dmat transformation matrix file
!----------------------------------------------------------------------
    ! Filename
    ilbl=index(ain,'.inp')
    admat=ain(1:ilbl-1)//'.trans'

    ! Open the dmat transformation matrix file
    call freeunit(unit)
    open(unit,file=admat,form='formatted',status='unknown')

!----------------------------------------------------------------------
! Write the dmat transformation matrix file
! Note that here we must use the indices of the selected adiabats    
!----------------------------------------------------------------------
    do i=1,nsta
       do j=1,nsta
          write(unit,'(2(2x,i2),2x,F15.10)') isel(i),isel(j),adt(i,j)
       enddo
    enddo
    
!----------------------------------------------------------------------
! Close the dmat transformation matrix file
!----------------------------------------------------------------------
    close(unit)
    
    return
    
  end subroutine write_dmat_trans
  
!######################################################################
  
end module pacher

!######################################################################
! detparsing: routines for the parsing of the determinant lists
!######################################################################

module detparsing

  implicit none

contains
  
!######################################################################
  
  subroutine rddetfiles

    use constants
    use channels
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
! Optional truncation of the wavefunctions by norm
!-----------------------------------------------------------------------
    if (ltruncate) then

       ! Sort the coefficient and spinorbital index arrays in order of
       ! decreasing absolute coefficient value
       call sort_det_arrays

       ! Truncation of the wavefunctions
       call truncate_det_arrays
       
    endif
    
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


       !do k=1,ndet_ref(i)
       !   call rdinp(idet)
       !   read(keyword(1),*) c_ref(k,i)
       !   do n=2,inkw
       !      read(keyword(n),*) det_ref(n-1,k,i)
       !   enddo          
       !enddo
       
       do k=1,ndet_ref(i)
          read(idet,*) c_ref(k,i),(det_ref(n,k,i),n=1,nmo_ref)
       enddo
          
       close(idet)
    enddo

    ! Displaced geometry
    do i=1,nsta_disp
       open(idet,file=adetdisp(i),form='formatted',status='old')

       !do k=1,ndet_disp(i)
       !   call rdinp(idet)
       !   read(keyword(1),*) c_disp(k,i)          
       !   do n=2,inkw
       !      read(keyword(n),*) det_disp(n-1,k,i)
       !   enddo          
       !enddo

       do k=1,ndet_disp(i)
          read(idet,*) c_disp(k,i),(det_disp(n,k,i),n=1,nmo_ref)
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

  subroutine sort_det_arrays

    use constants
    use bdglobal
    use utils
    
    implicit none

    integer               :: i,k
    integer, allocatable  :: indx(:),iswapvec(:,:)
    real(dp), allocatable :: cabs(:),fswapvec(:)
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(indx(maxdet))
    allocate(cabs(maxdet))
    allocate(fswapvec(maxdet))
    allocate(iswapvec(max(nalpha,nbeta),maxdet))
    
!----------------------------------------------------------------------
! Sort the coefficient and spin orbital indec arrays in order of
! decreasing absolute value
!----------------------------------------------------------------------
    ! Disp. geometry wavefunctions
    do i=1,nsta_disp

       ! Skip if the current wavefunction has a norm less than
       ! the cutoff value
       if (norm_disp(i).le.normcut) cycle
       
       ! Sort the absolute coefficient values
       indx=0
       cabs(:)=abs(c_disp(:,i))
       call dsortindxa1('D',ndet_disp(i),cabs(1:ndet_disp(i)),&
            indx(1:ndet_disp(i)))

       ! Rearrange the coefficient vector
       fswapvec=0.0d0
       do k=1,ndet_disp(i)
          fswapvec(k)=c_disp(indx(k),i)
       enddo
       c_disp(:,i)=fswapvec(:)

       ! Rearrange the alpha spin orbital index array
       iswapvec=0
       do k=1,ndet_disp(i)
          iswapvec(1:nalpha,k)=iocca_disp(:,indx(k),i)
       enddo
       iocca_disp(:,:,i)=iswapvec(1:nalpha,:)

       ! Rearrange the beta spin orbital index array
       iswapvec=0
       do k=1,ndet_disp(i)
          iswapvec(1:nbeta,k)=ioccb_disp(:,indx(k),i)
       enddo
       ioccb_disp(:,:,i)=iswapvec(1:nbeta,:)
       
    enddo

    ! Ref. geometry wavefunctions
    do i=1,nsta_ref

       ! Skip if the current wavefunction has a norm less than
       ! the cutoff value
       if (norm_ref(i).le.normcut) cycle
       
       ! Sort the absolute coefficient values
       indx=0
       cabs(:)=abs(c_ref(:,i))
       call dsortindxa1('D',ndet_ref(i),cabs(1:ndet_ref(i)),&
            indx(1:ndet_ref(i)))

       ! Rearrange the coefficient vector
       fswapvec=0.0d0
       do k=1,ndet_ref(i)
          fswapvec(k)=c_ref(indx(k),i)
       enddo
       c_ref(:,i)=fswapvec(:)

       ! Rearrange the alpha spin orbital index array
       iswapvec=0
       do k=1,ndet_ref(i)
          iswapvec(1:nalpha,k)=iocca_ref(:,indx(k),i)
       enddo
       iocca_ref(:,:,i)=iswapvec(1:nalpha,:)

       ! Rearrange the beta spin orbital index array
       iswapvec=0
       do k=1,ndet_ref(i)
          iswapvec(1:nbeta,k)=ioccb_ref(:,indx(k),i)
       enddo
       ioccb_ref(:,:,i)=iswapvec(1:nbeta,:)
       
    enddo
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(indx)
    deallocate(cabs)
    deallocate(fswapvec)
    deallocate(iswapvec)
    
    return
    
  end subroutine sort_det_arrays

!######################################################################

  subroutine truncate_det_arrays

    use constants
    use channels
    use bdglobal
    
    implicit none

    integer              :: i,k,ndet
    integer, allocatable :: ndet0_disp(:),ndet0_ref(:)
    real(dp)             :: normsq

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(ndet0_disp(nsta_disp))
    allocate(ndet0_ref(nsta_ref))
    ndet0_disp=0
    ndet0_ref=0

!----------------------------------------------------------------------
! Save the old numbers of determinants
!----------------------------------------------------------------------
    ndet0_disp=ndet_disp
    ndet0_ref=ndet_ref
    
!----------------------------------------------------------------------
! Re-set the ndet_disp and ndet_ref arrays
!----------------------------------------------------------------------
    ! Disp. geometry wavefunctions
    do i=1,nsta_disp
       
       ! Skip if the current wavefunction has a norm less than
       ! the cutoff value
       if (norm_disp(i).le.normcut) cycle

       ! Get the truncated no. determinants
       ndet=0
       normsq=0.0d0
       do k=1,ndet_disp(i)
          normsq=normsq+c_disp(k,i)**2
          if (sqrt(normsq).gt.normcut) then
             ndet=k
             exit
          endif
       enddo
       ndet_disp(i)=ndet
       
    enddo

    ! Ref. geometry wavefunctions
    do i=1,nsta_ref
       
       ! Skip if the current wavefunction has a norm less than
       ! the cutoff value
       if (norm_ref(i).le.normcut) cycle

       ! Get the truncated no. determinants
       ndet=0
       normsq=0.0d0
       do k=1,ndet_ref(i)
          normsq=normsq+c_ref(k,i)**2
          if (sqrt(normsq).gt.normcut) then
             ndet=k
             exit
          endif
       enddo
       ndet_ref(i)=ndet

    enddo

!----------------------------------------------------------------------
! Be cautious and zero the elements of the coefficient and spin
! orbital index arrays that we are discarding
!----------------------------------------------------------------------
    ! Disp. geometry wavefunctions
    do i=1,nsta_disp
       if (norm_disp(i).le.normcut) cycle
       c_disp(ndet_disp(i)+1:,i)=0.0d0
       iocca_disp(:,ndet_disp(i)+1,i)=0
       ioccb_disp(:,ndet_disp(i)+1,i)=0
    enddo

    ! Ref. geometry wavefunctions
    do i=1,nsta_ref
       if (norm_ref(i).le.normcut) cycle
       c_ref(ndet_ref(i)+1:,i)=0.0d0
       iocca_ref(:,ndet_ref(i)+1,i)=0
       ioccb_ref(:,ndet_ref(i)+1,i)=0
    enddo

!----------------------------------------------------------------------
! Reset the norms
!----------------------------------------------------------------------
    ! Disp. geometry wavefunctions
    do i=1,nsta_disp
       norm_disp(i)=sqrt(sum(c_disp(:,i)**2))
    enddo

    ! Ref. geometry wavefunctions
    do i=1,nsta_ref
       norm_ref(i)=sqrt(sum(c_ref(:,i)**2))
    enddo
    
!----------------------------------------------------------------------
! Output the truncation information to the log file
!----------------------------------------------------------------------
    ! Ref. geometry wavefunctions
    write(ilog,'(/,47a)') ('-',i=1,47)
    write(ilog,'(a)') '   Ref. State  |     NDet    |     NDet'
    write(ilog,'(a)') '     |I>       |   Original  |   Truncated'
    write(ilog,'(47a)') ('-',i=1,47)
    do i=1,nsta_ref
       write(ilog,'(5x,i2,8x,a1,3x,i8,2x,a1,3x,i8)') &
            i,'|',ndet0_ref(i),'|',ndet_ref(i)
    enddo
    write(ilog,'(47a)') ('-',i=1,47)

    ! Disp. geometry wavefunctions
    write(ilog,'(/,47a)') ('-',i=1,47)
    write(ilog,'(a)') '  Disp. State  |     NDet    |     NDet'
    write(ilog,'(a)') '     |I>       |   Original  |   Truncated'
    write(ilog,'(47a)') ('-',i=1,47)
    do i=1,nsta_disp
       write(ilog,'(5x,i2,8x,a1,3x,i8,2x,a1,3x,i8)') &
            i,'|',ndet0_disp(i),'|',ndet_disp(i)
    enddo
    write(ilog,'(47a)') ('-',i=1,47)
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(ndet0_disp)
    deallocate(ndet0_ref)
    
    return
    
  end subroutine truncate_det_arrays
    
!######################################################################
  
end module detparsing

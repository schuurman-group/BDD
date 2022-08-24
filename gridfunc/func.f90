module func

  implicit none

contains

!######################################################################

  subroutine calc_func

    use constants
    use utils
    use sysinfo
    use gridglobal
    
    implicit none

    integer(8)            :: ntotal,i
    integer               :: m,m1,s1,s2
    integer               :: indx(nfuncmode),iswap(nfuncmode)
    real(dp)              :: q(nfuncmode)
    real(dp), allocatable :: func(:,:,:)

!----------------------------------------------------------------------
! Sort the direct product (sub) grid modes to be in ascending order
!----------------------------------------------------------------------
    call i4sortindxa1('A',nfuncmode,funcmode,indx)

    do m=1,nfuncmode
       iswap(m)=funcmode(indx(m))
    enddo

    funcmode=iswap

!----------------------------------------------------------------------
! Total number of points on the (sub) direct product grid
!----------------------------------------------------------------------
    ntotal=1
    do m1=1,nfuncmode
       m=funcmode(m1)
       ntotal=ntotal*ndvr(m)
    enddo

!----------------------------------------------------------------------
! Make the output directory
!----------------------------------------------------------------------
    call mkoutdir
    
!----------------------------------------------------------------------
! Compute the function value at the grid points
!----------------------------------------------------------------------
    select case(ifunc)
       
    case(1) ! Projector onto an adiabatic state
       call calc_adproj(ntotal)
       
    case(2) ! Adiabatic state excitation operator
       call calc_adexci(ntotal)
       
    case(3) ! Rigid rotor energies
       call calc_roten(ntotal)
       
    end select
    
    return
    
  end subroutine calc_func

!######################################################################

  subroutine calc_adproj(ntotal)

    use constants
    use sysinfo
    use gridglobal

    implicit none

    integer(8), intent(in) :: ntotal
    integer(8)             :: i
    integer                :: s1,s2
    real(dp)               :: q(nfuncmode)
    real(dp), allocatable  :: func(:,:,:)

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(func(nsta,nsta,ntotal))
    func=0.0d0

!----------------------------------------------------------------------
! Compute the matrix representation of the function on the
! (sub) direct product grid
!----------------------------------------------------------------------
    ! Loop over grid points
    do i=1,ntotal

       ! Mode values at the current grid point
       call mode_values(i,q)

       ! Value of the elements of the diabatic state representation
       ! of the adiabatic state projector at the current grid point
       func(:,:,i)=adiabatic_projector(q)
       
    enddo

!----------------------------------------------------------------------
! Write the elements of the matrix representation of the function
! to disk
!----------------------------------------------------------------------
    ! Loop over pairs of states
    do s1=1,nsta
       do s2=s1,nsta

          ! Write the current element of the matrix representation
          ! of the function to disk
          call wrfunc_1element(s1,s2,func(s1,s2,:),ntotal)

       enddo
    enddo
       
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(func)
    
    return
    
  end subroutine calc_adproj

!######################################################################

  subroutine calc_adexci(ntotal)

    use constants
    use sysinfo
    use gridglobal
    
    implicit none

    integer(8), intent(in) :: ntotal
    integer(8)             :: i
    integer                :: s1,s2
    real(dp)               :: q(nfuncmode)
    real(dp), allocatable  :: func(:,:,:)
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(func(nsta,nsta,ntotal))
    func=0.0d0

!----------------------------------------------------------------------
! Compute the matrix representation of the function on the
! (sub) direct product grid
!----------------------------------------------------------------------
    ! Loop over grid points
    do i=1,ntotal

       ! Mode values at the current grid point
       call mode_values(i,q)

       ! Value of the elements of the diabatic state representation
       ! of the adiabatic state excitation operator at the current
       ! grid point
       func(:,:,i)=adiabatic_excitation(q)
       
    enddo

!----------------------------------------------------------------------
! Write the elements of the matrix representation of the function
! to disk
!----------------------------------------------------------------------
    ! Loop over pairs of states
    !do s1=1,nsta
    !   do s2=s1,nsta
    !
    !      ! Write the current element of the matrix representation
    !      ! of the function to disk
    !      call wrfunc_1element(s1,s2,func(s1,s2,:),ntotal)
    !
    !   enddo
    !enddo

    do s1=1,nsta
       do s2=1,nsta
    
          ! Write the current element of the matrix representation
          ! of the function to disk
          call wrfunc_1element(s1,s2,func(s1,s2,:),ntotal)
    
       enddo
    enddo
    
    return
    
  end subroutine calc_adexci

!######################################################################

  subroutine calc_roten(ntotal)

    use constants
    use sysinfo
    use gridglobal

    implicit none

    integer(8), intent(in) :: ntotal
    integer(8)             :: i
    integer                :: K
    real(dp)               :: q(nfuncmode)
    real(dp), allocatable  :: func(:,:)

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(func(2*Jval+1,ntotal))
    func=0.0d0

!----------------------------------------------------------------------
! Compute the matrix representation of the function on the
! (sub) direct product grid
!----------------------------------------------------------------------
    ! Loop over grid points
    do i=1,ntotal

       ! Mode values at the current grid point
       call mode_values(i,q)

       ! Value of the rigid rotor energy at the current grid point
       func(:,i)=rigid_rotor_energy(q)
       
    enddo

!----------------------------------------------------------------------
! Write the rigid rotor energies to disk
!----------------------------------------------------------------------
    ! Loop over K values
    do K=-Jval,Jval
       call wrfunc_1element(Jval,K,func(K+Jval+1,i),ntotal)
    enddo
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(func)
    
    return
    
  end subroutine calc_roten
  
!######################################################################

  subroutine mode_values(i,q)

    use constants
    use gridglobal
    
    implicit none

    integer(8), intent(in) :: i
    integer(8)             :: j,jj
    integer                :: m1,m,indx
    real(dp), intent(out)  :: q(nfuncmode)

!----------------------------------------------------------------------
! For the i'th point on the (sub) direct product grid, determine the
! mode values. Adapated from the decoder subroutine in the Quantics
! code.
!----------------------------------------------------------------------
    j=i-1
    do m1=1,nfuncmode
       m=funcmode(m1)
       jj=j/ndvr(m)
       indx=j-jj*ndvr(m)+1
       q(m1)=grid(indx,m)
       j=jj
    enddo
    
    return
    
  end subroutine mode_values
  
!######################################################################

  function adiabatic_projector(q) result(proj)

    use constants
    use sysinfo
    use potfuncs
    use extfunc
    use gridglobal
    
    implicit none

    integer              :: m,m1,i,j,iproj
    real(dp), intent(in) :: q(nfuncmode)
    real(dp)             :: proj(nsta,nsta)
    real(dp)             :: q1(nmodes)
    real(dp)             :: adt(nsta,nsta)
        
!----------------------------------------------------------------------
! Normal mode coordinates in the full space
!----------------------------------------------------------------------
    q1=0.0d0
    do m1=1,nfuncmode
       m=funcmode(m1)
       q1(m)=q(m1)
    enddo

!----------------------------------------------------------------------
! Compute the ADT matrix at the grid point q
!----------------------------------------------------------------------    
    if (idiabfunc == 1) then
       adt=adtmatrix(q1)
    else
       adt=extfunc_adtmatrix(q1)
    endif

!----------------------------------------------------------------------
! Diabatic representation of the projector onto the adiabatic state
! of interest
!----------------------------------------------------------------------
    ! Adiabatic state index
    iproj=funcsta(1)

    ! Form the matrix representation of the projector
    do i=1,nsta
       do j=1,nsta
          proj(j,i)=adt(j,iproj)*adt(i,iproj)
       enddo
    enddo

    return
    
  end function adiabatic_projector

!######################################################################

  function adiabatic_excitation(q) result(proj)

    use constants
    use sysinfo
    use potfuncs
    use extfunc
    use gridglobal
    
    implicit none

    integer              :: m,m1,n,i,j,iexci1,iexci2
    real(dp), intent(in) :: q(nfuncmode)
    real(dp)             :: proj(nsta,nsta)
    real(dp)             :: q1(nmodes)
    real(dp)             :: adt(nsta,nsta)
        
!----------------------------------------------------------------------
! Normal mode coordinates in the full space
!----------------------------------------------------------------------
    q1=0.0d0
    do m1=1,nfuncmode
       m=funcmode(m1)
       q1(m)=q(m1)
    enddo

!----------------------------------------------------------------------
! Compute the ADT matrix at the grid point q
!----------------------------------------------------------------------    
    if (idiabfunc == 1) then
       adt=adtmatrix(q1)
    else
       adt=extfunc_adtmatrix(q1)
    endif

!----------------------------------------------------------------------
! Diabatic representation of the adiabatic state excitation operator
!----------------------------------------------------------------------
    ! Adiabatic state indices
    iexci1=funcsta(1)
    iexci2=funcsta(2)
    
    ! Form the matrix representation of the excitation operator
    !do i=1,nsta
    !   do j=1,nsta
    !      proj(j,i)=adt(j,iexci1)*adt(i,iexci2) &
    !           +adt(i,iexci1)*adt(j,iexci2)
    !   enddo
    !enddo

    do i=1,nsta
       do j=1,nsta
          proj(j,i)=adt(j,iexci1)*adt(iexci2,i)
       enddo
    enddo
    
    return
    
  end function adiabatic_excitation
  
!######################################################################

  function rigid_rotor_energy(q) result(EJK)

    use constants
    use iomod
    use utils
    use sysinfo
    use gridglobal
    
    implicit none

    integer               :: m,m1
    integer               :: hdim
    real(dp), intent(in)  :: q(nfuncmode)
    real(dp)              :: EJK(2*Jval+1)
    real(dp)              :: q1(nmodes)
    real(dp)              :: x(ncoo)
    real(dp)              :: A,B,C
    real(dp), allocatable :: hmat(:,:)
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    hdim=2*Jval+1

    allocate(hmat(hdim,hdim))
    hmat=0.0d0
    
!----------------------------------------------------------------------
! Normal mode coordinates in the full space
!----------------------------------------------------------------------
    q1=0.0d0
    do m1=1,nfuncmode
       m=funcmode(m1)
       q1(m)=q(m1)
    enddo

!----------------------------------------------------------------------
! Cartesian coordinates in a.u.
!----------------------------------------------------------------------
    x=xcoo0+matmul(nmcoo,q1)*ang2bohr

!----------------------------------------------------------------------
! Rotational constants
!----------------------------------------------------------------------
    call rotational_constants(x,A,B,C)
    
!-----------------------------------------------------------------------
! Construct the rigid rotor Hamiltonian matrix
!-----------------------------------------------------------------------
    call rigid_rotor_hamiltonian(hmat,A,B,C,hdim)

!-----------------------------------------------------------------------
! Diagonalise the rigid rotor Hamiltonian matrix
!-----------------------------------------------------------------------
    call diag_matrix(hmat,EJK,hdim)

!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
    deallocate(hmat)
    
    return
    
  end function rigid_rotor_energy
    
!######################################################################

  subroutine rotational_constants(x1,A,B,C)

    use constants
    use iomod
    use sysinfo
    
    implicit none

    integer               :: i,j,error
    real(dp), intent(in)  :: x1(ncoo)
    real(dp), intent(out) :: A,B,C
    real(dp)              :: x(ncoo)
    real(dp)              :: itensor(3,3)
    real(dp)              :: iteig(3)
    real(dp)              :: work(9)
    real(dp)              :: totmass,com(3)
    real(dp)              :: xx,xy,xz,xm
    
!----------------------------------------------------------------------
! Shift to the centre of mass
!----------------------------------------------------------------------
    ! Total mass
    totmass=0.0d0
    do i=1,natm
       totmass=totmass+mass(i*3)
    enddo

    ! Center of mass
    com=0.0d0
    do i=1,natm
       do j=1,3
          com(j)=com(j)+mass(i*3)*x1(i*3-3+j)
       enddo
    enddo
    com=com/totmass

    ! Shift the origin to the centre of mass
    x=x1
    do i=1,natm
       do j=1,3
          x(i*3-3+j)=x(i*3-3+j)-com(j)
       enddo
    enddo
    
!----------------------------------------------------------------------
! Moment of intertia tensor in a.u.
!----------------------------------------------------------------------
    itensor=0.0d0
    do i=1,natm
       xx=x(i*3-2)
       xy=x(i*3-1)
       xz=x(i*3)
       xm=mass(i*3)
       itensor(1,1)=itensor(1,1)+xm*(xy**2+xz**2)
       itensor(2,2)=itensor(2,2)+xm*(xx**2+xz**2)
       itensor(3,3)=itensor(3,3)+xm*(xx**2+xy**2)
       itensor(1,2)=itensor(1,2)-xm*xx*xy
       itensor(1,3)=itensor(1,3)-xm*xx*xz
       itensor(2,3)=itensor(2,3)-xm*xy*xz
    enddo
    itensor(2,1)=itensor(1,2)
    itensor(3,1)=itensor(1,3)
    itensor(3,2)=itensor(2,3)

    ! Convert to atomic units of mass
    itensor=itensor*mu2me

!-----------------------------------------------------------------------
! Diagonalise the moment of inertia tensor
!-----------------------------------------------------------------------
    call dsyev('V','U',3,itensor,3,iteig,work,9,error)
    
    if (error.ne.0) then
       errmsg='Diagonalisation of the moment of inertia tensor in &
            rotational_energy failed'
       call error_control
    endif

!-----------------------------------------------------------------------
! Rotational constants in a.u.
!-----------------------------------------------------------------------
    A=1.0d0/(2.0d0*iteig(1))
    B=1.0d0/(2.0d0*iteig(2))
    C=1.0d0/(2.0d0*iteig(3))
    
    return
    
  end subroutine rotational_constants

!######################################################################

  subroutine rigid_rotor_hamiltonian(hmat,A,B,C,hdim)

    use constants
    use gridglobal
    
    implicit none

    integer, intent(in)   :: hdim
    integer               :: K,Kp,J,i,ip
    real(dp), intent(out) :: hmat(hdim,hdim)
    real(dp), intent(in)  :: A,B,C

!-----------------------------------------------------------------------
! Initialisation
!-----------------------------------------------------------------------
    hmat=0.0d0
    
!-----------------------------------------------------------------------
! Diagonal elements
!-----------------------------------------------------------------------
    J=Jval
    
    i=0
    do K=-J,J
       i=i+1
       hmat(i,i)=0.5d0*(A+B)*(J*(J+1)-K**2) + C*K**2
    enddo

!-----------------------------------------------------------------------
! Off-diagonal elements
!-----------------------------------------------------------------------
    do K=-J,J-1
       i=K+J+1
       do Kp=K+1,J
          ip=Kp+J+1
          if (Kp.eq.K+2) then
             hmat(i,ip)=(B-A)*sqrt(dble((J-K)*(J-K-1)*(J+K+1)*(J+K+2)))
             hmat(i,ip)=0.25d0*hmat(i,ip)
             hmat(ip,i)=hmat(i,ip)
          else if (Kp.eq.K-2) then
             hmat(i,ip)=(B-A)*sqrt(dble((J+K)*(J+K-1)*(J-K+1)*(J-K+2)))
             hmat(i,ip)=0.25d0*hmat(i,ip)
             hmat(ip,i)=hmat(i,ip)
          endif
       enddo
    enddo
    
    return
    
  end subroutine rigid_rotor_hamiltonian
  
!######################################################################
  
  subroutine mkoutdir

    use constants
    use gridglobal
    
    implicit none

    character(len=80) :: dirname

    ! Output directory name
    select case(ifunc)
    case(1) ! Projector onto an adiabatic state
       write(dirname,'(a,i0)') 'adproj',funcsta(1)
    case(2) ! Adiabatic state excitation operator
       write(dirname,'(a,2(i0))') 'adexci',funcsta(1),funcsta(2)
    case(3) ! Rotational energies
       write(dirname,'(a,i0)') 'EJ',Jval
    end select

    ! Clean up
    call system('rm -rf '//trim(dirname)//' 2>/dev/null')

    ! Make the output directory
    call system('mkdir '//trim(dirname))
    
    return
    
  end subroutine mkoutdir
  
!######################################################################
  
  subroutine wrfunc_1element(indx1,indx2,func,ntotal)

    use constants
    use iomod
    use sysinfo
    use gridglobal
    
    implicit none

    integer, intent(in)    :: indx1,indx2
    integer(8), intent(in) :: ntotal
    integer                :: unit,i
    real(dp), intent(in)   :: func(ntotal)
    character(len=80)      :: filename,stem
    
!----------------------------------------------------------------------
! Open the output file
!----------------------------------------------------------------------
    call freeunit(unit)

    select case(ifunc)
    case(1) ! Projector onto an adiabatic state
       write(stem,'(a,i0)') 'adproj',funcsta(1)
       write(filename,'(a,i0,a,i0,a)') &
            trim(stem)//'/'//trim(stem)//'_',indx1,'_',indx2,'.dat'
    case(2) ! Adiabatic state excitation operator
       write(stem,'(a,2(i0))') 'adexci',funcsta(1),funcsta(2)
       write(filename,'(a,i0,a,i0,a)') &
            trim(stem)//'/'//trim(stem)//'_',indx1,'_',indx2,'.dat'
    case(3) ! Rotational energies
       write(stem,'(a,i0)') 'EJ',indx1
       write(filename,'(a,i0,a)') &
            trim(stem)//'/'//trim(stem)//'_K',indx2,'.dat'
    end select

    open(unit,file=filename,form='unformatted',status='unknown')
    
!----------------------------------------------------------------------
! Write the function to file
!----------------------------------------------------------------------
    do i=1,ntotal
       write(unit) func(i)
    enddo
    
!----------------------------------------------------------------------
! Close the output file
!----------------------------------------------------------------------
    close(unit)
    
    return
    
  end subroutine wrfunc_1element
    
!######################################################################
  
end module func

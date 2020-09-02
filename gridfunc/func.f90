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
       
    case(3) ! Rotational energies
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
    do s1=1,nsta
       do s2=s1,nsta
          
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
    real(dp)               :: q(nfuncmode)
    real(dp), allocatable  :: func(:)

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(func(ntotal))
    func=0.0d0

!----------------------------------------------------------------------
! Compute the matrix representation of the function on the
! (sub) direct product grid
!----------------------------------------------------------------------
    ! Loop over grid points
    do i=1,ntotal

       ! Mode values at the current grid point
       call mode_values(i,q)

       ! Value of the rotational energy at the current grid point
       func(i)=rotational_energy(q)
       
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
    adt=adtmatrix(q1)

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
    use gridglobal
    
    implicit none

    integer              :: m,m1,i,j,iexci1,iexci2
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
    adt=adtmatrix(q1)

!----------------------------------------------------------------------
! Diabatic representation of the adiabatic state excitation operator
!----------------------------------------------------------------------
    ! Adiabatic state indices
    iexci1=funcsta(1)
    iexci2=funcsta(2)
    
    ! Form the matrix representation of the excitation operator
    do i=1,nsta
       do j=1,nsta
          proj(j,i)=adt(j,iexci1)*adt(i,iexci2) &
               +adt(i,iexci1)*adt(j,iexci2)
       enddo
    enddo
       
    return
    
  end function adiabatic_excitation
  
!######################################################################

  function rotational_energy(q) result(EJ)

    use constants
    use iomod
    use sysinfo
    use gridglobal
    
    implicit none

    integer              :: m,m1,i,j,error
    real(dp), intent(in) :: q(nfuncmode)
    real(dp)             :: EJ
    real(dp)             :: q1(nmodes)
    real(dp)             :: x(ncoo)
    real(dp)             :: itensor(3,3)
    real(dp)             :: iteig(3)
    real(dp)             :: work(9)
    real(dp)             :: totmass,com(3)
    real(dp)             :: xx,xy,xz,xm
    real(dp)             :: A,B,C
    
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
          com(j)=com(j)+mass(i*3)*x(i*3-3+j)
       enddo
    enddo
    com=com/totmass

    ! Shift the origin to the centre of mass
    do i=1,natm
       do j=1,3
          x(i*3-3+j)=x(i*3-3+j)-com(j)
       enddo
    enddo

    com=0.0d0
    do i=1,natm
       do j=1,3
          com(j)=com(j)+mass(i*3)*x(i*3-3+j)
       enddo
    enddo
    com=com/totmass
    
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
    
  end function rotational_energy
    
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
  
  subroutine wrfunc_1element(s1,s2,func,ntotal)

    use constants
    use iomod
    use sysinfo
    use gridglobal
    
    implicit none

    integer, intent(in)    :: s1,s2
    integer(8), intent(in) :: ntotal
    integer                :: unit
    real(dp), intent(in)   :: func(ntotal)
    character(len=80)      :: filename,stem
    
!----------------------------------------------------------------------
! Open the output file
!----------------------------------------------------------------------
    call freeunit(unit)

    select case(ifunc)
    case(1) ! Projector onto an adiabatic state
       write(stem,'(a,i0)') 'adproj',funcsta(1)
    case(2) ! Adiabatic state excitation operator
       write(stem,'(a,2(i0))') 'adexci',funcsta(1),funcsta(2)
    case(3) ! Rotational energies
       write(stem,'(a,i0)') 'EJ',Jval
    end select
    
    write(filename,'(a,i0,a,i0,a)') &
         trim(stem)//'/'//trim(stem)//'_',s1,'_',s2,'.dat'
    
    open(unit,file=filename,form='unformatted',status='unknown')

!----------------------------------------------------------------------
! Write the function to file
!----------------------------------------------------------------------
    write(unit) func
    
!----------------------------------------------------------------------
! Close the output file
!----------------------------------------------------------------------
    close(unit)
    
    return
    
  end subroutine wrfunc_1element
    
!######################################################################
  
end module func

!######################################################################
! cartgrad: routines to write the gradient and non-adiabatic coupling
!           vectors to xyz files for visualisation
!######################################################################
module cartgrad

  implicit none
  
contains

!######################################################################

  subroutine write_cartgrad

    use constants
    use sysinfo
    use parameters
    use iomod
    
    implicit none

    integer           :: s,s1,s2,i,j,iout
    real(dp)          :: grad(ncoo)
    real(dp)          :: norm
    logical           :: found
    character(len=60) :: outfile
    
!----------------------------------------------------------------------
! Create the cartgrad directory
!----------------------------------------------------------------------
    ! Works with intel
    inquire(directory='cartgrad',exist=found)

    ! Works with gfortran
    !inquire(file='cartgrad/.',exist=found)

    if (found) then
       call system('rm -rf cartgrad/*')
    else
       call system('mkdir cartgrad')
    endif

!----------------------------------------------------------------------
! Gradients
!----------------------------------------------------------------------
    call freeunit(iout)

    ! Loop over states
    do s=1,nsta

       ! Open the output file
       write(outfile,'(a,i0,a)') 'cartgrad/grad_',s,'.xyz'
       open(iout,file=outfile,form='formatted',status='unknown')

       ! Gradient in Cartesian coordinates
       grad=matmul(transpose(coonm), coeff1(:,s,s,1))

       ! Normalise
       norm=sqrt(dot_product(grad,grad))
       grad=grad/norm

       ! Write the gradient to file
       write(iout,'(i0,/,a,x,F9.6)') natm,'norm:',norm
       do i=1,natm
          write(iout,'(a,6(2x,F12.7))') atlbl(i),&
               (xcoo0(j)/ang2bohr,j=i*3-2,i*3),&
               (grad(j),j=i*3-2,i*3)
       enddo
       
       ! Close the output file
       close(iout)
       
    enddo

!----------------------------------------------------------------------
! Non-adiabatic couplings
!----------------------------------------------------------------------
    call freeunit(iout)

    ! Loop over pairs of states
    do s1=1,nsta-1
       do s2=s1+1,nsta

          ! Open the output file
          write(outfile,'(2(a,i0),a)') 'cartgrad/nact_',s1,'_',s2,&
               '.xyz'

          open(iout,file=outfile,form='formatted',status='unknown')
          
          ! Non-adiabatic coupling in Cartesian coordinates
          grad=matmul(transpose(coonm), coeff1(:,s1,s2,1))
          
          ! Normalise
          norm=sqrt(dot_product(grad,grad))
          grad=grad/norm

          ! Write the gradient to file
          write(iout,'(i0,/,a,x,F9.6)') natm,'norm:',norm
          do i=1,natm
             write(iout,'(a,6(2x,F12.7))') atlbl(i),&
                  (xcoo0(j)/ang2bohr,j=i*3-2,i*3),&
                  (grad(j),j=i*3-2,i*3)
          enddo
          
          ! Close the output file
          close(iout)
          
       enddo
    enddo
    
    return
    
  end subroutine write_cartgrad
    
!######################################################################
  
end module cartgrad

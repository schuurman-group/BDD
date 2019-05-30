module utils

  implicit none

contains

!######################################################################
! nlines: returns the number of lines in a formatted file
!######################################################################
  
  function nlines(filename)

    use constants
    use iomod
    use parsemod
    
    implicit none

    integer          :: nlines,unit
    character(len=*) :: filename

    call freeunit(unit)
    open(unit,file=filename,form='formatted',status='old')

    nlines=0
5   read(unit,*,end=10)
    nlines=nlines+1
    goto 5
10  continue

    close(unit)
    
    return
    
  end function nlines
  
!######################################################################

  subroutine dsortindxa1(order,ndim,arrin,indx)

    use constants

    implicit none
    
    character(1), intent(in) :: order
    integer, intent(in) :: ndim
    real(dp), dimension(ndim), intent(in)   :: arrin
    integer, dimension(ndim), intent(inout) :: indx
    
    integer :: i,l,ir,indxt,j
    real(dp) :: q

!!$ The subroutine is taken from the NR p233, employs heapsort.

    if (ndim.eq.1) then
       indx(1)=1
       return
    endif
       
    do i= 1,ndim
       indx(i)=i
    end do
    
    l=ndim/2+1
    ir=ndim
    
    if(order .eq. 'D') then
       
10     continue
       if(l .gt. 1) then
          l=l-1
          indxt=indx(l)
          q=arrin(indxt)
       else
          indxt=indx(ir)
          q=arrin(indxt)
          indx(ir)=indx(1)
          ir=ir-1
          if(ir .eq. 1) then
             indx(1)=indxt
             return
          end if
       end if
       
       i=l
       j=l+l
       
20     if(j .le. ir) then
          if(j .lt. ir) then
             if(arrin(indx(j)) .gt. arrin(indx(j+1))) j=j+1 !
          end if
          if(q .gt. arrin(indx(j))) then !
             indx(i)=indx(j)
             i=j
             j=j+j
          else
             j=ir+1
          end if
          go to 20
       end if
       indx(i)=indxt
       go to 10
       
    elseif(order .eq. 'A') then
       
100    continue
       if(l .gt. 1) then
          l=l-1
          indxt=indx(l)
          q=arrin(indxt)
       else
          indxt=indx(ir)
          q=arrin(indxt)
          indx(ir)=indx(1)
          ir=ir-1
          if(ir .eq. 1) then
             indx(1)=indxt
             return
          end if
       end if
       
       i=l
       j=l+l
       
200    if(j .le. ir) then
          if(j .lt. ir) then
             if(arrin(indx(j)) .lt. arrin(indx(j+1))) j=j+1 !
          end if
          if(q .lt. arrin(indx(j))) then !
             indx(i)=indx(j)
             i=j
             j=j+j
          else
             j=ir+1
          end if
          go to 200
       end if
       indx(i)=indxt
       go to 100
       
    end if
       
  end subroutine dsortindxa1

!######################################################################

  subroutine i8sortindxa1(order,ndim,arrin,indx)

    use constants

    implicit none
    
    character(1), intent(in) :: order
    integer*8, intent(in) :: ndim
    integer*8, dimension(ndim), intent(in)    :: arrin
    integer*8, dimension(ndim), intent(inout) :: indx
    
    integer*8 :: i,l,ir,indxt,j
    integer*8 :: q

!!$ The subroutine is adapted from the NR p233, employs heapsort.
    
    if (ndim.eq.1) return

    do i= 1,ndim
       indx(i)=i
    end do
    
    l=ndim/2+1
    ir=ndim
    
    if(order .eq. 'D') then
       
10     continue
       if(l .gt. 1) then
          l=l-1
          indxt=indx(l)
          q=arrin(indxt)
       else
          indxt=indx(ir)
          q=arrin(indxt)
          indx(ir)=indx(1)
          ir=ir-1
          if(ir .eq. 1) then
             indx(1)=indxt
             return
          end if
       end if
       
       i=l
       j=l+l
       
20     if(j .le. ir) then
          if(j .lt. ir) then
             if(arrin(indx(j)) .gt. arrin(indx(j+1))) j=j+1 !
          end if
          if(q .gt. arrin(indx(j))) then !
             indx(i)=indx(j)
             i=j
             j=j+j
          else
             j=ir+1
          end if
          go to 20
       end if
       indx(i)=indxt
       go to 10
       
    elseif(order .eq. 'A') then
       
100    continue
       if(l .gt. 1) then
          l=l-1
          indxt=indx(l)
          q=arrin(indxt)
       else
          indxt=indx(ir)
          q=arrin(indxt)
          indx(ir)=indx(1)
          ir=ir-1
          if(ir .eq. 1) then
             indx(1)=indxt
             return
          end if
       end if
       
       i=l
       j=l+l
       
200    if(j .le. ir) then
          if(j .lt. ir) then
             if(arrin(indx(j)) .lt. arrin(indx(j+1))) j=j+1 !
          end if
          if(q .lt. arrin(indx(j))) then !
             indx(i)=indx(j)
             i=j
             j=j+j
          else
             j=ir+1
          end if
          go to 200
       end if
       indx(i)=indxt
       go to 100
       
    end if
       
  end subroutine i8sortindxa1

!######################################################################

  subroutine invert_matrix(mat,invmat,dim,pseudo)

    use constants
    use channels
    use iomod
    
    implicit none

    integer                      :: dim,info,i
    real(dp), dimension(dim,dim) :: mat,invmat,tmp,umat,vtmat,&
                                    sigmaplus
    real(dp), dimension(dim)     :: sigma
    real(dp), dimension(5*dim)   :: work
    real(dp), parameter          :: thrsh=1e-10_dp
    logical                      :: pseudo
    
!----------------------------------------------------------------------    
! SVD of the input matrix
!----------------------------------------------------------------------
    tmp=mat
    call dgesvd('A','A',dim,dim,tmp,dim,sigma,umat,dim,vtmat,dim,&
         work,5*dim,info)

    ! Exit if the SVD failed
    if (info.ne.0) then
       errmsg='SVD failed in subroutine invert_matrix'
       call error_control
    endif

!----------------------------------------------------------------------
! Moore-Penrose inverse
!----------------------------------------------------------------------
    sigmaplus=0.0d0
    pseudo=.false.
    do i=1,dim
       if (abs(sigma(i)).lt.thrsh) then
          sigmaplus(i,i)=0.0d0
          pseudo=.true.
       else
          sigmaplus(i,i)=1.0d0/sigma(i)
       endif
    enddo

    invmat=matmul(transpose(vtmat),matmul(sigmaplus,transpose(umat)))
    
    return
    
  end subroutine invert_matrix

!######################################################################

  subroutine sqrt_matrix(mat,sqrtmat,dim)

    use constants
    use channels
    use iomod
    
    implicit none

    integer                      :: dim,info,i
    real(dp), dimension(dim,dim) :: mat,sqrtmat,eigvec,dmat
    real(dp), dimension(dim)     :: lambda
    real(dp), dimension(3*dim)   :: work

!----------------------------------------------------------------------
! Diagonalisation of the input matrix
!----------------------------------------------------------------------
    eigvec=mat
    call dsyev('V','U',dim,eigvec,dim,lambda,work,3*dim,info)

    ! Exit if the diagonalisation failed
    if (info.ne.0) then
       errmsg='Diagonalisation failed in subroutine sqrt_matrix'
       call error_control
    endif

!----------------------------------------------------------------------
! Square root of the input matrix
!----------------------------------------------------------------------
    dmat=0.0d0
    do i=1,dim
       dmat(i,i)=sqrt(lambda(i))
    enddo

    sqrtmat=matmul(eigvec,matmul(dmat,transpose(eigvec)))

    return
    
  end subroutine sqrt_matrix

!######################################################################

  subroutine diag_matrix(mat,eigval,dim)

    use constants
    use iomod
    
    implicit none

    integer, intent(in)          :: dim
    integer                      :: e2,error
    real(dp), dimension(dim,dim) :: mat
    real(dp), dimension(dim)     :: eigval
    real(dp), dimension(dim,dim) :: tmp
    real(dp), dimension(3*dim)   :: work

    tmp=mat
    
    call dsyev('V','U',dim,mat,dim,eigval,work,3*dim,error)

    if (error.ne.0) then
       errmsg='Diagonalisation failed in subroutine diag_matrix'
       call error_control
    endif
    
    return
    
  end subroutine diag_matrix
    
!######################################################################
      
end module utils

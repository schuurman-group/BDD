module gridglobal

  use constants

  implicit none
  
  ! DVR information
  real(dp), allocatable :: dvrpar(:,:)
  integer, allocatable  :: idvr(:)
  integer, allocatable  :: ndvr(:)
  integer, parameter    :: maxdvrpar=10
  
  ! Function type
  integer               :: ifunc

  ! Electronic state indices
  integer               :: funcsta(2)

  ! Mode indices
  integer               :: nfuncmode
  integer, allocatable  :: funcmode(:)
  
end module gridglobal

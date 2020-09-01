module gridglobal

  use constants

  implicit none
  
  ! DVR basis information
  real(dp), allocatable :: dvrpar(:,:)
  integer, allocatable  :: idvr(:)
  integer, allocatable  :: ndvr(:)
  integer, parameter    :: maxdvrpar=10

  ! DVR grids
  integer               :: maxgrid
  real(dp), allocatable :: grid(:,:)
  
  ! Function type
  integer               :: ifunc

  ! Electronic state indices
  integer               :: funcsta(2)

  ! Mode indices
  integer               :: nfuncmode
  integer, allocatable  :: funcmode(:)

  ! Rotational energies
  integer               :: Jval
  
end module gridglobal

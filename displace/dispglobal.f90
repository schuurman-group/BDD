module dispglobal

  use constants
  
  implicit none
  
  ! Cut information
  integer              :: icut
  integer              :: npnts
  integer, allocatable :: cut_mask(:,:)
  real(dp)             :: dq
  
end module dispglobal

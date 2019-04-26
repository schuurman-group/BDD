module pltglobal

  use constants
  
  implicit none
  
  ! Plotting parameters
  integer               :: mplt,npnts,si,sf,surftyp
  real(dp)              :: qi,qf,ei,ef
  real(dp), allocatable :: surf(:,:)
  logical               :: leps
  
end module pltglobal

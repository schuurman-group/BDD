module pltglobal

  use constants
  
  implicit none
  
  ! Plotting parameters
  integer               :: mplt,npnts,si,sf,surftyp
  real(dp)              :: qi,qf,ei,ef
  real(dp), allocatable :: surf(:,:)
  logical               :: leps

  ! Ab initio diabatic potential values
  integer               :: ndat,nabinit
  real(dp), allocatable :: wq0(:)
  real(dp), allocatable :: wdisp(:,:,:)
  real(dp), allocatable :: qvec(:,:)
  real(dp), allocatable :: abinit(:,:,:)
  integer, allocatable  :: iabinit(:)
  
end module pltglobal

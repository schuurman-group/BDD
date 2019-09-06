module pltglobal

  use constants
  
  implicit none
  
  ! Plotting parameters
  integer               :: mplt,mplt2,npnts,si,sf,surftyp
  real(dp)              :: qi,qf,ei,ef
  real(dp), allocatable :: surf(:,:)
  logical               :: ldiagcut
  logical               :: leps

  ! Ab initio diabatic potential values
  integer               :: ndat,nabinit
  real(dp), allocatable :: wq0(:)
  real(dp), allocatable :: wdisp(:,:,:)
  real(dp), allocatable :: qvec(:,:)
  real(dp), allocatable :: abinit(:,:,:)
  integer, allocatable  :: iabinit(:)

  ! Diabatic dipole surfaces
  integer               :: dipsta1,dipsta2
  real(dp), allocatable :: ddisp(:,:,:,:)
  logical               :: ldip

  ! Diabatic coupling surfaces
  integer               :: dcpsta1,dcpsta2
  
end module pltglobal

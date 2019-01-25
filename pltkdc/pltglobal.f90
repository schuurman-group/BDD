module pltglobal

  use constants
  
  implicit none

  ! Vertical excitation energies
  real(dp), allocatable :: e0(:)
  
  ! Coupling coeffients
  real(dp), allocatable :: kappa(:,:)
  real(dp), allocatable :: lambda(:,:,:)
  real(dp), allocatable :: gamma(:,:,:)
  real(dp), allocatable :: mu(:,:,:,:)

  ! Plotting parameters
  integer               :: mplt,npnts,si,sf,surftyp
  real(dp)              :: qi,qf,ei,ef
  real(dp), allocatable :: surf(:,:)
  logical               :: leps
  
end module pltglobal

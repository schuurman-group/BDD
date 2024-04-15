module parameters

  use constants

  implicit none

  ! Vertical excitation energies
  real(dp), allocatable :: e0(:)
  
  ! Coupling coeffients
  real(dp), allocatable :: kappa(:,:)
  real(dp), allocatable :: lambda(:,:,:)
  real(dp), allocatable :: gamma(:,:,:)
  real(dp), allocatable :: mu(:,:,:,:)
  real(dp), allocatable :: iota(:,:)
  real(dp), allocatable :: tau(:,:,:)
  real(dp), allocatable :: epsilon(:,:)
  real(dp), allocatable :: xi(:,:,:)

  ! Dipole matrix expansion coefficients
  real(dp), allocatable :: dip0(:,:,:)
  real(dp), allocatable :: dip1(:,:,:,:)
  real(dp), allocatable :: dip2(:,:,:,:,:)
  real(dp), allocatable :: dip3(:,:,:,:)
  real(dp), allocatable :: dip4(:,:,:,:)

  ! New coupling coefficients
  integer               :: order1
  real(dp), allocatable :: coeff1(:,:,:,:)
  real(dp), allocatable :: coeff2(:,:,:,:)
  
end module parameters

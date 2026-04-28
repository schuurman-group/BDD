module parameters

  use constants

  implicit none

  ! Vertical excitation energies
  real(dp), allocatable :: e0(:)

  ! Off-diagonal zeroth-order constants (populated by $reexpand;
  ! zero by default). Symmetric, with zero diagonal.
  real(dp), allocatable :: e0_off(:,:)

  ! Coupling coeffients
  integer               :: order1
  real(dp), allocatable :: coeff1(:,:,:,:)
  real(dp), allocatable :: coeff2(:,:,:,:)

  ! Dipole matrix expansion coefficients
  real(dp), allocatable :: dip0(:,:,:)
  real(dp), allocatable :: dip1(:,:,:,:)
  real(dp), allocatable :: dip2(:,:,:,:,:)
  real(dp), allocatable :: dip3(:,:,:,:)
  real(dp), allocatable :: dip4(:,:,:,:)

end module parameters

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
  
end module parameters

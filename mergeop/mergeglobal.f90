module mergeglobal

  use constants

  implicit none

  ! KDC file names
  character(len=120) :: abinA,abinB
  
  ! Subsystem dimensions
  integer :: nmodesA,nmodesB
  integer :: nstaA,nstaB
  integer :: ncooA,ncooB
  integer :: natmA,natmB

  ! Frequencies
  real(dp), allocatable :: freqA(:),freqB(:)
  
  ! Vertical excitation energies
  real(dp), allocatable :: e0A(:),e0B(:)
  
  ! Coupling coeffients
  integer               :: order1A,order1B
  real(dp), allocatable :: coeff1A(:,:,:,:),coeff1B(:,:,:,:)
  real(dp), allocatable :: coeff2A(:,:,:,:),coeff2B(:,:,:,:)

  ! Dipole fitting
  logical               :: ldipA,ldipB
  real(dp), allocatable :: dip0A(:,:,:),dip0B(:,:,:)
  real(dp), allocatable :: dip1A(:,:,:,:),dip1B(:,:,:,:)
  real(dp), allocatable :: dip2A(:,:,:,:,:),dip2B(:,:,:,:,:)
  real(dp), allocatable :: dip3A(:,:,:,:),dip3B(:,:,:,:)
  real(dp), allocatable :: dip4A(:,:,:,:),dip4B(:,:,:,:)
  
  ! Symmetry information
  integer, allocatable :: coeff1_maskA(:,:,:,:),coeff1_maskB(:,:,:,:)
  integer, allocatable :: coeff2_maskA(:,:,:,:),coeff2_maskB(:,:,:,:)
  integer, allocatable :: dip0_maskA(:,:,:),dip0_maskB(:,:,:)
  integer, allocatable :: dip1_maskA(:,:,:,:),dip1_maskB(:,:,:,:)
  integer, allocatable :: dip2_maskA(:,:,:,:,:),dip2_maskB(:,:,:,:,:)
  integer, allocatable :: dip3_maskA(:,:,:,:),dip3_maskB(:,:,:,:)
  integer, allocatable :: dip4_maskA(:,:,:,:),dip4_maskB(:,:,:,:)
    
end module mergeglobal

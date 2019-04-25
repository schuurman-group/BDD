module kdcglobal

  use constants
  
  implicit none

  ! Number of blockdiag files
  integer                         :: nfiles

  ! Filenames
  character(len=120)              :: setfile
  character(len=120), allocatable :: bdfiles(:)
  logical                         :: lsetfile
  
  ! Model type and mode <-> file mappings
  integer, allocatable            :: coefftyp(:)
  integer, allocatable            :: findx1m(:,:)
  integer, allocatable            :: findx2m(:,:,:)
  integer, allocatable            :: nfiles1m(:)
  integer, allocatable            :: nfiles2m(:,:)
  integer                         :: maxfiles1m,maxfiles2m
  logical                         :: lbilinear
  
  ! Adiabatic potential at Q0
  real(dp), allocatable           :: q0pot(:)
  
  ! Diabatic potential matrices at the displaced geometries
  real(dp), allocatable           :: diabpot(:,:,:)

  ! Normal mode coordinates
  real(dp), allocatable           :: qvec(:,:)
  
  ! Coupling coeffients
  real(dp), allocatable           :: kappa(:,:)
  real(dp), allocatable           :: lambda(:,:,:)
  real(dp), allocatable           :: gamma(:,:,:)
  real(dp), allocatable           :: mu(:,:,:,:)
  real(dp), allocatable           :: iota(:,:)
  real(dp), allocatable           :: tau(:,:,:)
  real(dp), allocatable           :: epsilon(:,:)
  real(dp), allocatable           :: xi(:,:,:)
  
  ! Parameterisation algorithm
  integer                         :: ialgor
  
end module kdcglobal

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
  logical                         :: lbilinear
  
  ! Adiabatic potential at Q0
  real(dp), allocatable           :: q0pot(:)
  
  ! Diabatic potential matrices at the displaces geometries
  real(dp), allocatable           :: diabpot(:,:,:)

  ! Normal mode coordinates
  real(dp), allocatable           :: qvec(:,:)
  
  ! Coupling coeffients
  real(dp), allocatable           :: kappa(:,:)
  real(dp), allocatable           :: lambda(:,:,:)
  real(dp), allocatable           :: gamma(:,:,:)
  real(dp), allocatable           :: mu(:,:,:,:)
  
end module kdcglobal

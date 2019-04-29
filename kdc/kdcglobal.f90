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

  ! Diabatic dipoles
  real(dp), allocatable           :: diabdip(:,:,:,:)
  logical                         :: ldipfit

  ! Adiabatic dipoles at Q0
  real(dp), allocatable           :: q0dip(:,:)
  
  ! Parameterisation algorithm
  integer                         :: ialgor

  ! RMSDs of the normal equations fit
  real(dp)                        :: rmsd
  real(dp), allocatable           :: rmsd1m(:)
  real(dp), allocatable           :: rmsd2m(:,:)
  
end module kdcglobal

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
  
  ! Zeroth-order parameter shifts
  real(dp), allocatable           :: shift0(:)
  logical                         :: lshift
  
  ! Diabatic potential matrices at the displaced geometries
  real(dp), allocatable           :: diabpot(:,:,:)

  ! Normal mode coordinates
  real(dp), allocatable           :: qvec(:,:)

  ! Diabatic dipoles
  real(dp), allocatable           :: diabdip(:,:,:,:)
  logical                         :: ldipfit

  ! Parameterisation algorithm
  integer                         :: ialgor
  integer                         :: iweight
  real(dp)                        :: wfac
  
  ! RMSDs of the normal equations fit
  real(dp)                        :: rmsd
  real(dp), allocatable           :: rmsd1m(:)
  real(dp), allocatable           :: rmsd2m(:,:)
  
end module kdcglobal

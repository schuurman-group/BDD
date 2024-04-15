module kdcglobal

  use constants
  
  implicit none

  ! Program used to generate the BDD potentials
  ! iprog = 1 <-> BDD
  !         2 <-> GRaCI
  integer                         :: iprog
  
  ! Number of BDD output files
  integer                         :: nfiles

  ! Number of displaced geometries
  integer                         :: ngeom
  
  ! Filenames
  character(len=120)              :: setfile
  character(len=120), allocatable :: bdfiles(:)
  logical                         :: lsetfile
  
  ! Model type and mode <-> file mappings
  integer, allocatable            :: coefftyp(:)
  integer, allocatable            :: findx1m(:,:)
  integer, allocatable            :: findx2m(:,:,:)
  integer, allocatable            :: ngeom1m(:)
  integer, allocatable            :: ngeom2m(:,:)
  integer                         :: maxfiles1m,maxfiles2m
  logical                         :: lbilinear
  
  ! Adiabatic potential at Q0
  real(dp), allocatable           :: q0pot(:)
  
  ! Zeroth-order parameter shifts
  real(dp), allocatable           :: shift0(:,:)
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

  ! Label stem to use when parsing GRaCI output
  character(len=60)               :: lblstem

  ! Block diagonalisation variables
  integer                         :: nbd(2)
  integer                         :: maxbd
  integer, allocatable            :: ibd(:,:)
  real(dp), allocatable           :: Tmat(:,:,:)
  logical                         :: lblockdiag
  
end module kdcglobal

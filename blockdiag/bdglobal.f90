module bdglobal

  use constants
  use import_gamess
  
  implicit none

  ! Dimensions
  integer                         :: nsta,nsta_disp,nsta_ref,&
                                     maxdet,nmo_ref,nmo_disp,&
                                     nel
  integer, allocatable            :: ndet_ref(:),ndet_disp(:)
  
  ! Filenames
  character(len=120), allocatable :: adetref(:),adetdisp(:)
  character(len=120)              :: amosref,amosdisp
  logical                         :: lbinary
  
  ! Multigrid internal gamess derived data types
  type(gam_structure)             :: gam_ref,gam_disp

  ! Wavefunction arrays
  integer, allocatable            :: det_ref(:,:,:),det_disp(:,:,:)
  real(dp), allocatable           :: c_ref(:,:),c_disp(:,:)

  ! Alpha and beta string arrays
  integer                         :: nalpha,nbeta
  integer, allocatable            :: iocca_ref(:,:,:),iocca_disp(:,:,:),&
                                     ioccb_ref(:,:,:),ioccb_disp(:,:,:)
    
  ! Wavefunction norms
  real(dp), allocatable           :: norm_ref(:),norm_disp(:)

  ! Reference geometry transformation matrix
  real(dp), allocatable           :: reftrans(:,:)
  character(len=120)              :: areftrans
  logical                         :: lreftrans,lrdreftrans
  
  ! MO overlaps
  real(dp), allocatable           :: smo(:,:)

  ! Electronic state overlaps
  integer                         :: ioverlap
  real(dp), allocatable           :: spsi(:,:)
  
  ! Adiabatic potential matrix
  real(dp), allocatable           :: Vmat(:,:),Vmat1(:,:)

  ! Diabatic potential matrix
  real(dp), allocatable           :: Wmat(:,:)
  logical                         :: ldiabpot

  ! Dipole matrix
  real(dp), allocatable           :: ddip(:,:,:)
  real(dp), allocatable           :: adip(:,:,:),adip1(:,:,:)
  logical                         :: ldipole
  
  ! ADT matrix
  real(dp), allocatable           :: adt(:,:)

  ! Hadamard screening threshold
  real(dp)                        :: dthresh

  ! Wavefunction norm cutoff
  real(dp)                       :: normcut
  logical                        :: ltruncate
   
end module bdglobal

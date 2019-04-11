module bdglobal

  use constants
  use import_gamess
  
  implicit none

  ! Dimensions
  integer                         :: nsta,maxdet,nmo_ref,nmo_disp,&
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
  integer*8, allocatable          :: ilbla_ref(:,:),ilbla_disp(:,:),&
                                     ilblb_ref(:,:),ilblb_disp(:,:)
  integer                         :: na_ref,na_disp,nb_ref,nb_disp
  integer, allocatable            :: ia_ref(:,:),ib_ref(:,:),&
                                     ia_disp(:,:),ib_disp(:,:)
  integer, allocatable            :: stringa_ref(:,:),stringb_ref(:,:),&
                                     stringa_disp(:,:),stringb_disp(:,:)
  real(dp), allocatable           :: afac(:,:),bfac(:,:)
  
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
  real(dp), allocatable           :: Vmat(:,:)

  ! Diabatic potential matrix
  real(dp), allocatable           :: Wmat(:,:)
  logical                         :: ldiabpot
  
  ! ADT matrix
  real(dp), allocatable           :: adt(:,:)
  
end module bdglobal

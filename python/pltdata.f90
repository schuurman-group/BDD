!######################################################################
! pltdata: module to store ab initio data read from the binary
!          parameter file (for use by the f2py wrapper routines)
!######################################################################

module pltdata

  use constants

  implicit none

  ! Ab initio data
  integer                         :: ndat
  real(dp), allocatable           :: wq0(:)
  real(dp), allocatable           :: wdisp(:,:,:)
  real(dp), allocatable           :: qvec(:,:)

  ! Diabatic dipole ab initio data
  real(dp), allocatable           :: ddisp(:,:,:,:)
  logical                         :: ldip

end module pltdata

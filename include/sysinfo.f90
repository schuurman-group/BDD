module sysinfo

  use constants
  
  implicit none
  
  ! System information
  integer                       :: ncoo,nmodes,natm,nsta
  integer, allocatable          :: atnum(:)
  real(dp), allocatable         :: xcoo0(:),mass(:)
  character(len=2), allocatable :: atlbl(:)

  ! Normal modes
  real(dp), allocatable         :: nmcoo(:,:),coonm(:,:)
  real(dp), allocatable         :: freq(:)
  character(len=3), allocatable :: nmlab(:)
    
  ! Filenames
  character(len=120)            :: freqfile
  
  ! QC calculation info
  integer                       :: freqtyp

end module sysinfo

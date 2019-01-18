module constants
  
  save
  
  integer, parameter     :: dp=selected_real_kind(8)
  integer, parameter     :: lng=selected_int_kind(16)
  
  real(dp), parameter    :: rzero=0._dp
  real(dp), parameter    :: rone=1._dp
  real(dp), parameter    :: pi=3.14159265358979_dp
  real(dp), parameter    :: ang2bohr=1.88972612d0
  real(dp), parameter    :: invcm2ev=1.23985e-4_dp
  real(dp), parameter    :: eh2ev=27.2113845d0
  real(dp), parameter    :: c_au=137.0359991d0
  complex(dp), parameter :: ci=(0._dp,1._dp)
  complex(dp), parameter :: czero=(0._dp,0._dp)
  complex(dp), parameter :: cone=(1._dp,0._dp)
  
end module constants

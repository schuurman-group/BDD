#-----------------------------------------------------------------------
# Compiler-specific Fortran flags
#-----------------------------------------------------------------------

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")

  set(BDD_Fortran_FLAGS
    -cpp
    -ffixed-line-length-none
    -ffree-line-length-none
    -fbacktrace
  )
  set(BDD_PP_DEFS GFORTRAN)

elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")

  set(BDD_Fortran_FLAGS
    -cpp
    -traceback
    "-diag-disable=8290,8291,10448"
  )
  set(BDD_Fortran_FLAGS_F90 -free)
  set(BDD_PP_DEFS "")

else()
  message(WARNING "Untested Fortran compiler: ${CMAKE_Fortran_COMPILER_ID}")
  set(BDD_Fortran_FLAGS -cpp)
  set(BDD_PP_DEFS "")
endif()

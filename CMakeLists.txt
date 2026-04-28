cmake_minimum_required(VERSION 3.14)
project(BDD LANGUAGES Fortran)

#-----------------------------------------------------------------------
# Compiler flags
#-----------------------------------------------------------------------
include(cmake/CompilerFlags.cmake)

# Place .mod files in a central directory
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
include_directories(${CMAKE_BINARY_DIR}/modules)

# Default to Release build
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

#-----------------------------------------------------------------------
# External libraries
#-----------------------------------------------------------------------
find_package(LAPACK REQUIRED)
find_package(OpenMP COMPONENTS Fortran)

#-----------------------------------------------------------------------
# Helper function: apply project-wide compile options to a target
#-----------------------------------------------------------------------
function(bdd_target_setup target)
  target_compile_options(${target} PRIVATE ${BDD_Fortran_FLAGS})
  if(BDD_PP_DEFS)
    target_compile_definitions(${target} PRIVATE ${BDD_PP_DEFS})
  endif()
  if(OpenMP_Fortran_FOUND)
    target_link_libraries(${target} PRIVATE OpenMP::OpenMP_Fortran)
  endif()
  # Apply free-form flag (e.g. Intel -free) only to .f90 sources
  if(BDD_Fortran_FLAGS_F90)
    get_target_property(_srcs ${target} SOURCES)
    foreach(_src ${_srcs})
      if(_src MATCHES "\\.f90$")
        set_source_files_properties(${_src} PROPERTIES
          COMPILE_OPTIONS "${BDD_Fortran_FLAGS_F90}")
      endif()
    endforeach()
  endif()
endfunction()

#-----------------------------------------------------------------------
# Library: bdd_include (constants, channels, sysinfo)
#-----------------------------------------------------------------------
add_library(bdd_include STATIC
  include/constants.f90
  include/channels.f90
  include/sysinfo.f90
)
bdd_target_setup(bdd_include)

#-----------------------------------------------------------------------
# Library: bdd_iomodules (iomod, parsemod)
#-----------------------------------------------------------------------
add_library(bdd_iomodules STATIC
  iomodules/iomod.f90
  iomodules/parsemod.f90
)
bdd_target_setup(bdd_iomodules)
target_link_libraries(bdd_iomodules PRIVATE bdd_include)

#-----------------------------------------------------------------------
# Library: bdd_utilities (timingmod, utils)
#-----------------------------------------------------------------------
add_library(bdd_utilities STATIC
  utilities/timingmod.f90
  utilities/utils.f90
)
bdd_target_setup(bdd_utilities)
target_link_libraries(bdd_utilities PRIVATE bdd_include bdd_iomodules)

#-----------------------------------------------------------------------
# Library: bdd_ioqc
#-----------------------------------------------------------------------
add_library(bdd_ioqc STATIC
  ioqc/ioqc.f90
)
bdd_target_setup(bdd_ioqc)
target_link_libraries(bdd_ioqc PRIVATE bdd_include bdd_iomodules bdd_utilities)

#-----------------------------------------------------------------------
# Library: bdd_symmetry
#-----------------------------------------------------------------------
add_library(bdd_symmetry STATIC
  symmetry/symmetry.f90
)
bdd_target_setup(bdd_symmetry)
target_link_libraries(bdd_symmetry PRIVATE bdd_include bdd_iomodules)

#-----------------------------------------------------------------------
# Library: bdd_potfuncs (parameters, potfuncs)
#-----------------------------------------------------------------------
add_library(bdd_potfuncs STATIC
  potfuncs/parameters.f90
  potfuncs/potfuncs.f90
)
bdd_target_setup(bdd_potfuncs)
target_link_libraries(bdd_potfuncs PRIVATE
  bdd_include bdd_iomodules bdd_symmetry)

#-----------------------------------------------------------------------
# Library: bdd_opt
#-----------------------------------------------------------------------
add_library(bdd_opt STATIC
  opt/opt.f90
  opt/geomopt.f90
)
bdd_target_setup(bdd_opt)
target_link_libraries(bdd_opt PRIVATE bdd_include bdd_potfuncs)

#-----------------------------------------------------------------------
# Library: bdd_multi (accuracy, printing, timer, lapack, math,
#                     os_integral_operators, gamess_internal,
#                     import_gamess, dgefa, dgedi)
#
# Note: import_gamess.f90 uses Fortran include statements to pull in
# _common.f90 files from the multi/ directory, so we add that
# directory as an include path.
#-----------------------------------------------------------------------
add_library(bdd_multi STATIC
  multi/accuracy.f90
  multi/printing.f90
  multi/timer.f90
  multi/lapack.f90
  multi/math.f90
  multi/os_integral_operators.f90
  multi/gamess_internal.f90
  multi/import_gamess.f90
  multi/dgefa.f
  multi/dgedi.f
)
bdd_target_setup(bdd_multi)
target_include_directories(bdd_multi PRIVATE ${CMAKE_SOURCE_DIR}/multi)

#-----------------------------------------------------------------------
# Library: bdd_kdc_common (kdcglobal, opermod -- shared by kdc and
#                          mergeop targets)
#-----------------------------------------------------------------------
add_library(bdd_kdc_common STATIC
  kdc/kdcglobal.f90
  kdc/opermod.f90
)
bdd_target_setup(bdd_kdc_common)
target_link_libraries(bdd_kdc_common PRIVATE
  bdd_include bdd_iomodules bdd_symmetry bdd_potfuncs)

#=======================================================================
# Executables
#=======================================================================

#-----------------------------------------------------------------------
# blockdiag.x
#-----------------------------------------------------------------------
add_executable(blockdiag.x
  blockdiag/bdglobal.f90
  blockdiag/adt.f90
  blockdiag/mooverlaps.f90
  blockdiag/wfoverlaps.f90
  blockdiag/tamura.f90
  blockdiag/pacher.f90
  blockdiag/detparsing.f90
  blockdiag/blockdiag.f90
)
bdd_target_setup(blockdiag.x)
target_link_libraries(blockdiag.x PRIVATE
  bdd_multi bdd_include bdd_iomodules bdd_utilities bdd_ioqc
  LAPACK::LAPACK)

#-----------------------------------------------------------------------
# displace.x
#-----------------------------------------------------------------------
add_executable(displace.x
  displace/dispglobal.f90
  displace/displace.f90
)
bdd_target_setup(displace.x)
target_link_libraries(displace.x PRIVATE
  bdd_include bdd_iomodules bdd_utilities bdd_ioqc bdd_symmetry
  LAPACK::LAPACK)

#-----------------------------------------------------------------------
# kdc.x
#-----------------------------------------------------------------------
add_executable(kdc.x
  kdc/cartgrad.f90
  kdc/parinfo.f90
  kdc/nmeqmod.f90
  kdc/transform.f90
  kdc/kdc.f90
)
bdd_target_setup(kdc.x)
target_link_libraries(kdc.x PRIVATE
  bdd_kdc_common
  bdd_include bdd_iomodules bdd_utilities bdd_ioqc
  bdd_symmetry bdd_potfuncs bdd_opt
  LAPACK::LAPACK)

#-----------------------------------------------------------------------
# pltkdc.x
#-----------------------------------------------------------------------
add_executable(pltkdc.x
  pltkdc/pltglobal.f90
  pltkdc/pltkdc.f90
)
bdd_target_setup(pltkdc.x)
target_link_libraries(pltkdc.x PRIVATE
  bdd_include bdd_iomodules bdd_symmetry bdd_potfuncs
  LAPACK::LAPACK)

#-----------------------------------------------------------------------
# mergeop.x
#-----------------------------------------------------------------------
add_executable(mergeop.x
  mergeop/mergeglobal.f90
  mergeop/mergeop.f90
)
bdd_target_setup(mergeop.x)
target_link_libraries(mergeop.x PRIVATE
  bdd_kdc_common
  bdd_include bdd_iomodules bdd_symmetry bdd_potfuncs
  LAPACK::LAPACK)

#-----------------------------------------------------------------------
# Install executables into bin/
#-----------------------------------------------------------------------
install(TARGETS blockdiag.x displace.x kdc.x pltkdc.x mergeop.x
  RUNTIME DESTINATION bin)

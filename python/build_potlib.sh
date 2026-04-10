#!/bin/bash
#
# Build the bdd_potlib f2py module
#
# Usage: ./build_potlib.sh [FC]
#
# FC: Fortran compiler (default: gfortran)
#

set -e

FC=${1:-gfortran}

SRCDIR=$(cd "$(dirname "$0")/.." && pwd)
BUILDDIR=$(mktemp -d)
WRAPDIR=$(mktemp -d)

trap "rm -rf $BUILDDIR $WRAPDIR" EXIT

# Compiler flags
if [[ "$FC" == *"gfortran"* ]]; then
    FFLAGS="-cpp -ffree-line-length-none -fPIC -DGFORTRAN -O3"
elif [[ "$FC" == *"ifx"* ]] || [[ "$FC" == *"ifort"* ]]; then
    FFLAGS="-cpp -free -fPIC -O3"
else
    FFLAGS="-cpp -fPIC -O3"
fi

# Python/numpy paths
PYTHON_INC=$(python3 -c "import sysconfig; print(sysconfig.get_path('include'))")
NUMPY_INC=$(python3 -c "import numpy; print(numpy.get_include())")
F2PY_SRC=$(python3 -c "import numpy.f2py; import os; print(os.path.join(os.path.dirname(numpy.f2py.__file__), 'src'))")
PYTHON_LIB=$(python3 -c "import sysconfig; print(sysconfig.get_config_var('LDLIBRARY'))")
PYTHON_LIBDIR=$(python3 -c "import sysconfig; print(sysconfig.get_config_var('LIBDIR'))")
EXT_SUFFIX=$(python3 -c "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX'))")

echo "Building bdd_potlib with $FC"
echo "Python include: $PYTHON_INC"
echo "NumPy include:  $NUMPY_INC"
echo "Extension:      $EXT_SUFFIX"

#-----------------------------------------------------------------------
# Step 1: Pre-compile all Fortran dependency modules to object files
#-----------------------------------------------------------------------
echo ""
echo "Step 1: Compiling Fortran sources..."

DEPS=(
    "$SRCDIR/include/constants.f90"
    "$SRCDIR/include/channels.f90"
    "$SRCDIR/include/sysinfo.f90"
    "$SRCDIR/iomodules/iomod.f90"
    "$SRCDIR/symmetry/symmetry.f90"
    "$SRCDIR/potfuncs/parameters.f90"
    "$SRCDIR/potfuncs/potfuncs.f90"
    "$SRCDIR/python/pltdata.f90"
    "$SRCDIR/python/bdd_potlib.f90"
)

OBJS=()
for src in "${DEPS[@]}"; do
    obj="$BUILDDIR/$(basename "${src%.f90}.o")"
    $FC $FFLAGS -J"$BUILDDIR" -I"$BUILDDIR" -c "$src" -o "$obj"
    OBJS+=("$obj")
done

#-----------------------------------------------------------------------
# Step 2: Generate the f2py C wrapper and Fortran wrapper files
#         (no compilation, just code generation)
#-----------------------------------------------------------------------
echo "Step 2: Generating f2py wrappers..."

F2PY_INPUT=(
    "$SRCDIR/include/constants.f90"
    "$SRCDIR/include/channels.f90"
    "$SRCDIR/include/sysinfo.f90"
    "$SRCDIR/iomodules/iomod.f90"
    "$SRCDIR/symmetry/symmetry.f90"
    "$SRCDIR/potfuncs/parameters.f90"
    "$SRCDIR/python/pltdata.f90"
    "$SRCDIR/python/bdd_potlib.f90"
)

python3 -m numpy.f2py \
    "${F2PY_INPUT[@]}" \
    -m bdd_potlib \
    --build-dir "$WRAPDIR" \
    only: init get_dimensions get_e0 get_freq \
          calc_adiab_1d calc_adiab_diag \
          calc_diab_1d calc_diab_diag \
          calc_diabcp_1d calc_diabcp_diag \
          calc_dip_1d calc_dip_diag \
          get_abinit_pot :

# Find the generated files
C_WRAPPER=$(find "$WRAPDIR" -name "bdd_potlibmodule.c" | head -1)
F90_WRAPPER=$(find "$WRAPDIR" -name "bdd_potlib-f2pywrappers2.f90" | head -1)

if [ -z "$C_WRAPPER" ]; then
    echo "Error: f2py failed to generate C wrapper"
    exit 1
fi

#-----------------------------------------------------------------------
# Step 3: Compile the f2py support files, C wrapper, and F90 wrapper
#-----------------------------------------------------------------------
echo "Step 3: Compiling wrappers..."

# fortranobject.c (f2py runtime support)
gcc -fPIC -c "$F2PY_SRC/fortranobject.c" \
    -I"$F2PY_SRC" -I"$NUMPY_INC" -I"$PYTHON_INC" \
    -o "$BUILDDIR/fortranobject.o"

# C wrapper
gcc -fPIC -c "$C_WRAPPER" \
    -I"$F2PY_SRC" -I"$NUMPY_INC" -I"$PYTHON_INC" \
    -o "$BUILDDIR/bdd_potlibmodule.o"

# Fortran 90 wrapper (if it exists)
if [ -n "$F90_WRAPPER" ]; then
    $FC $FFLAGS -J"$BUILDDIR" -I"$BUILDDIR" -c "$F90_WRAPPER" \
        -o "$BUILDDIR/bdd_potlib-f2pywrappers2.o"
    OBJS+=("$BUILDDIR/bdd_potlib-f2pywrappers2.o")
fi

#-----------------------------------------------------------------------
# Step 4: Link everything into a shared library
#-----------------------------------------------------------------------
echo "Step 4: Linking..."

OUTPUT="bdd_potlib${EXT_SUFFIX}"

$FC -shared \
    "$BUILDDIR/bdd_potlibmodule.o" \
    "$BUILDDIR/fortranobject.o" \
    "${OBJS[@]}" \
    -llapack -lblas \
    -o "$OUTPUT"

echo ""
echo "Build successful. Module file:"
ls -la "$OUTPUT"
echo ""
echo "Usage: python pltkdc.py -f <paramfile.dat> -m <mode> [options]"

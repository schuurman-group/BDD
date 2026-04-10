# BDD

Tools for the construction of vibronic coupling Hamiltonians
using block diagonalisation diabatisation:

(1) displace: creation of the cuts along normal modes required for the
              determination of the parameters of the vibronic coupling
              Hamiltonian using either finite differences or fitting.

(2) blockdiag: block diagonalisation diabatisation code.

(3) kdc: calculation of the parameters of the vibronic coupling
         Hamiltonian using the output of the blockdiag program.

(4) pltkdc: plotting of the vibronic coupling Hamiltonian potentials
            using the output of the kdc program.


## Prerequisites

- A Fortran compiler (gfortran or Intel ifx/ifort)
- CMake (>= 3.14)
- LAPACK and BLAS libraries
- OpenMP (optional)

## Installation

Build all programs using the install script:

    ./install_BDD

Build a single program:

    ./install_BDD kdc

Executables are installed to `bin/`.

### Choosing the Fortran compiler

The compiler can be selected by setting the `FC` environment variable
before the first build:

    rm -rf build
    FC=ifx ./install_BDD

Or equivalently via the CMake variable:

    rm -rf build
    ./install_BDD "" -DCMAKE_Fortran_COMPILER=ifx

### Manual CMake workflow

    mkdir build && cd build
    cmake ..
    cmake --build .
    cmake --install . --prefix ..

## Python pltkdc

A Python version of `pltkdc` is provided in the `python/` directory.
It uses matplotlib for plotting and calls the Fortran potential
evaluation routines via f2py.

### Prerequisites

- Python 3
- NumPy (with f2py)
- matplotlib
- LAPACK and BLAS libraries
- A Fortran compiler (gfortran or Intel ifx/ifort)

### Building the f2py module

    cd python
    ./build_potlib.sh

To use a specific Fortran compiler:

    ./build_potlib.sh ifx

This produces a shared library (`bdd_potlib*.so`) in the `python/`
directory.

### Usage

    ./pltkdc.py -f <paramfile.dat> -m <mode> [options]

Options:

    -f FILE            Binary parameter file (.dat)
    -m MODE            Normal mode index for the cut
    -m2 MODE2          Second mode for a diagonal cut
    -xrange QI QF      Coordinate range (default: -7 7)
    -yrange EI EF      Energy/function range
    -npnts N           Number of grid points (default: 1000)
    -si I              Lowest state index to plot
    -sf F              Highest state index to plot
    -adiab             Plot adiabatic potentials (default)
    -diab              Plot diabatic potentials
    -dip S1 S2         Plot diabatic dipole for states S1, S2
    -diabcp S1 S2      Plot diabatic coupling for states S1, S2
    -pdf FILE          Output PDF file name (default: auto)

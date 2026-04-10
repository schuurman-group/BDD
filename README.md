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

### Using the Intel Fortran compiler

Remove any existing build directory and pass the compiler to CMake:

    rm -rf build
    ./install_BDD "" -DCMAKE_Fortran_COMPILER=ifx

### Manual CMake workflow

    mkdir build && cd build
    cmake ..
    cmake --build .
    cmake --install . --prefix ..

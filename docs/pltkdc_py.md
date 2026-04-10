---
layout: default
title: pltkdc.py
---

# pltkdc.py

Python replacement for the Fortran pltkdc program. Uses matplotlib
for plotting and f2py-compiled Fortran routines for potential
evaluation. Model potentials are plotted as lines with ab initio
reference data overlaid as black dots. Plots are displayed on screen
and saved to PDF.

## Prerequisites

- Python 3
- NumPy
- Matplotlib
- gfortran (or Intel Fortran compiler)
- The `bdd_potlib` shared library (see [Build](#build) below)

## Build

The f2py interface must be built before use. From the `python/`
directory:

    ./build_potlib.sh

For Intel Fortran:

    ./build_potlib.sh ifort

This produces a `bdd_potlib` shared library that pltkdc.py imports
at runtime. Run pltkdc.py from the same directory as the shared
library, or add its location to `PYTHONPATH`.

## Usage

    ./pltkdc.py -f <paramfile> -m <mode> [options]

### Required arguments

    -f FILE          Binary parameter file (.dat) produced by kdc
    -m MODE          Normal mode index for the cut

### Optional arguments

    -m2 MODE2        Second mode for a diagonal cut
    -xrange QI QF    Coordinate range (default: -7 7)
    -yrange EI EF    Energy/function range (default: auto)
    -npnts N         Number of grid points (default: 1000)
    -si I            Lowest state index to plot (default: 1)
    -sf F            Highest state index to plot (default: all)
    -adiab           Plot adiabatic potentials (default)
    -diab            Plot diabatic potentials
    -dip S1 S2       Plot diabatic dipole for states S1, S2
    -diabcp S1 S2    Plot diabatic coupling for states S1, S2
    -pdf FILE        Output PDF file name (default: auto)

## Plot types

**Adiabatic potentials** (default): eigenvalues of the diabatic
potential matrix as a function of the normal coordinate.

**Diabatic potentials** (`-diab`): diagonal elements of the
diabatic potential matrix.

**Diabatic couplings** (`-diabcp S1 S2`): off-diagonal element
W(S1,S2) of the diabatic potential matrix.

**Diabatic dipoles** (`-dip S1 S2`): diabatic dipole matrix element
between states S1 and S2, with x, y, z components plotted
separately.

## Diagonal cuts

When `-m2` is given, both modes are displaced simultaneously along
the diagonal: Q(m) = Q(m2) = coord / sqrt(2). The x-axis shows the
combined coordinate.

## Output

Plots are displayed on screen interactively before being saved.
PDF files are written to the current directory with auto-generated
names:

    plot_q<m>.pdf                              1D cut
    plot_q<m>_q<m2>.pdf                        diagonal cut
    plot_q<m>_s<S1>_s<S2>_diabcp.pdf           diabatic coupling
    plot_q<m>_s<S1>_s<S2>_dip.pdf              diabatic dipole

A custom filename can be specified with `-pdf`.

## See also

[pltkdc](pltkdc) -- Fortran version using gnuplot.

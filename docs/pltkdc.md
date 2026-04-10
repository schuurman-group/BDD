---
layout: default
title: pltkdc
---

# pltkdc

Plotting of the vibronic coupling Hamiltonian potentials along 1D
normal mode cuts. Model potentials are plotted as lines with ab
initio reference data overlaid as points. Output is a gnuplot
script that is executed automatically.

## Usage

    pltkdc.x -f <paramfile> -m <mode> [options]

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
    -eps             Also write EPS output
    -h               Print help and exit

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

Gnuplot scripts are written to the current directory:

    plot_q<m>.gnu                              1D cut
    plot_q<m>_q<m2>.gnu                        diagonal cut
    plot_q<m>_s<S1>_s<S2>_diabcp.gnu           diabatic coupling
    plot_q<m>_s<S1>_s<S2>_dip.gnu              diabatic dipole

EPS files with the same base name are written when `-eps` is used.

## See also

[pltkdc.py](pltkdc_py) -- Python replacement using matplotlib.

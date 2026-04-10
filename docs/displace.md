---
layout: default
title: displace
---

# displace

Creation of displaced geometries along mass- and frequency-weighted
normal modes for the determination of the parameters of a vibronic
coupling Hamiltonian.

## Usage

    displace.x <input_file>

The `.inp` extension is appended automatically if not given. A log
file with the same base name and a `.log` extension is created
alongside the input file.

## Input file format

The input file uses a keyword-based format. Keywords and values are
separated by whitespace. Lines beginning with `#` are treated as
comments.

### Required keywords

**$freqfile**

Quantum chemistry output file containing the normal modes and
vibrational frequencies:

    $freqfile = freq.out

The file format is detected automatically. Supported formats:

| Format             | Detection                                  |
|--------------------|--------------------------------------------|
| Gaussian           | "Entering Gaussian System" in first line   |
| CFOUR              | Directory containing a file named `FCM`    |
| Hessian            | "Hessian" in first line                    |
| Turbomole (aoforce)| Recognised from aoforce output format      |
| ORCA               | ORCA Hessian file format                   |

**$cut**

Type of cut, step size, and number of points:

    $cut = <type> , <stepsize> , <npoints>

`type` is one of:

- `1mode` (or `all_1d`) -- 1D cuts along each normal mode.
- `2mode` -- diagonal 2D cuts for all symmetry-allowed bilinear
  coupling terms. Requires `$point_group`.
- `2mode_ondiag` (or `gamma_diag_2d`) -- diagonal 2D cuts for
  on-diagonal bilinear coupling terms only. Requires `$point_group`.

`stepsize` is the displacement increment in dimensionless normal
coordinates.

`npoints` is the number of displacement points in each direction
(positive and negative) from the reference geometry.

Example: a 1D cut with step size 0.5 and 10 points in each
direction:

    $cut = 1mode , 0.5 , 10

### Keywords required for 2D cuts

The following keywords are required when `$cut` is `2mode` or
`2mode_ondiag`.

**$point_group**

Point group of the molecule:

    $point_group = c2v

Used to determine which bilinear coupling coefficients are non-zero
by symmetry, and therefore which 2D cuts are needed.

Supported point groups: C1, Cs, Ci, C2, C2v, C2h, D2, D2h.

**$state_sym**

Symmetry labels of the electronic states:

    $state_sym
    A1 1
    B2 2
    A1 3
    $end

Each line contains a symmetry label followed by a state index.

## Output

All displaced geometries are written to a `geoms/` directory
(created automatically, overwriting any existing contents).

**Reference geometry:**

    geoms/q0.xyz

**1D cuts along mode *n*:**

    geoms/q<n>r.xyz     positive displacements (Q = 0, dq, 2*dq, ...)
    geoms/q<n>l.xyz     negative displacements (Q = 0, -dq, -2*dq, ...)

Each file is a concatenated XYZ file containing `npoints + 1`
structures.

**2D diagonal cuts along modes *n1* and *n2* (n1 < n2):**

    geoms/q<n1>r_q<n2>r.xyz     both modes displaced positively
    geoms/q<n1>l_q<n2>l.xyz     both modes displaced negatively

For diagonal cuts the two normal coordinates are each displaced by
the grid value divided by sqrt(2).

All coordinates are in Angstrom.

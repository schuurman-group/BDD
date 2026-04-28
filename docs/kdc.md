---
layout: default
title: kdc
---

# kdc

Calculation of the parameters of a vibronic coupling Hamiltonian.
The coupling coefficients are determined by fitting to ab initio
diabatic potential values, and an MCTDH or MultiQD operator file is
written.

## Usage

    kdc.x <input_file>

The `.inp` extension is appended automatically if not given. The
following output files are created using the input file base name:

- `.log` -- log file with fit information and parameter analysis
- `.op` or `.sop` -- MCTDH or MultiQD operator file
- `.dat` -- binary parameter file (used by pltkdc for plotting)

## Input file format

The input file uses a keyword-based format. Keywords begin with `$`
and values follow an `=` sign. Sections are terminated by `$end`.
Lines beginning with `#` are treated as comments.

### Required keywords

**$freqfile**

Quantum chemistry output file containing the normal modes and
vibrational frequencies:

    $freqfile = freq.out

The file format is detected automatically (Gaussian, CFOUR,
Turbomole, ORCA, or Hessian file).

**$point_group**

Point group of the molecule:

    $point_group = c2v

Used to determine which coupling coefficients are zero by symmetry.

Supported point groups: C1, Cs, Ci, C2, C2v, C2h, D2, D2h.

**$q0_ener**

Adiabatic energies at the reference geometry (in Hartree):

    $q0_ener
    -78.267243 1
    -77.999755 2
    -77.966886 3
    $end

Each line contains an energy value followed by a state index.

**$state_sym**

Symmetry labels of the electronic states:

    $state_sym
    A1 1
    B1 2
    A1 3
    $end

Each line contains a symmetry label followed by a state index.

**$bdfiles**

Blockdiag or GRaCI output files containing the diabatic potential
values. Files can be listed directly:

    $bdfiles
    mode1.out
    mode2.out
    $end

or read from a set file:

    $bdfiles = setfile.txt

where the set file containes the list of files to be read.

When using a set file, individual ab initio points can be excluded
from the fit using the `rm` keyword:

    mode1.out rm 5 10 15
    mode2.out
    mode3.out rm 2

### Optional keywords

**$order**

Order of the one-mode Taylor expansions (default: 6):

    $order = 4

**$opfile**

Operator file format. `mctdh` (default) or `multiqd`:

    $opfile = multiqd

MCTDH format produces a `.op` file. MultiQD format produces a
`.sop` file.

**$opstates**

Electronic states to include in the operator file (MultiQD format
only). By default all states are included:

    $opstates = 1 , 2 , 3

**$weight**

Enable weighted normal equation fitting. An optional scaling factor
can be given (default: 1.0):

    $weight = 0.5

**$shifts**

Zeroth-order shifts to the diabatic potential matrix elements (in
eV):

    $shifts
    0.1   1  2
    -0.05  2  3
    $end

Each line contains a shift value followed by two state indices. The
shift is applied symmetrically.

**$reexpand**

After the standard fit completes, re-expand all diabatic potential
matrix elements about the minimum of an adiabatic potential surface
v_s(Q), in normal-mode coordinates, instead of about the reference
geometry Q = 0. By default the minimum of adiabatic state 1 is used:

    $reexpand

To re-expand about the minimum of a different adiabatic state, give
the state index as an argument:

    $reexpand = 1

The re-expansion is exact and analytical: the model is a polynomial
of finite total degree, so the change of expansion centre stays
within the same monomial basis (no re-fitting). The bilinear
coefficients are unchanged; the 1-mode coefficients and the
diagonal energies are updated; new off-diagonal zeroth-order
constants `eps_<s2>_<s1>` are emitted into the operator file because
W_{ij}(Qmin) is generally non-zero for i≠j.

The minimum of the *adiabatic* v_s(Q) is the standard physical
ground-state geometry. The on-diagonal 1st-order coupling
coefficients of state s vanish at the new origin only when the
diabatic and adiabatic representations coincide at Qmin (e.g. when
state s is decoupled from all others by symmetry, as enforced by
`$blockdiag = {s}, {...}`). When state s mixes with others at Qmin,
the on-diagonal `tau1` coefficients reflect the gradient of the
diabatic W_{ss}(Q) at the adiabatic minimum, which is generally
non-zero.

A diagnostic block is written to the log file after the
re-expansion, reporting the diabatic and adiabatic gradient norms at
Qmin, the maximum residual coupling W_{s,j}(Qmin), the adiabatic
eigenvalues at Qmin, and a comparison `||Qmin_diab − Qmin_adiab||_∞`.
This allows the user to verify that the optimisation converged to
the expected minimum.

Notes:

- Qmin is located via a BFGS optimisation starting from Q=0 on the
  requested adiabatic surface, with a Hellmann-Feynman analytical
  gradient. The line search tolerates near-converged numerical
  noise from `dsyev` (treats `|grad| < 100 × tol_g` as converged
  even on Armijo failure). If convergence fails outright the
  original (Q=0) expansion is kept and a warning is written to the
  log file.
- The Cartesian geometry corresponding to Qmin is written to
  `reexpand_qmin.xyz` in the working directory.
- `$reexpand` is currently incompatible with diabatic dipole fitting
  (`$dip_sym` / `$q0_dipole`); attempting to combine them is rejected
  at input parsing.

**$blockdiag**

Block diagonalise the diabatic potential matrix. States are grouped
into blocks using braces:

    $blockdiag = {1, 2}, {3, 4, 5}

States within different blocks are decoupled by a unitary
transformation applied at each geometry. See: L. S. Cederbaum,
J. Schirmer, and H.-D. Meyer, *Block diagonalisation of Hermitian
matrices*, J. Phys. A: Math. Gen. **22**, 2427 (1989).

**$blockdiag_algorithm**

Algorithm used for block diagonalisation. `svd` (default) or
`invsqrt`:

    $blockdiag_algorithm = svd

- `svd` -- singular value decomposition (more numerically robust)
- `invsqrt` -- inverse square root method

**$cartgrad**

Write gradient and non-adiabatic coupling vectors to XYZ files for
visualisation:

    $cartgrad

No argument is needed.

### Dipole-related keywords

The following keywords enable fitting of diabatic dipole matrix
element expansions. Both `$dip_sym` and `$q0_dipole` must be given.

**$dip_sym**

Symmetry labels of the dipole operator components:

    $dip_sym
    B1 x
    B1 y
    A1 z
    $end

**$q0_dipole**

Diabatic dipole matrix elements at the reference geometry (in
atomic units). Each line contains the x, y, z components followed
by two state indices:

    $q0_dipole
    0.0  0.0  0.5  1  1
    0.1  0.2  0.3  1  2
    $end

## Output

### Operator file

The operator file (`.op` or `.sop`) contains the vibronic coupling
Hamiltonian in a format readable by MCTDH or MultiQD. It includes:

- Normal mode frequencies
- Vertical excitation energies
- One-mode coupling coefficients (up to the specified order)
- Two-mode (bilinear) coupling coefficients
- Dipole expansion coefficients (if dipole fitting is enabled)

Only symmetry-allowed parameters are written. Parameters with
magnitude below 5 x 10^-4 eV are omitted.

### Log file

The log file contains:

- Effective frequencies for each mode and state
- Frequency-weighted first-order coupling coefficients
- Spectroscopic importance ranking of normal modes
- Fit RMSD values
- Warnings for any non-zero parameters that are forbidden by
  symmetry

### Binary data file

The binary `.dat` file contains all system parameters and ab initio
data in Fortran unformatted format. It is read by the pltkdc
program for plotting model potentials against the ab initio
reference data.

## Example input

    $freqfile = freq.out

    $point_group = c2v

    $order = 4

    $q0_ener
    -78.267243 1
    -77.999755 2
    -77.966886 3
    $end

    $state_sym
    A1 1
    B1 2
    A1 3
    $end

    $bdfiles
    q1.out
    q2.out
    q3.out
    $end

    $opfile = mctdh

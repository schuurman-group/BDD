---
layout: default
title: mergeop
---

# mergeop

Merging of two kdc-generated Hamiltonian models into a single
MCTDH operator file. The two sets of electronic states are combined
into a block-diagonal Hamiltonian with energies referenced to the
global minimum.

## Usage

    mergeop.x <paramfile_A> <paramfile_B>

Both arguments are binary parameter files (`.dat`) produced by
separate kdc calculations.

## Input requirements

The two parameter files must be compatible:

- Same number of normal modes
- Same number of atoms and Cartesian coordinates
- Same one-mode expansion order
- Same normal mode frequencies
- Same diabatic dipole fitting flag (both on or both off)

The program exits with an error if any of these checks fail.

## Merging strategy

The merged Hamiltonian has N_A + N_B electronic states, where N_A
and N_B are the numbers of states in each input file.

- Coupling coefficients from file A occupy the (1:N_A, 1:N_A)
  block.
- Coupling coefficients from file B occupy the
  (N_A+1:N_A+N_B, N_A+1:N_A+N_B) block.
- All cross-coupling elements between the two sets of states are
  zero.
- Vertical excitation energies are shifted so that the global
  minimum across both sets of states is the energy zero.
- Dipole expansion coefficients (if present) are merged in the
  same block-diagonal fashion.

## Output

    merged.op

An MCTDH operator file written to the current directory containing
the complete merged Hamiltonian.

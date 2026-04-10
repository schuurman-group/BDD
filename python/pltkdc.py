#!/usr/bin/env python3
"""
pltkdc.py - Plot vibronic coupling Hamiltonian potentials

Python replacement for the Fortran pltkdc program, using matplotlib
for plotting and f2py-compiled Fortran routines for potential
evaluation.

Usage:
    python pltkdc.py -f FILE -m MODE [options]

Options mirror those of the Fortran pltkdc program.
"""

import argparse
import sys
import numpy as np

try:
    import bdd_potlib as potlib
except ImportError:
    sys.exit(
        "Error: bdd_potlib module not found.\n"
        "Build it with: ./build_potlib.sh (in the python/ directory)"
    )

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def parse_args():
    p = argparse.ArgumentParser(
        description="Plot model Hamiltonian potentials from KDC output"
    )
    p.add_argument("-f", required=True, metavar="FILE",
                   help="Binary parameter file (.dat)")
    p.add_argument("-m", required=True, type=int, metavar="MODE",
                   help="Plotting mode (normal mode index)")
    p.add_argument("-m2", type=int, default=None, metavar="MODE2",
                   help="Second mode for diagonal cut")
    p.add_argument("-xrange", nargs=2, type=float, default=[-7.0, 7.0],
                   metavar=("QI", "QF"),
                   help="Coordinate range (default: -7 7)")
    p.add_argument("-yrange", nargs=2, type=float, default=None,
                   metavar=("EI", "EF"),
                   help="Energy/function range")
    p.add_argument("-npnts", type=int, default=1000,
                   help="Number of grid points (default: 1000)")
    p.add_argument("-si", type=int, default=None,
                   help="Index of lowest state to plot")
    p.add_argument("-sf", type=int, default=None,
                   help="Index of highest state to plot")
    p.add_argument("-adiab", action="store_const", const=1, dest="surftyp",
                   help="Plot adiabatic potentials (default)")
    p.add_argument("-diab", action="store_const", const=2, dest="surftyp",
                   help="Plot diabatic potentials")
    p.add_argument("-dip", nargs=2, type=int, default=None,
                   metavar=("S1", "S2"),
                   help="Plot diabatic dipole for states S1, S2")
    p.add_argument("-diabcp", nargs=2, type=int, default=None,
                   metavar=("S1", "S2"),
                   help="Plot diabatic coupling for states S1, S2")
    p.add_argument("-pdf", default=None, metavar="FILE",
                   help="Output PDF file name (default: auto)")
    p.set_defaults(surftyp=1)
    return p.parse_args()


def main():
    args = parse_args()

    # Determine surface type from mutually exclusive options
    surftyp = args.surftyp
    dipsta1, dipsta2 = 0, 0
    dcpsta1, dcpsta2 = 0, 0

    if args.dip is not None:
        surftyp = 3
        dipsta1, dipsta2 = args.dip
    if args.diabcp is not None:
        surftyp = 4
        dcpsta1, dcpsta2 = args.diabcp

    # Diagonal cut flag
    ldiag = args.m2 is not None
    mode = args.m
    mode2 = args.m2 if ldiag else 0

    qi, qf = args.xrange
    npnts = args.npnts

    # Initialise the Fortran library
    potlib.init(args.f)

    # Get system dimensions
    nmodes, nsta, ndat, ldip = potlib.get_dimensions()

    # Compute model surfaces
    if surftyp == 1:
        # Adiabatic potentials
        if ldiag:
            qgrid, vpot = potlib.calc_adiab_diag(
                mode, mode2, nsta, npnts, qi, qf)
        else:
            qgrid, vpot = potlib.calc_adiab_1d(
                mode, nsta, npnts, qi, qf)

    elif surftyp == 2:
        # Diabatic potentials
        if ldiag:
            qgrid, vpot = potlib.calc_diab_diag(
                mode, mode2, nsta, npnts, qi, qf)
        else:
            qgrid, vpot = potlib.calc_diab_1d(
                mode, nsta, npnts, qi, qf)

    elif surftyp == 3:
        # Diabatic dipoles
        if ldiag:
            qgrid, vdip = potlib.calc_dip_diag(
                mode, mode2, dipsta1, dipsta2, npnts, qi, qf)
        else:
            qgrid, vdip = potlib.calc_dip_1d(
                mode, dipsta1, dipsta2, npnts, qi, qf)

    elif surftyp == 4:
        # Diabatic couplings
        if ldiag:
            qgrid, vcoup = potlib.calc_diabcp_diag(
                mode, mode2, dcpsta1, dcpsta2, npnts, qi, qf)
        else:
            qgrid, vcoup = potlib.calc_diabcp_1d(
                mode, dcpsta1, dcpsta2, npnts, qi, qf)

    # Get ab initio reference data
    nab, qab, vab = potlib.get_abinit_pot(
        mode, mode2 if ldiag else 0, ldiag, surftyp,
        dipsta1 if surftyp == 3 else dcpsta1,
        dipsta2 if surftyp == 3 else dcpsta2,
        ndat, nsta)

    # Trim ab initio arrays to actual number of points
    qab = qab[:nab]
    vab = vab[:, :, :nab]

    # State ranges
    si = args.si if args.si is not None else 1
    sf = args.sf if args.sf is not None else nsta

    # Build the x-axis label
    xlabel = f"$Q_{{{mode}}}$"
    if ldiag:
        xlabel = f"$Q_{{{mode}}}$ $Q_{{{mode2}}}$"

    # Determine output PDF filename
    if args.pdf is not None:
        pdfname = args.pdf
    else:
        pdfname = build_filename(mode, mode2, ldiag, surftyp,
                                 dipsta1, dipsta2, dcpsta1, dcpsta2)

    # Plot
    if surftyp == 1 or surftyp == 2:
        plot_potentials(qgrid, vpot, qab, vab, si, sf, nsta,
                        qi, qf, args.yrange, xlabel, surftyp,
                        pdfname)
    elif surftyp == 3:
        plot_dipoles(qgrid, vdip, qab, vab, dipsta1, dipsta2,
                     qi, qf, args.yrange, xlabel, pdfname,
                     ldip)
    elif surftyp == 4:
        plot_coupling(qgrid, vcoup, qab, vab, dcpsta1, dcpsta2,
                      qi, qf, args.yrange, xlabel, pdfname)

    print(f"Plot written to {pdfname}")


def build_filename(mode, mode2, ldiag, surftyp,
                   dipsta1, dipsta2, dcpsta1, dcpsta2):
    """Build default output filename."""
    base = f"plot_q{mode}"
    if ldiag:
        base += f"_q{mode2}"

    if surftyp == 3:
        base += f"_s{dipsta1}_s{dipsta2}_dip"
    elif surftyp == 4:
        base += f"_s{dcpsta1}_s{dcpsta2}_diabcp"

    return base + ".pdf"


def plot_potentials(qgrid, vpot, qab, vab, si, sf, nsta,
                    qi, qf, yrange, xlabel, surftyp, pdfname):
    """Plot adiabatic or diabatic potentials."""
    fig, ax = plt.subplots(figsize=(6, 6))

    # State index range (convert to 0-based for array slicing)
    si0 = si - 1
    sf0 = sf

    # Auto y-range
    vmin = np.min(vpot[:, si0:sf0])
    vmax = np.max(vpot[:, si0:sf0])

    if yrange is not None:
        ei, ef = yrange
    else:
        ei = min(-0.2, 0.98 * vmin)
        ef = 1.02 * vmax

    # Plot model potentials (lines)
    for s in range(si0, sf0):
        ax.plot(qgrid, vpot[:, s], linewidth=3)

    # Plot ab initio data (points)
    for s in range(si0, sf0):
        ax.plot(qab, vab[s, s, :], "o", color="black",
                markersize=3, zorder=5)

    ax.set_xlim(qi, qf)
    ax.set_ylim(ei, ef)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Energy (eV)")

    label = "Adiabatic" if surftyp == 1 else "Diabatic"
    ax.set_title(f"{label} potentials")

    fig.tight_layout()

    plt.show()

    with PdfPages(pdfname) as pdf:
        pdf.savefig(fig)

    plt.close(fig)


def plot_coupling(qgrid, vcoup, qab, vab, s1, s2,
                  qi, qf, yrange, xlabel, pdfname):
    """Plot a diabatic coupling element."""
    fig, ax = plt.subplots(figsize=(6, 6))

    if yrange is not None:
        ei, ef = yrange
    else:
        vmin = np.min(vcoup)
        vmax = np.max(vcoup)
        ei = min(-0.1, 0.98 * vmin)
        ef = 1.02 * vmax

    # Model coupling (line)
    ax.plot(qgrid, vcoup, linewidth=3)

    # Ab initio data (points)
    ax.plot(qab, vab[s1 - 1, s2 - 1, :], "o", color="black",
            markersize=3, zorder=5)

    ax.set_xlim(qi, qf)
    ax.set_ylim(ei, ef)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Energy (eV)")
    ax.set_title(f"Diabatic coupling W({s1},{s2})")

    fig.tight_layout()

    plt.show()

    with PdfPages(pdfname) as pdf:
        pdf.savefig(fig)

    plt.close(fig)


def plot_dipoles(qgrid, vdip, qab, vab, s1, s2,
                 qi, qf, yrange, xlabel, pdfname, ldip):
    """Plot diabatic dipole surfaces."""
    fig, ax = plt.subplots(figsize=(6, 6))

    comp_labels = ["x", "y", "z"]

    if yrange is not None:
        ei, ef = yrange
    else:
        vmin = np.min(vdip[:, :3])
        vmax = np.max(vdip[:, :3])
        ei = min(-0.1, 0.98 * vmin)
        ef = 1.02 * vmax

    # Model dipole surfaces (lines)
    for c in range(3):
        ax.plot(qgrid, vdip[:, c], linewidth=3, label=comp_labels[c])

    # Ab initio data (points)
    if ldip:
        for c in range(3):
            ax.plot(qab, vab[c, c, :], "o", color="black",
                    markersize=3, zorder=5)

    ax.set_xlim(qi, qf)
    ax.set_ylim(ei, ef)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("D (a.u.)")
    ax.set_title(f"Diabatic dipole ({s1},{s2})")
    ax.legend()

    fig.tight_layout()

    plt.show()

    with PdfPages(pdfname) as pdf:
        pdf.savefig(fig)

    plt.close(fig)


if __name__ == "__main__":
    main()

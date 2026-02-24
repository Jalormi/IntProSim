# pulsar_sim/cli.py
import os
import sys
import numpy as np
import argparse
import psrchive

from pulsar_sim.io import readparfile, write_samplefile
from pulsar_sim.dm_series import make_dmseries, make_dips
from pulsar_sim.scat_series import make_scatseries
from pulsar_sim.fits_writer import makefits


def parse_arguments():
    parser = argparse.ArgumentParser(description="Code to create pulsar fits files")
    parser.add_argument(
        "-E", "--ephem", default="Downloads/B1937+21.txt", help="Ephemeris file"
    )
    parser.add_argument("-Ms", "--mjd_start", type=float, default=58900.0)
    parser.add_argument("-Mf", "--mjd_finish", type=float, default=59000.0)
    parser.add_argument("-fl", "--freq_low", type=float, default=1300.0)
    parser.add_argument("-fh", "--freq_high", type=float, default=1500.0)
    parser.add_argument("-nc", "--nchan", type=int, default=128)
    parser.add_argument("-nb", "--nbin", type=int, default=1024)
    parser.add_argument("-ns", "--nsub", type=int, default=1)
    parser.add_argument("-ts", "--tsub", type=float, default=10.0)
    parser.add_argument("-N", "--nobs", type=int, default=10)
    parser.add_argument("-np", "--npol", type=int, default=1)
    parser.add_argument("-snr", "--snr", type=float, default=10)
    parser.add_argument("-s", "--spidx", type=float, default=-1.7)
    parser.add_argument("-fd", "--fluxd", type=float, default=1)
    parser.add_argument("-pbf", "--pbf", type=int, default=1)
    parser.add_argument("-i", "--simpfile", default=None)
    parser.add_argument("-dmf", "--dmfile", default=None)
    parser.add_argument("-noscat", "--noscat", action="store_true", default=False)
    parser.add_argument(
        "--nodm", action="store_true", default=False, help="Do not apply DM smearing"
    )
    parser.add_argument(
        "--onlydm",
        action="store_true",
        default=False,
        help="Apply only DM smearing (no scattering)",
    )
    parser.add_argument(
        "--onlyscat",
        action="store_true",
        default=False,
        help="Apply only scattering (no DM smearing)",
    )
    parser.add_argument("-v", "--verbose", action="store_true", default=True)
    return parser


def main():
    args = parse_arguments().parse_args()

    verbose = args.verbose
    parfile = args.ephem

    if verbose:
        print("Reading parameter file...")

    ar_psr, _, _, _ = readparfile(parfile)

    if verbose:
        print(f"Pulsar name: {ar_psr}")

    simpfile = args.simpfile or "Downloads/simpulse_input.txt"

    if not os.path.isfile(simpfile):
        print("\nError: Input simpulse file missing. Creating sample then exiting.")
        write_samplefile()
        exit(0)

    if args.ephem == None:
        print("\nError: Parameter file is not given. Exiting.\n")
        exit(0)

    if args.simpfile == None:
        simpfile = (
            "simulation_results/result_simpulse_input/J1713_simpulse_input_gm.txt"
        )
    else:
        simpfile = args.simpfile

    if not os.path.isfile(str(simpfile)):
        print("\nError: Input file is not given.")
        print("A sample file simpulse_input_J1939.txt is written out for reference.")
        print("Exiting.")
        write_samplefile()
        exit(0)

    # ---- Prevent conflicting switches ----
    if args.onlydm and args.onlyscat:
        print("ERROR: Cannot use --onlydm and --onlyscat together.")
        sys.exit()

    nsub = args.nsub
    tsub = args.tsub
    npol = args.npol
    nchan = args.nchan
    nbin = args.nbin
    freq_low = args.freq_low
    freq_high = args.freq_high
    epoch_start = args.mjd_start
    epoch_stop = args.mjd_finish
    nobs = args.nobs
    snr = args.snr
    spidx = args.spidx
    fluxd = args.fluxd
    pbf = args.pbf
    onlyscat = args.onlydm
    onlydm = args.onlydm
    nodm = args.nodm
    noscat = args.noscat
    if not args.noscat:
        if isinstance(pbf, int):
            if pbf > 4:
                print("Given PBF number %d doesn't exist! Exiting." % pbf)
                exit(0)
    # Reads the pulsar name from parameter file
    ar_psr, _, _, _ = readparfile(parfile)

    # Removes the simulation parameters out file, if already exists
    if os.path.isfile(ar_psr + "_dmalptscat.dm"):
        os.remove(ar_psr + "_dmalptscat.dm")

    # Creates DM timeseries or read from input file
    if args.dmfile == None:

        # Generates the epochs for simulation
        epochs = np.linspace(epoch_start, epoch_stop, nobs, endpoint=True, dtype=float)
        dDM = np.zeros(nobs)

        # Creates the DM timeseries
        dDM = make_dmseries(nobs, simpfile)

        # Creates tscat and alpha series
        tscat = np.zeros(nobs)
        alpha = np.zeros(nobs)

    else:
        # Reads epochs and delta DMs from the input file
        epochs = np.loadtxt(args.dmfile, usecols=0, dtype=float, comments="#", ndmin=1)
        dDM = np.loadtxt(args.dmfile, usecols=1, dtype=float, comments="#", ndmin=1)
        nobs = np.size(dDM)

        # Creates tscat and alpha series
        tscat = np.zeros(nobs)
        alpha = np.zeros(nobs)

    tscat, alpha = make_scatseries(dDM, simpfile)
    dmevents = make_dips(epochs, dDM, ar_psr, simpfile)
    dDM += dmevents
    # Runs the makefits() to generate simulated archives
    for obs in range(nobs):
        makefits(
            parfile,
            nsub,
            tsub,
            npol,
            nchan,
            nbin,
            freq_low,
            freq_high,
            snr,
            epochs[obs],
            dDM[obs],
            tscat[obs],
            alpha[obs],
            spidx,
            fluxd,
            pbf,
            simpfile,
            args.noscat,
            args.nodm,
            args.onlydm,
            args.onlyscat,
            verbose=True,
        )

    # for printing ephemeris file, MJD range, frequency range, no. of observations
    if verbose:
        print("Starting pulsar simulation...")
        print(f"Ephemeris file: {args.ephem}")
        print(f"MJD range: {args.mjd_start} – {args.mjd_finish}")
        print(f"Frequency range: {args.freq_low} – {args.freq_high} MHz")
        print(f"Number of observations: {args.nobs}")
    # for reading the ephemeris file and creating DM series with a set of observations
    if verbose:
        print(f"Reading ephemeris from {parfile}")
    if verbose:
        print(
            f"Creating DM series with {nobs} observations between {epoch_start} and {epoch_stop}"
        )


if __name__ == "__main__":
    main()

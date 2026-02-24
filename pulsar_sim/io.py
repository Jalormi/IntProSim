import numpy as np
import psrchive
from numpy import exp
from astropy.coordinates import SkyCoord
import astropy.units as u


def readparfile(parfile):
    """
    Reads the input par file and collects pulsar name (PSRJ), right ascension (RAJ),
    Declination (DECJ), and dispersion measure (DM).
    """
    ar_psr = None
    ar_raj = None
    ar_decj = None
    ar_dm = 0.0
    elong = 0.0
    elat = 0.0

    with open(parfile, "r") as parf:
        for line in parf:
            line = line.split()
            if len(line) < 2:
                continue
            if line[0] == "PSRJ":
                ar_psr = line[1]
            elif line[0] == "RAJ":
                ar_raj = line[1]
            elif line[0] == "DECJ":
                ar_decj = line[1]
            elif line[0] == "DM":
                ar_dm = float(line[1])

    # If RAJ/DECJ not found, try ELONG/ELAT conversion
    if ar_raj is None or ar_decj is None:
        with open(parfile, "r") as parf:
            for line in parf:
                line = line.split()
                if len(line) < 2:
                    continue
                if line[0] == "ELONG":
                    elong = float(line[1])
                elif line[0] == "ELAT":
                    elat = float(line[1])
        from astropy.coordinates import SkyCoord
        import astropy.units as u

        coo = SkyCoord(
            lon=elong * u.degree, lat=elat * u.degree, frame="geocentrictrueecliptic"
        )
        coo_tmp = coo.transform_to("fk5")
        coo_conv = coo_tmp.to_string("hmsdms", sep=":")
        ar_raj = str(coo_conv.split(" ")[0])
        ar_decj = str(coo_conv.split(" ")[1])

    if ar_psr is None or ar_raj is None or ar_decj is None:
        print("\nError: Pulsar name and/or position are not given. Exiting.\n")
        exit(0)

    return ar_psr, ar_raj, ar_decj, ar_dm


def write_samplefile():
    f = open("simpulse_input.txt", "a")
    f.write("# all the values are referenced to 1 GHz\n")
    f.write(
        "# gauss_no. / phase(0-1) / amplitude / w50(0-0.5) / Pscaling / Hscaling / Wscaling\n"
    )
    f.write("gauss1 0.5 1.0 0.05 0.0 0.0 0.0\n")
    f.write(
        "# lorentz_no. / phase (0-1) / amplitude / gamma (width) / Xscaling / Ascaling / Gscaling\n"
    )
    f.write("#lorentz1 0.1 1.0 0.03 0.0 0.0 0.0\n\n")
    f.write("# scattering time reference\n")
    f.write("# scat_time (ms) / alpha / rms / mode / param\n")
    f.write("scat 0.005 -4.4 0.4 stable tscat\n")
    f.write(
        '# mode can be "stable" - same alpha/tscat for all epochs; "dm" - same variation as DM timeseries in alpha/tscat;\n'
    )
    f.write(
        '# "random" - vary alpha/tscat randomly with the given rms; and "powerlaw" - Uses the given alpha to create a tscat series\n'
    )
    f.write(
        '# param can be "tscat", "alpha" or "none". Mode powerlaw will only work with "tscat"\n\n'
    )
    f.write("# DM spectral model\n")
    f.write("#     amplitude   spectral_index\n")
    f.write("dmsp    1e-6        -2.0\n")
    f.close()

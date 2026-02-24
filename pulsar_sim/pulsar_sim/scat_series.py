import numpy as np
from numpy import exp
from pulsar_sim.utils import make_timeseries


# Generates the scattering timeseris.
def make_scatseries(dmseries, simpfile):
    """
    Creates a scattering time series based on the arguments from the
    input file. The line should be in the format
    "scat <tscat> <alpha> <rms> <mode> <param>". The value of tscat
    should be at 1 GHz reference. The value of alpha should be negative.
    mode can be
    "stable" - same alpha/tscat for all epochs;
    "dm" - same variation as DM timeseries in alpha/tscat;
    "random" - vary alpha/tscat randomly with the given rms
    "powerlaw" - Use the given alpha to create a tscat series
    param can be "tscat", "alpha" or "none". Mode powerlaw will only
    work with "tscat". The variables are read as positional arguments,
    so be careful with that.
    """

    # Reads the different scattering variables from the input file.
    tsc_ref = []
    alp_ref = []
    param_rms = []
    mode = []
    param_ref = []

    for line in open(simpfile):

        if line.startswith("scat"):
            data = line.split()
            tsc_ref.append(float(data[1]))
            alp_ref.append(float(data[2]))
            param_rms.append(float(data[3]))
            mode.append(data[4])
            param_ref.append(data[5])
    tsc_ref = np.asarray(tsc_ref, dtype=float)
    alp_ref = np.asarray(alp_ref, dtype=float)
    param_rms = np.asarray(param_rms, dtype=float)

    # Creates a timeseries in tscat and alpha same as dm variation.
    # Only one of those can be varied w.r.t dm, not both.
    if mode[0] == "dm":
        medDM = np.median(dmseries)

        if param_ref[0] == "tscat":
            tt = dmseries / np.sum(dmseries)
            tt = tt * (np.std(tt) * (1.0 / np.max(tt)))
            tscat = tsc_ref + (tsc_ref * tt / 10.0)
            tscat = tsc_ref + (dmseries / medDM) / 100.0
            tscat += np.min(tscat) * -1.0
            alpha = alp_ref * np.ones(np.size(dmseries), dtype=float)

        if param_ref[0] == "alpha":
            tt = dmseries / np.sum(dmseries)
            tt = tt * (np.std(tt) * (1.0 / np.max(tt)))
            alpha = alp_ref + (alp_ref * tt / 100.0)
            tscat = tsc_ref * np.ones(np.size(dmseries), dtype=float)

        return (tscat, alpha)

    # Creates a timeseries of tscat and alpha with the same input values.
    if mode[0] == "stable":
        tscat = tsc_ref * np.ones(np.size(dmseries), dtype=float)
        alpha = alp_ref * np.ones(np.size(dmseries), dtype=float)

        return (tscat, alpha)

    # Creates a tscat or dm variable timeseries with a given rms
    if mode[0] == "random":

        if param_ref[0] == "none" and param_ref != "tscat" and param_ref != "alpha":
            print(
                '\nERROR: Wrong combination of mode "%s" with param "%s".'
                % (mode[0], param_ref[0])
            )
            print("Change the variables accordingly to re-run simulator.\n")
            exit(0)

        if param_ref[0] == "tscat":

            if (tsc_ref - param_rms * 5.0) < 0.0:
                param_rms = tsc_ref / 5.0
                print(
                    "Warning: RMS tscat makes it go negative. Changing it to %f"
                    % (param_rms)
                )
            tscat = np.random.normal(tsc_ref, param_rms, np.size(dmseries))
            alpha = alp_ref * np.ones(np.size(dmseries), dtype=float)

            return (tscat, alpha)

        if param_ref[0] == "alpha":
            tscat = tsc_ref * np.ones(np.size(dmseries), dtype=float)
            alpha = np.random.normal(alp_ref, param_rms, np.size(dmseries))
            # alpha = alp_ref * np.random.rand(alp_ref, param_rms, np.size(dmseries))

            return (tscat, alpha)

    # Creates a tscat timeseries based on the powerlaw
    if mode[0] == "powerlaw":
        tscat = make_timeseries(1e-8, alp_ref, np.size(dmseries))
        tscat += np.min(tscat) * -1
        tscat /= 10.0
        tscat += tsc_ref
        alpha = alp_ref * np.ones(np.size(dmseries), dtype=float)

        return (tscat, alpha)

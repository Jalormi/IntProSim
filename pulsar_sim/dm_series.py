# Creates DM series based on the input
import numpy as np
from numpy import exp
from pulsar_sim.utils import make_timeseries


def make_dmseries(n, simpfile):
    """
    Creates a DM time series based on the amplitude and power given in dmsp
    line of the input file.
    """
    amp = 1e-6
    power = -2.0
    for line in open(simpfile):
        if line.startswith("dmsp"):
            data = line.split()
            try:
                amp = float(data[1])
            except:
                print("Amplitude is not given. Using 1e-6.")
            try:
                power = float(data[2])
            except:
                print("Spectral index is not given. Using -2.0")

    dmseries = np.zeros(n)

    dmseries = make_timeseries(amp, power, n)

    return dmseries


def make_dips(epochs, dmseries, psr, simpfile):
    epoch_dip = []
    amplitude = []
    rtime = []
    param = []
    for line in open(simpfile):
        if line.startswith("edip"):
            data = line.split()
            epoch_dip.append(float(data[1]))
            amplitude.append(float(data[2]))
            rtime.append(float(data[3]))
            param.append(str(data[4]))

    if not (len(epoch_dip) == 0):
        x = np.linspace(
            0, (epochs[-1] - epochs[0]), np.size(epochs), endpoint=True, dtype=float
        )
        tada = np.absolute(epochs - epoch_dip)
        dipx = tada.argmin()
        dip = np.zeros(np.size(epochs))
        if "gauss" in param[0]:
            dip = np.roll(amplitude * ((1.0 - np.exp(-x / rtime))), dipx)
            f = open(psr + "_amplitudes.txt", "w")
            for i in range(np.size(epochs)):
                f.write("%f %f\n" % (epochs[i], dip[i] - amplitude))
            f.close()
            return np.zeros(np.size(epochs))
        if param[0] == "dm":
            dip = np.roll(amplitude * ((1.0 - np.exp(-x / rtime))), dipx)
            return dip - amplitude
    if len(epoch_dip) == 0:
        return np.zeros(np.size(epochs))

import numpy as np
from numpy import exp


def bandshape(nchan):
    """
    A stupid way to create a bandshape. Could be improved
    to have a better bandshape so that the S/N across the
    band resembles that.
    """

    sqband = np.zeros(nchan, dtype=float)

    for i in range(nchan):

        if i > 0.02 * nchan:
            sqband[i] = 1

    chans = np.linspace(0.0, nchan, nchan, dtype=int)
    chans = chans.astype(float)
    phase = nchan / 2.0
    sigma = nchan * 0.1
    gauss = exp(-chans / 15.0)
    gauss = np.flip(gauss)
    bandshape = convolve_funcs(gauss, sqband)
    bandshape += np.random.normal(0, 0.05, nchan)
    bandshape /= np.median(bandshape)
    bandshape = (10**bandshape) * 10.0

    return bandshape


def convolve_funcs(func1, func2):
    """
    Algorithm implementing circular convolution of two functions.
    Currently using the FFT based convolution as it is way faster
    than the usual way.
    """

    # Checking the size of the arrays and zero padding if one of
    # them is shorter than the other.
    n1 = np.size(func1)
    n2 = np.size(func2)

    if n2 < n1:
        zeros = np.zeros(n1 - n2).tolist()
        func2 = func2.tolist() + zeros
        func2 = np.asarray(func2, dtype=float)

    if n1 < n2:
        zeros = np.zeros(n2 - n1).tolist()
        func1 = func1.tolist() + zeros
        func1 = np.asarray(func1, dtype=float)

    # Taking the fft of both the functions
    f1_fft = np.fft.fft(func1)
    f2_fft = np.fft.fft(func2)

    # Multiplying the fft of two functions
    m = f1_fft * f2_fft

    # Gives out the final convolved function
    y = abs(np.fft.ifft(m))

    return y


# Creates a power law spectrum.
def generate_spectrum(freqs, amp, beta):

    return amp * freqs**beta


# Creates a timeseries from the power spectrum.
def make_timeseries(amp, beta, ndays):
    """
    Creates a timeseries using the given spectral index (beta) and amplitude
    (amp). Random noise is also added to the spectrum.
    """
    nfreqs = 10000
    freqs = np.linspace(1e-3, 2, nfreqs)
    powersp = generate_spectrum(freqs, amp, beta) * np.random.rand(len(freqs))
    fouriersp = np.sqrt(powersp) * np.exp(1j * np.pi * np.random.rand(len(powersp)))
    Tseriesfft = np.concatenate(([0], fouriersp, np.conj(fouriersp[::-nfreqs])))
    Tseries = np.real(np.fft.ifft(Tseriesfft))[1 : nfreqs + 1]
    Tseries = Tseries[:-100]
    Tseries = Tseries[100:]
    freqs = freqs[:-100]
    freqs = freqs[100:]

    phase_points = np.linspace(1e-3, 1, ndays)
    Tseries = np.interp(phase_points, freqs, Tseries)
    return Tseries

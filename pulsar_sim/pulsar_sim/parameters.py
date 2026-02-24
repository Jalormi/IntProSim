"""
Default parameters for the pulsar_sim package.
"""

DEFAULTS = {
    "ephem": "Downloads/B1937+21.txt",
    "mjd_start": 58900.0,
    "mjd_finish": 59000.0,
    "freq_low": 1300.0,
    "freq_high": 1500.0,
    "nchan": 128,
    "nbin": 1024,
    "nsub": 1,
    "tsub": 10.0,
    "nobs": 10,
    "npol": 1,
    "snr": 10.0,
    "spidx": -1.7,
    "fluxd": 1.0,
    "pbf": 1,
    "simpfile": "Downloads/simpulse_input.txt",
    "dmfile": None,
    "verbose": True,
    "noscat": False,
    "nodm": False,
    "onlydm": False,
    "onlyscat": False,
}

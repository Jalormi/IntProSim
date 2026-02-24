# pulsar_sim/__init__.py
"""
pulsar_sim package — pulsar simulation tools.
"""

from .io import readparfile, write_samplefile
from .dm_series import make_dmseries, make_dips
from .scat_series import make_scatseries
from .profiles import make_profile
from .fits_writer import makefits
from .cli import parse_arguments, main

__all__ = [
    "readparfile",
    "write_samplefile",
    "make_dmseries",
    "make_dips",
    "make_scatseries",
    "make_profile",
    "makefits",
    "parse_arguments",
    "main",
]

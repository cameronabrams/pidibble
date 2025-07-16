# Author: Cameron F. Abrams <cfa22@drexel.edu>

from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("pidibble")
except PackageNotFoundError:
    __version__ = "unknown"

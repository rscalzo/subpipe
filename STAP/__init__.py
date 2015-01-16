from .STAP_unzip import unzip
from .STAP_WCS import WCS
from .STAP_remap import remap
from .STAP_SEx import SEx
from .STAP_hotpants import hotpants
from .STAP_flag import classify
from .STAP_xref import xref
from .STAP_roundint import roundint
from .STAP_flatcombine import flatcombine
from .STAP_maskcombine import maskcombine
from .STAP_zeropoint import zeropoint
from .STAP_unzip import unzip
from .STAP_maskcombine import maskcombine
from .STAP_makethumbs import makethumbs
# from .STAP_getDSS import getDSS
# from .STAP_getmag import getmag
__all__ = ["WCS", "remap", "SEx", "hotpants", "flag", "xref", "roundint",
           "flatcombine", "zeropoint", "unzip", "maskcombine","makethumbs"]

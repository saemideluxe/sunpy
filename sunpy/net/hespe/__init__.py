"""Access to the HESPE

The HESPE database contains a lot of pre-calculated maps from RHESSI.
The database is accessible at http://hsp.cs.technik.fhnw.ch/browser in the
browser. This package provides access with python to all the data on HESPE.

"""

#
# TODO
# - Add support for Spectra and Spectrogram classes
# - Cache-management?
#

from __future__ import absolute_import

__authors__ = ["Samuel Stachelski"]
__email__ = "samuel.vonstachelski@fhnw.ch"


from sunpy.net.hespe.hespe import get_flareevents_between
from sunpy.net.hespe.hespe import get_maps_of_flareevent
from sunpy.net.hespe.hespe import get_lightcurve_of_flareevent
from sunpy.net.hespe.hespe import group_maps_by_energy

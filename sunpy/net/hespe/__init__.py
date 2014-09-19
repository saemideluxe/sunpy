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

__authors__ = ["Samuel Stachelski"]
__email__ = "samuel.vonstachelski@fhnw.ch"

from __future__ import absolute_import

from sunpy.net.hespe.hespe import *

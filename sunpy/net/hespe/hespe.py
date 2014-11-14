# -*- coding: utf-8 -*-
# Author: Samuel Stachelski <samuel.vonstachelsi@fhnw.ch>

from __future__ import absolute_import
from __future__ import division

import os
import urllib
import json
import time
from datetime import datetime
import numpy
import pandas

import sunpy.map
from sunpy.cm import cm
from sunpy.lightcurve import RHESSISummaryLightCurve
import astropy.io.fits as fits


# URL-constants
_FLARE_DATA_BASE_URL = "http://hsp.cs.technik.fhnw.ch/webdbi-0.6-SNAPSHOT/getE\
ventAllData/"
_IMAGE_DATABASE_BASE_URL = "http://hsp.cs.technik.fhnw.ch/flaredata/"
_EVENT_LIST_BASE_URL = "http://hsp.cs.technik.fhnw.ch/webdbi-0.6-SNAPSHOT/filt\
eredEvents/"
_EVENT_LIST_URL = "http://hsp.cs.technik.fhnw.ch/webdbi-0.6-SNAPSHOT/filteredE\
vents/-1/-1/-1/-1/null/null/-1/0/0/startTime/DESC"

# HESPE uses it's own time system, which is seconds since 1.1.1979.
# To work with python's datetime, the HESPE-time is converted to UNIX-time.
# This constant is the number of seconds between 1.1.1970 and 1.1.1979.
_HESPE_TIME_TO_UNIX_TIME = 283996800

# keymaps for url/json-url access
ALG_KEYMAP = {"uvsmooth": "quicklook",
              "visclean": "quicklook_visclean",
              "memnjit": "quicklook_memnjit",
              "bpmap": "quicklook_bpmap", }
MAP_TYPE_KEYMAP = {"photon": "visBags",
                   "electron": "electron_visBags"}


def _get_flare_data(flare_id, ivs_type="coarse"):
    """Returns a dictionary with information (including URLs to maps)
    about the specified event.
    """
    flare_data_url = "%s%d" % (_FLARE_DATA_BASE_URL, flare_id)
    url_obj = urllib.urlopen(flare_data_url)
    if url_obj.code != 200:
        raise Exception("Address \"%s\" returned HTTP-Code %d" %
                        (flare_data_url, url_obj.code))
    data = json.load(url_obj)

    # return only the selected ivs_type
    found_ivs_type = False
    for i in data:
        if i["ivsType"] == ivs_type:
            data = i
            found_ivs_type = True
            break

    if not found_ivs_type:
        raise Exception("Unsupported ivsType: " + ivs_type)

    return data


def _rewrite_url_to_fits_url(url):
    # The returned event data from _get_flare_data() contains only URLs to the
    # *.png files, but the most time there exists also a *.fits file
    """Returns the given URL with ".png" replace by ".fits" at the end."""
    if not url.endswith(".png"):
        raise Exception("Expected a png-file: %s" % (url))
    return "".join([_IMAGE_DATABASE_BASE_URL, url[:-3], "fits"])


def _hespeTimeToUTC(time):
    """Converts a HESPE time stamp to a python datetime object (UTC). HESPE 
    uses internally the number of seconds since the 1.1.1979.

    Parameters
    ----------
    time : numeric type
        Time in seconds since 1.1.1979 (HESPE time stamp)
    """
    return datetime.utcfromtimestamp(time + _HESPE_TIME_TO_UNIX_TIME)


def _download_hespe_fits_as_map(url, energy_l=0, energy_h=0):
    """Downloads a FITS file (which should contain a map) from the given URL
    and converts it to a sunpy Map object.
    The downloaded FITS file is saved in the converted version on the local
    disk for caching.
    energy_l and energy_h can not be read from the FITS file itself. These
    values must be passed as arguments.
    """

    download_dir = sunpy.config.get('downloads', 'download_dir') + "/HESPE/"

    try:
        os.mkdir(download_dir)
    except OSError:
        pass

    filename = url[url.rfind("/") + 1:]
    path = download_dir + filename

    if os.path.exists(path):
        # file is cached on disk
        fits_img = fits.open(path)
    else:
        # file is not cached, download and converting necessary
        fits_img = fits.open(url)

        # extract seconds since 1.1.1979 from filename
        pos = url.find("_")
        _time_string = url[pos + 1:url.find("_", pos+1)]
        _time = float(_time_string)
        dt = _hespeTimeToUTC(_time)

        # set attributes for sunpy-compatability
        fits_img[1].header.append(("energy_l", energy_l))
        fits_img[1].header.append(("energy_h", energy_h))
        fits_img[1].header.append(("telescop", "RHESSI"))
        fits_img[1].header.append(("instrume", "RHESSI"))

        fits_img[1].header.append(("crval1", fits_img[1].data.XC[0]))
        fits_img[1].header.append(("crval2", fits_img[1].data.YC[0]))
        fits_img[1].header.append(("naxis1", fits_img[1].data[0][0].shape[0]))
        fits_img[1].header.append(("naxis2", fits_img[1].data[0][0].shape[1]))
        center_x = (fits_img[1].data[0][0].shape[0] + 1) / 2.0
        center_y = (fits_img[1].data[0][0].shape[1] + 1) / 2.0
        fits_img[1].header.append(("crpix1", center_x))
        fits_img[1].header.append(("crpix2", center_y))

        fits_img[1].header.append(("cdelt1", fits_img[1].data.DX[0]))
        fits_img[1].header.append(("cdelt2", fits_img[1].data.DY[0]))
        fits_img[1].header.append(("date-obs", str(dt)))
        fits_img[1].header.append(("ctype1", fits_img[1].data.XUNITS[0]))
        fits_img[1].header.append(("ctype2", fits_img[1].data.YUNITS[0]))
        fits_img[1].header.append(("crota1", fits_img[1].data.ROLL_ANGLE[0]))
        fits_img[1].header.append(("crota2", fits_img[1].data.ROLL_ANGLE[0]))

        fits_img.writeto(path)

    sunpy_map = sunpy.map.rhessi.RHESSIMap(fits_img[1].data[0][0],
                                           fits_img[1].header)
    sunpy_map.cmap = cm.cmlist.get('rhessi')

    return sunpy_map


def _generate_parameter_url(min_date=-1,
                            max_date=-1,
                            min_duration=-1,
                            max_duration=-1,
                            min_goes="null",
                            max_goes="null",
                            flare_id=-1,
                            start_nmbr=0,
                            total_nmbr=0,
                            sort_parameter="startTime",
                            sortDirection="DESC"):
    """HESPE provides a filter-URL to get a list of events which match the
    filter criteria. This function takes any of the provided filter arguments
    and returns the corresponding URL.
    """

    return "%s%s/%s/%s/%s/%s/%s/%s/%s/%s/%s/%s" % (
        _EVENT_LIST_BASE_URL, min_date, max_date, min_duration, max_duration,
        min_goes, max_goes, flare_id, start_nmbr, total_nmbr, sort_parameter,
        sortDirection
    )


def get_maps_of_flareevent(flare_id,
                      ivs_type="coarse",
                      ray_type="photon",
                      img_alg="visclean"):
    """Downloads all maps of a flare event.

    Parameters
    ----------
    flare_id : int
        HESPE-ID of the flare, from which the data should be obtained
    ivs_type : str
        The type of visualization with which the map was generated.  May be
        "coarse", "fixed" or "fine". This parameter describes it was integrated
        over the raw-data.
    ray_type : str
        The type of rays with which the map was generated. May be "photon" or
        "electron".  Electron maps are available for all events.
    img_alg : str
        The image algorithm with which the map was generated.  May be
        "uvsmooth", "visclean", "memnjit" or "bpmap".

    Returns
    ----------
    A list containing sunpy.map.rhessi.RHESSIMap objects for the given flare

    Examples
    --------
    >>> maps = get_maps_of_flare(13060315, ivs_type="fine")
    """

    flare_data = _get_flare_data(flare_id, ivs_type)
    img_alg = ALG_KEYMAP[img_alg]

    map_list = list()
    n = len(flare_data[MAP_TYPE_KEYMAP[ray_type]])
    c = 1
    for i in flare_data[MAP_TYPE_KEYMAP[ray_type]]:
        quicklook_url = i[img_alg]
        energy_l = i["lowerEnergy"]
        energy_h = i["upperEnergy"]
        fits_url = _rewrite_url_to_fits_url(quicklook_url)
        sunpy_map = _download_hespe_fits_as_map(fits_url, energy_l, energy_h)
        map_list.append(sunpy_map)
        print("downloaded %d of %d" % (c, n))
        c += 1

    return map_list


def get_lightcurve_of_flareevent(flare_id):
    """This returns a lightcurve object for the specified flare event"""
    flare_data = _get_flare_data(flare_id)

    # download lightcurves data
    url = _rewrite_url_to_fits_url(flare_data["lightcurves"])
    lc_fits = fits.open(url)
    time_line = map(_hespeTimeToUTC, lc_fits[1].data)

    # create list of labels
    labels = []
    for label in lc_fits[2].data:
        labels.append(str(label[0]) + " - " + str(label[1]) + " keV")

    # create lightcurve object with pandas DataFrame
    data_frame = pandas.DataFrame(lc_fits[0].data,
                                  index=time_line,
                                  columns=labels)
    ligthcurves = sunpy.lightcurve.RHESSISummaryLightCurve.create(data_frame)

    return ligthcurves


def get_flareevents_between(start, end):
    """Returns a list which contains the IDs of all flare events in the
    specified time range.

    Parameters
    ----------
    start : datetime
        Start time to filter from.
    end : datetime
        End time to filter to.

    Examples
    --------
    >>> list = get_flares_between(datetime(2013, 1, 1), datetime(2013, 3, 1))
    """

    # Convert to unix-time. URL needs time stamp in milliseconds
    start_unix = int(time.mktime(start.timetuple()) * 1000)
    end_unix = int(time.mktime(end.timetuple()) * 1000)

    server_filter_url = _generate_parameter_url(min_date=start_unix,
                                                max_date=end_unix)

    url_obj = urllib.urlopen(server_filter_url)
    event_list = json.load(url_obj)
    reduced_event_list = map(lambda x: x["flareId"], event_list)

    return reduced_event_list


def group_maps_by_energy(map_list):
    """Returns a dictionary with an entry for each energy level. The keys are 
    strings which describe the energy level. The values are CompositeMaps where
    as all maps of the same energy level are alpha blended over each other.
    TODO: make similar method and group by time-intervalls

    Parameters
    ----------
    map_list : iterable
        rhessi maps from a flare event
        
    Returns
    -------
    dict: (str: CompositeMap)

    """
    dic = {}
    
    for i in map_list:
        _key = (i.meta["energy_l"], i.meta["energy_h"])
        if _key not in dic:
            dic[_key] = [i]
        else:
            dic[_key].append(i)

    for i in dic:
        comp_map = sunpy.map.CompositeMap()
        n = len(dic[i])
        alpha = 1.0 / n
        for j in dic[i]:
            comp_map.add_map(j, zorder=0, alpha=alpha)
        dic[i] = comp_map

    return dic

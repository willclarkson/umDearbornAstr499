"""
Toolbox for the Photometry Pipeline
2016-03-09, mommermiscience@gmail.com
"""
from __future__ import print_function
from __future__ import division

# Photometry Pipeline
# Copyright (C) 2016-2018  Michael Mommert, mommermiscience@gmail.com

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see
# <http://www.gnu.org/licenses/>.

import sys
try:
    from past.utils import old_div
except ImportError:
    print('Module future not found. Please install with: pip install future')
    sys.exit()

import math
import numpy as np

# only import if Python3 is used
if sys.version_info > (3, 0):
    from future import standard_library
    standard_library.install_aliases()
    from builtins import range


# TIME AND DATE

def jd_to_gregorian(jd, is_mjd=False):
    """ convert a julian date into a gregorian data """
    if is_mjd:
        mjd = jd
    else:
        mjd = jd - 2400000.5

    MJD0 = 2400000.5  # 1858 November 17, 00:00:00 hours

    modf = math.modf
    a = int(mjd+MJD0+0.5)
    b = int(old_div((a-1867216.25), 36524.25))
    c = a + b - int(modf(old_div(b, 4))[1]) + 1525

    d = int(old_div((c-122.1), 365.25))
    e = 365*d + int(modf(old_div(d, 4))[1])
    f = int(old_div((c-e), 30.6001))

    day = int(c - e - int(30.6001*f))
    month = int(f - 1 - 12*int(modf(old_div(f, 14))[1]))
    year = int(d - 4715 - int(modf(old_div((7+month), 10))[1]))
    fracofday = mjd - math.floor(mjd)
    hour = int(math.floor(fracofday * 24.0))
    minute = int(math.floor(((fracofday*24.0)-hour)*60.))
    second = int(math.floor(((((fracofday*24.0)-hour)*60.)-minute)*60.))

    return (year, month, day, hour, minute, second)


def dateobs_to_jd(date):
    """convert a string of the format YYYY-MM-DDTHH:MM:SS into a julian
        date; 'T' is used as a separator between date and time
    """
    if 'T' in date:
        date = date.split('T')
    if ' ' in date:
        date = date.split(' ')
    time = date[1].split(':')
    date = date[0].split('-')

    # check if date is yyyy-mm-dd or dd-mm-yyyy
    if len(date[2]) == 4 and len(date[0]) < 3:
        date = date[::-1]

    a = (14 - float(date[1]))//12
    y = float(date[0]) + 4800 - a
    m = float(date[1]) + 12*a - 3
    return float(date[2]) + ((153*m + 2)//5) + 365*y + y//4 - y//100 \
        + y//400 - 32045.5 + old_div(float(time[0]), 24.) + old_div(float(time[1]), 1440.) \
        + old_div(float(time[2]), 86400.)


def jd_to_fractionalyear(jd, is_mjd=False):
    """ convert a julian date into a fractional year, e.g., 2000.123456 """
    if is_mjd:
        jd += 2400000.5
    date = jd_to_gregorian(jd)
    year = date[0]+old_div(date[1], 12.)+old_div(date[2], 365.) + \
        old_div(date[3], 8760.)+old_div(date[4], 525600.)
    return year


def fractionalyear_to_jd(date):
    """ convert a fractional year into a julian date """
    jd_jan1 = dateobs_to_jd('%4d-01-01T00:00:00' % math.floor(date))
    return jd_jan1 + 365*(date-math.floor(date))


# ASTROMATIC tools

def read_scamp_output():
    """ routine to read in the 'scamp.xml' file """
    raw = open('scamp_output.xml', 'r').readlines()
    headers, hdr_idx, data = {}, 0, []
    read_this, idx = False, 0
    while idx < len(raw):
        # read header
        if read_this and raw[idx].find('<FIELD name=') > -1:
            headers[raw[idx][raw[idx].find('<FIELD name')+13:
                             raw[idx].find('" datatype')]] = hdr_idx
            hdr_idx += 1
        # read data
        # new data line
        if read_this and raw[idx].find('<TR>') > -1:
            this_data = []
        # flush data line
        if read_this and raw[idx].find('</TR>') > -1:
            data.append(np.hstack(this_data))
        # actually read data line
        if read_this and raw[idx].find('<TD>') > -1:
            line = raw[idx].replace('</TD>', '<TD>').split('<TD>')
            for item in line:
                if len(item.strip()) > 0 and item.find('\n') == -1:
                    this_data.append(item)
        # control reading
        # activate reading
        if not read_this and \
           raw[idx].find('<TABLE ID="Fields" name="Fields">') > -1:
            read_this = True
        # deactivate reading
        if read_this and raw[idx].find('</TABLEDATA></DATA>') > -1:
            read_this = False
        idx += 1

    # check if data rows have same length as header
    for i in range(len(data)):
        if len(headers) != len(data[i]):
            raise (RuntimeError,
                   ('data and header lists from SCAMP output file have '
                    'different lengths for image %s; do the FITS files have the '
                    'OBJECT keyword populated?') % data[i][headers['Catalog_Name']])
    return (headers, data)


# PP tools

def get_binning(header, obsparam):
    """ derive binning from image header
        use obsparam['binning'] keywords, unless both keywords are set to 1
        return: tuple (binning_x, binning_y)"""

    if (isinstance(obsparam['binning'][0], int) and
            isinstance(obsparam['binning'][1], int)):
        binning_x = obsparam['binning'][0]
        binning_y = obsparam['binning'][1]
    elif '#' in obsparam['binning'][0]:
        if '#blank' in obsparam['binning'][0]:
            binning_x = float(header[obsparam['binning'][0].
                                     split('#')[0]].split()[0])
            binning_y = float(header[obsparam['binning'][1].
                                     split('#')[0]].split()[1])
        elif '#x' in obsparam['binning'][0]:
            binning_x = float(header[obsparam['binning'][0].
                                     split('#')[0]].split('x')[0])
            binning_y = float(header[obsparam['binning'][1].
                                     split('#')[0]].split('x')[1])
        elif '#_' in obsparam['binning'][0]:
            binning_x = float(header[obsparam['binning'][0].
                                     split('#')[0]].split('_')[0])
            binning_y = float(header[obsparam['binning'][1].
                                     split('#')[0]].split('_')[1])
        elif '#CH#' in obsparam['binning'][0]:
            # only for RATIR
            channel = header['INSTRUME'].strip()[1]
            binning_x = float(header[obsparam['binning'][0].
                                     replace('#CH#', channel)])
            binning_y = float(header[obsparam['binning'][1].
                                     replace('#CH#', channel)])
    else:
        binning_x = header[obsparam['binning'][0]]
        binning_y = header[obsparam['binning'][1]]

    return (binning_x, binning_y)


def skycenter(catalogs, ra_key='ra_deg', dec_key='dec_deg'):
    """derive center position and radius from catalogs"""
    from astropy.coordinates import SkyCoord
    from astropy import units as u

    # using percentiles instead of min/max to get better handle
    # on outliers
    min_ra = min([np.percentile(cat[ra_key], 1)
                  for cat in catalogs])
    max_ra = max([np.percentile(cat[ra_key], 99)
                  for cat in catalogs])
    min_dec = min([np.percentile(cat[dec_key], 1)
                   for cat in catalogs])
    max_dec = max([np.percentile(cat[dec_key], 99)
                   for cat in catalogs])

    ra, dec = (np.rad2deg(np.angle(np.exp(1j*np.deg2rad(min_ra)) +
                                   np.exp(1j*np.deg2rad(max_ra)))),
               np.rad2deg(np.angle(np.exp(1j*np.deg2rad(min_dec)) +
                                   np.exp(1j*np.deg2rad(max_dec)))))

    lower_left = SkyCoord(ra=min_ra, dec=min_dec, frame='icrs', unit='deg')
    upper_right = SkyCoord(ra=max_ra, dec=max_dec, frame='icrs', unit='deg')

    rad = lower_left.separation(upper_right).deg/2

    return ra, dec, rad


# miscellaneous tools

def if_val_in_dict(target_val, dic):
    """check if a value appears in a nested dict structure"""
    result = False
    for key, val in dic.items():
        if type(val) is dict:
            if if_val_in_dict(target_val, val):
                result = True
        elif type(val) is list:
            if target_val in val:
                result = True
        else:
            if target_val == val:
                result = True
    return result

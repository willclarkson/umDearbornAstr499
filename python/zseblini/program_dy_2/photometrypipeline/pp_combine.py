#!/usr/bin/env python3

""" PP_COMBINE - combine frames based on wcs
    v1.0: 2017-10-03, mommermiscience@gmail.com
"""
from __future__ import print_function, division

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


import numpy
import os
import sys
import shutil
import logging
import subprocess
import argparse
import shlex
import time
from astropy.io import fits
from past.utils import old_div
from astroquery.jplhorizons import Horizons

# pipeline-specific modules
import _pp_conf
import toolbox

# create a portable DEVNULL
# necessary to prevent subprocess.PIPE and STDOUT from clogging if
# Source Extractor runs for too long
try:
    from subprocess import DEVNULL  # Py3
except ImportError:
    import os  # Py2
    DEVNULL = open(os.devnull, 'wb')

# only import if Python3 is used
if sys.version_info > (3, 0):
    from builtins import str

# setup logging
logging.basicConfig(filename=_pp_conf.log_filename,
                    level=_pp_conf.log_level,
                    format=_pp_conf.log_formatline,
                    datefmt=_pp_conf.log_datefmt)


def combine(filenames, obsparam, comoving, targetname,
            manual_rates, combine_method, keep_files,
            backsub=False, display=True, diagnostics=True):
    """
    image combination wrapper
    output: diagnostic properties
    """

    # start logging
    logging.info('starting image combination with parameters: %s' %
                 (', '.join([('%s: %s' % (var, str(val))) for
                             var, val in list(locals().items())])))

    # check if images have been run through pp_prepare
    try:
        midtime_jd = fits.open(filenames[0], verify='silentfix',
                               ignore_missing_end=True)[0].header['MIDTIMJD']
    except KeyError:
        raise KeyError(('%s image header incomplete, have the data run ' +
                        'through pp_prepare?') % filenames[0])
        return None

    # adopt first frame as reference frame
    hdulist = fits.open(filenames[0])
    header = hdulist[0].header
    refdate = float(header['MIDTIMJD'])
    # read out ra and dec from header
    if obsparam['radec_separator'] == 'XXX':
        ref_ra_deg = float(header[obsparam['ra']])
        ref_dec_deg = float(header[obsparam['dec']])
        if obsparam['telescope_keyword'] == 'UKIRTWFCAM':
            ref_ra_deg = ref_ra_deg/24.*360. - 795/3600.
            ref_dec_deg -= 795/3600.
    else:
        ra_string = header[obsparam['ra']].split(
            obsparam['radec_separator'])
        dec_string = header[obsparam['dec']].split(
            obsparam['radec_separator'])
        ref_ra_deg = 15.*(float(ra_string[0]) +
                          old_div(float(ra_string[1]), 60.) +
                          old_div(float(ra_string[2]), 3600.))
        ref_dec_deg = (abs(float(dec_string[0])) +
                       old_div(float(dec_string[1]), 60.) +
                       old_div(float(dec_string[2]), 3600.))
        if dec_string[0].find('-') > -1:
            ref_dec_deg = -1 * ref_dec_deg

    if obsparam['telescope_keyword'] == 'UKIRTWFCAM':
        ref_ra_deg = ref_ra_deg/24.*360.

    if obsparam['telescope_keyword'] == "UKIRTWFCAM":
        ref_ra_deg -= float(header['TRAOFF'])/3600
        ref_dec_deg -= float(header['TDECOFF'])/3600

    hdulist.close()

    # modify individual frames if comoving == True
    if comoving:
        movingfilenames = []

        # sort filenames by MIDTIMJD
        mjds = []
        for filename in filenames:
            hdulist = fits.open(filename)
            mjds.append(float(hdulist[0].header['MIDTIMJD']))
        filenames = [filenames[i] for i in numpy.argsort(mjds)]

        for filename in filenames:
            movingfilename = filename[:filename.find('.fits')]+'_moving.fits'
            print('shifting %s -> %s' % (filename, movingfilename))
            logging.info('shifting %s -> %s' % (filename, movingfilename))

            # read out date and pointing information
            hdulist = fits.open(filename)
            header = hdulist[0].header
            date = hdulist[0].header['MIDTIMJD']
            data = hdulist[0].data
            hdulist.close()

            # use ephemerides from Horizons if no manual rates are provided
            if manual_rates is None:
                # call HORIZONS to get target coordinates
                obj = Horizons(targetname.replace('_', ' '), epochs=date,
                               location=str(obsparam['observatory_code']))
                try:
                    eph = obj.ephemerides()
                    n = len(eph)
                except ValueError:
                    print('Target (%s) not an asteroid' % targetname)
                    logging.warning('Target (%s) not an asteroid' % targetname)
                    n = None
                    time.sleep(0.5)
                if n is None or n == 0:
                    logging.warning('WARNING: No position from Horizons!' +
                                    'Name (%s) correct?' % targetname)
                    logging.warning('HORIZONS call: %s' % eph.url)
                    raise(ValueError, 'no Horizons ephemerides available')
                else:
                    logging.info('ephemerides for %s pulled from Horizons' %
                                 targetname)
                    logging.info('Horizons call: %s' %
                                 obj.uri)

                    target_ra, target_dec = eph[0]['RA'], eph[0]['DEC']

                # get image pointing from header
                if obsparam['radec_separator'] == 'XXX':
                    ra_deg = float(header[obsparam['ra']])
                    dec_deg = float(header[obsparam['dec']])
                    if obsparam['telescope_keyword'] == 'UKIRTWFCAM':
                        ra_deg = ra_deg/24.*360. - 795/3600.
                        dec_deg -= 795/3600.
                else:
                    ra_string = header[obsparam['ra']].split(
                        obsparam['radec_separator'])
                    dec_string = header[obsparam['dec']].split(
                        obsparam['radec_separator'])
                    ra_deg = 15.*(float(ra_string[0]) +
                                  old_div(float(ra_string[1]), 60.) +
                                  old_div(float(ra_string[2]), 3600.))
                    dec_deg = (abs(float(dec_string[0])) +
                               old_div(float(dec_string[1]), 60.) +
                               old_div(float(dec_string[2]), 3600.))
                    if dec_string[0].find('-') > -1:
                        dec_deg = -1 * dec_deg

                if filename == filenames[0]:
                    ref_offset_ra = target_ra - ref_ra_deg
                    ref_offset_dec = target_dec - ref_dec_deg

                offset_ra = target_ra - ref_ra_deg - ref_offset_ra
                offset_dec = target_dec - ref_dec_deg - ref_offset_dec

            else:
                # use manual rates (since they are provided)
                offset_ra = ((float(header['MIDTIMJD'])-refdate)*86400 *
                             float(manual_rates[0]))/3600
                offset_dec = ((float(header['MIDTIMJD'])-refdate)*86400 *
                              float(manual_rates[1]))/3600

            logging.info('offsets in RA and Dec: %f, %f arcsec' %
                         (offset_ra*3600, offset_dec*3600))

            crval1 = float(header['CRVAL1'])
            crval2 = float(header['CRVAL2'])

            # write new CRVALi keywords in different file
            new_hdu = fits.PrimaryHDU(data)
            new_hdu.header = header
            new_hdu.header['CRVAL1'] = (crval1-offset_ra,
                                        'updated in the moving frame of the object')
            new_hdu.header['CRVAL2'] = (crval2-offset_dec,
                                        'updated in the moving frame of the object')
            movingfilenames.append(movingfilename)
            new_hdu.writeto(movingfilename, overwrite=True,
                            output_verify='silentfix')

    if comoving:
        outfile_name = 'comove.fits'
        fileline = " ".join(movingfilenames)
        n_frames = len(movingfilenames)
    else:
        outfile_name = 'skycoadd.fits'
        fileline = " ".join(filenames)
        n_frames = len(filenames)

    # run swarp on all image catalogs using different catalogs
    commandline = (('swarp -combine Y -combine_type %s -delete_tmpfiles ' +
                    'Y -imageout_name %s -interpolate Y -subtract_back %s ' +
                    '-weight_type NONE -copy_keywords %s -write_xml N ' +
                    '-CENTER_TYPE MOST %s') %
                   ({'median': 'MEDIAN', 'average': 'AVERAGE',
                     'clipped': 'CLIPPED -CLIP_AMPFRAC 0.2 -CLIP_SIGMA 0.1 '}
                    [combine_method],
                    outfile_name,
                    {True: 'Y', False: 'N'}[backsub],
                    obsparam['copy_keywords'], fileline))

    logging.info('call SWARP as: %s' % commandline)
    print('running SWARP to combine {:d} frames...'.format(n_frames))

    try:
        swarp = subprocess.Popen(shlex.split(commandline),
                                 stdout=DEVNULL,
                                 stderr=DEVNULL,
                                 close_fds=True)
        # do not direct stdout to subprocess.PIPE:
        # for large FITS files, PIPE will clog, stalling
        # subprocess.Popen
    except Exception as e:
        print('SWARP call:', (e))
        logging.error('SWARP call:', (e))
        return None

    swarp.wait()
    print('done!')

    # remove files that are not needed anymore
    if not keep_files:
        if comoving:
            for filename in movingfilenames:
                os.remove(filename)

    # update combined image header
    total_exptime = 0
    for filename in filenames:
        hdulist = fits.open(filename)
        total_exptime += float(hdulist[0].header[obsparam['exptime']])

    hdulist = fits.open(outfile_name, mode='update')
    hdulist[0].header[obsparam['exptime']] = (total_exptime, 'PP: cumulative')
    hdulist[0].header['COMBO_N'] = (len(filenames), 'PP: N files combo')
    hdulist[0].header['COMBO_M'] = (combine_method, 'PP: combo method')
    hdulist[0].header['COMOVE'] = (str(comoving), 'PP: comoving?')
    hdulist.flush()

    return n_frames


if __name__ == '__main__':

    # command line arguments
    parser = argparse.ArgumentParser(description='image combination')
    parser.add_argument("-comoving", action="store_true",
                        help='combine in moving target frame')
    parser.add_argument("-targetname",
                        help='moving target name')
    parser.add_argument("-manual_rates", help='manual rates in arcsec/s',
                        nargs=2)
    parser.add_argument('-method',
                        help='combination method',
                        choices=['average', 'median', 'clipped'],
                        default='clipped')
    parser.add_argument("-backsub", action="store_true",
                        help='subtract background in each frame ')
    parser.add_argument("-keep_files", action="store_true",
                        help='keep intermediate files', default=False)
    parser.add_argument('images', help='images to process', nargs='+')

    args = parser.parse_args()
    comoving = args.comoving
    targetname = args.targetname
    manual_rates = args.manual_rates
    combine_method = args.method
    backsub = args.backsub
    keep_files = args.keep_files
    filenames = args.images

    # read telescope and filter information from fits headers
    # check that they are the same for all images
    instruments = []
    for filename in filenames:
        hdulist = fits.open(filename, ignore_missing_end=True,
                            verify='silentfix')
        header = hdulist[0].header
        for key in _pp_conf.instrument_keys:
            if key in header:
                instruments.append(header[key])

    if len(instruments) == 0:
        raise KeyError('cannot identify telescope/instrument; please update'
                       '_pp_conf.instrument_keys accordingly')

    # assign telescope parameters (telescopes.py)
    telescope = _pp_conf.instrument_identifiers[instruments[0]]
    obsparam = _pp_conf.telescope_parameters[telescope]

    if manual_rates is not None:
        comoving = True

    if comoving and targetname is None:
        targetname = header[obsparam['object']]

    # run image combination wrapper
    combination = combine(filenames, obsparam, comoving, targetname,
                          manual_rates, combine_method, keep_files,
                          backsub, display=True, diagnostics=True)

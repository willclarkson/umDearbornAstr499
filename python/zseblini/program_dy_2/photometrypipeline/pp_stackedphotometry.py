#!/usr/bin/env python3

""" PP_STACKEDPHOTOMETRY - wrapper to perform photometry on stacked images
    v1.0: 2017-10-19, mommermiscience@gmail.com
"""
from __future__ import print_function

import diagnostics as diag
import pp_combine
import pp_distill
import pp_calibrate
import pp_photometry
import pp_register
import pp_extract
import pp_prepare
from catalog import *
import _pp_conf


# Photometry Pipeline
# Copyright (C) 2016-2018 Michael Mommert, mommermiscience@gmail.com

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

import re
import os
import gc
import sys
try:
    import numpy as np
except ImportError:
    print('Module numpy not found. Please install with: pip install numpy')
    sys.exit()
import shutil
import logging
import subprocess
import argparse
import shlex
import time
try:
    from astropy.io import fits
except ImportError:
    print('Module astropy not found. Please install with: pip install astropy')
    sys.exit()

# only import if Python3 is used
if sys.version_info > (3, 0):
    from builtins import str
    from builtins import range

# pipeline-specific modules

# setup logging
logging.basicConfig(filename=_pp_conf.log_filename,
                    level=_pp_conf.log_level,
                    format=_pp_conf.log_formatline,
                    datefmt=_pp_conf.log_datefmt)

if __name__ == '__main__':

    # command line arguments
    parser = argparse.ArgumentParser(description='stacked photometry')
    parser.add_argument('-comoving', help='stack in moving target frame',
                        action='store_true')
    parser.add_argument("-targetname",
                        help='moving target name', default=None)
    parser.add_argument('-filter', help='filter name override',
                        default=None)
    parser.add_argument('-method',
                        help='combination method',
                        choices=['average', 'median', 'clipped'],
                        default='clipped')
    parser.add_argument('-fixed_aprad', help='fixed aperture radius (px)',
                        default=0)
    parser.add_argument('-snr',
                        help='SNR limit for detected sources',
                        default=1.5)
    parser.add_argument('-solar',
                        help='restrict to solar-color stars',
                        action="store_true", default=False)
    parser.add_argument('-reject',
                        help='schemas for target rejection',
                        nargs=1, default='pos')

    parser.add_argument('images', help='images to process',
                        nargs='+')
    args = parser.parse_args()
    comoving = args.comoving
    targetname = args.targetname
    man_filtername = args.filter
    combinemethod = args.method
    fixed_aprad = float(args.fixed_aprad)
    snr = float(args.snr)
    solar = args.solar
    #reject = args.reject[0]
    reject = args.reject
    filenames = args.images

    # use current directory as root directory
    rootdir = os.getcwd()

    # check if input filenames is actually a list
    if len(filenames) == 1:
        if filenames[0].find('.lst') > -1 or filenames[0].find('.list') > -1:
            filenames = [filename[:-1] for filename in open(filenames[0], 'r')
                         .readlines()]

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

    # ------------------- SKYCOADD

    # create skycoadd in current directory
    ppcombine_comoving = False
    #targetname = None
    manual_rates = None
    keep_files = False
    combination = pp_combine.combine(filenames, obsparam, ppcombine_comoving,
                                     targetname, manual_rates,
                                     combinemethod, keep_files,
                                     display=True, diagnostics=True)

    # create separate directory to analyze skycoadd data
    if os.path.exists('skycoadd/'):
        shutil.rmtree('skycoadd/')
    os.mkdir('skycoadd/')
    os.rename('skycoadd.fits', 'skycoadd/skycoadd.fits')
    os.chdir('skycoadd/')

    # # diagnostics and logging for skycoadd and comove go into respective dirs
    # _pp_conf.dataroot, _pp_conf.diagroot, \
    #     _pp_conf.index_filename, _pp_conf.reg_filename, _pp_conf.cal_filename, \
    #     _pp_conf.res_filename = _pp_conf.setup_diagnostics()

    # setup logging again
    logging.basicConfig(filename=_pp_conf.log_filename,
                        level=_pp_conf.log_level,
                        format=_pp_conf.log_formatline,
                        datefmt=_pp_conf.log_datefmt)

    logging.info('create skycoadd.fits from images: %s' % ','.join(filenames))
    logging.info('move skycoadd.fits into skycoadd/ directory')

    # prepare image
    preparation = pp_prepare.prepare(['skycoadd.fits'], obsparam,
                                     {}, keep_wcs=True,
                                     diagnostics=True, display=True)

    # run photometry (curve-of-growth analysis)
    source_minarea = obsparam['source_minarea']
    background_only = True
    target_only = False
    if fixed_aprad == 0:
        aprad = None  # force curve-of-growth analysis
        print('\n----- derive optimum photometry aperture\n')
        logging.info('----- derive optimum photometry aperture')
    else:
        aprad = fixed_aprad  # skip curve_of_growth analysis
        print('\n----- use fixed aperture radius (%5.2f)\n' % fixed_aprad)
        logging.info('----- use fixed aperture radius (%5.2f)' % fixed_aprad)

    phot = pp_photometry.photometry(['skycoadd.fits'], snr, source_minarea,
                                    aprad,
                                    None, background_only,
                                    target_only,
                                    telescope, obsparam, display=True,
                                    diagnostics=True)

    # data went through curve-of-growth analysis
    if phot is not None:
        aprad = phot['optimum_aprad']
    # a fixed aperture radius has been used
    else:
        aprad = fixed_aprad

    # run photometric calibration
    minstars = _pp_conf.minstars
    manualcatalog = None
    if man_filtername is None:
        man_filtername = False

    print('\n----- run photometric calibration\n')

    calibration = pp_calibrate.calibrate(['skycoadd.fits'], minstars,
                                         man_filtername,
                                         manualcatalog, obsparam, solar=solar,
                                         display=True,
                                         diagnostics=True)

    zp = calibration['zeropoints'][0]['zp']
    zp_err = calibration['zeropoints'][0]['zp_sig']

    logging.info('zeropoint derived from skycoadd.fits: %5.2f+-%4.2f' %
                 (zp, zp_err))

    os.chdir(rootdir)

    # ------------------- COMOVE

    hdulist = fits.open(filenames[0])
    if comoving and targetname is None:
        targetname = header[obsparam['object']]

    if comoving:

        # create comove in current directory
        ppcombine_comoving = True
        manual_rates = None
        keep_files = False
        combination = pp_combine.combine(filenames, obsparam,
                                         ppcombine_comoving,
                                         targetname, manual_rates,
                                         combinemethod, keep_files,
                                         display=True, diagnostics=True)

        # create separate directory to analyze skycoadd data
        if os.path.exists('comove/'):
            shutil.rmtree('comove/')
        os.mkdir('comove/')
        os.rename('comove.fits', 'comove/comove.fits')
        os.chdir('comove/')

        # diagnostics + logging for skycoadd and comove go into respective dirs
        # _pp_conf.dataroot, _pp_conf.diagroot, \
        #     _pp_conf.index_filename, _pp_conf.reg_filename, \
        #     _pp_conf.cal_filename, \
        #     _pp_conf.res_filename = _pp_conf.setup_diagnostics()

        # setup logging again
        logging.basicConfig(filename=_pp_conf.log_filename,
                            level=_pp_conf.log_level,
                            format=_pp_conf.log_formatline,
                            datefmt=_pp_conf.log_datefmt)

        logging.info('create comove.fits from images: %s' %
                     ','.join(filenames))
        logging.info('move comove.fits into comove/ directory')

        # prepare image
        preparation = pp_prepare.prepare(['comove.fits'], obsparam,
                                         {}, keep_wcs=True,
                                         diagnostics=True, display=True)

        # run photometry (curve-of-growth analysis)
        source_minarea = obsparam['source_minarea']
        background_only = False
        target_only = False

        print('\n----- use skycoadd optimum photometry aperture (%4.2f)\n' %
              aprad)
        phot = pp_photometry.photometry(['comove.fits'], snr, source_minarea,
                                        aprad,
                                        None, background_only,
                                        target_only,
                                        telescope, obsparam, display=True,
                                        diagnostics=True)

        # run photometric calibration (instrumental)
        minstars = _pp_conf.minstars
        man_filtername = obsparam['filter_translations'][hdulist[0].header['filter']]
        manualcatalog = None

        print('\n----- run photometric calibration\n')
        logging.info('use skycoadd.fits magnitude zeropoint: %5.2f+-%4.2f' %
                     (zp, zp_err))
        calibration = pp_calibrate.calibrate(['comove.fits'], minstars,
                                             man_filtername,
                                             manualcatalog, obsparam,
                                             magzp=(zp, zp_err), solar=solar,
                                             display=True,
                                             diagnostics=True)

    # distill target brightness from database
    #man_targetname = None
    if targetname != header[obsparam['object']]:
        man_targetname = targetname
    man_offset = [0, 0]
    fixed_targets_file = None
    posfile = None
    distillate = pp_distill.distill(calibration['catalogs'],
                                    man_targetname, man_offset,
                                    fixed_targets_file, posfile,
                                    rejectionfilter=reject,
                                    display=True, diagnostics=True)

    os.chdir(rootdir)

    logging.info('move comove/photometry*.dat to root directory')
    targets = np.array(list(distillate['targetnames'].keys()))
    target = targets[targets != 'control_star'][0]
    shutil.copyfile(('comove/photometry_%s.dat' %
                     target.translate(_pp_conf.target2filename)),
                    ('photometry_%s.dat' %
                     target.translate(_pp_conf.target2filename)))
    logging.info('move skycoadd.fits into skycoadd/ directory')

    print('\nDone!\n')
    logging.info('----- successfully done with this process ----')

    gc.collect()  # collect garbage; just in case, you never know...

#!/usr/bin/env python3

""" PPTOOL_PSFSUB - PSF subtraction tool
    v1.0: 2017-12-10, mommermiscience@gmail.com
"""
from __future__ import print_function
from __future__ import division

# Photometry Pipeline
# Copyright (C) 2016-2018  ichael Mommert, mommermiscience@gmail.com

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


from past.utils import old_div
import numpy
import os
import sys
import subprocess
import logging
import argparse
import time
from copy import deepcopy
from astropy.io import fits
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt
import callhorizons

# only import if Python3 is used
if sys.version_info > (3, 0):
    from builtins import str
    from builtins import range

# pipeline-specific modules
import _pp_conf
import pp_extract
from catalog import *
from toolbox import *
import diagnostics as diag

# setup logging
logging.basicConfig(filename=_pp_conf.log_filename,
                    level=_pp_conf.log_level,
                    format=_pp_conf.log_formatline,
                    datefmt=_pp_conf.log_datefmt)


def psfsubtraction(filenames, display=False, diagnostics=False):

    pass

# MAIN


if __name__ == '__main__':

    # define command line arguments
    parser = argparse.ArgumentParser(description='automated PSF subtraction')
    # parser.add_argument('-snr', help='sextractor SNR threshold for '+\
    #                     'photometry catalog', default=2)
    # parser.add_argument('-minarea', help='sextractor SNR threshold for '+\
    #                     'photometry catalog', default=0)
    # parser.add_argument('-aprad', help='aperture radius for photometry (px)',
    #                     default=None)
    # parser.add_argument('-target',
    #                     help='object name override (e.g., 2015_AB123)',
    #                     default=None)
    # parser.add_argument('-background_only',
    #                     help='find aperture for background only',
    #                     action="store_true")
    # parser.add_argument('-target_only', help='find aperture for target only',
    #                     action="store_true")
    parser.add_argument('images', help='images to process', nargs='+')

    args = parser.parse_args()
    # sex_snr = float(args.snr)
    # source_minarea = float(args.minarea)
    # aprad = float(args.aprad) if args.aprad is not None else None
    # manobjectname = args.target
    # background_only = args.background_only
    # target_only = args.target_only
    filenames = args.images

    # check if input filenames is actually a list
    if len(filenames) == 1:
        if filenames[0].find('.lst') > -1 or filenames[0].find('.list') > -1:
            filenames = [filename[:-1] for filename in open(filenames[0], 'r')
                         .readlines()]

    # obtain telescope information
    hdu = fits.open(filenames[0], ignore_missing_end=True)
    try:
        telescope = hdu[0].header['TEL_KEYW']
    except KeyError:
        print('ERROR: cannot find telescope keyword in image header;' +
              'has this image run through wcs_register?')
        sys.exit(0)
    obsparam = _pp_conf.telescope_parameters[telescope]

    if type(manobjectname) == str:
        manobjectname = manobjectname.translate(_pp_conf.target2filename)

    phot = psfsubtraction(filenames, display=True,
                          diagnostics=True)

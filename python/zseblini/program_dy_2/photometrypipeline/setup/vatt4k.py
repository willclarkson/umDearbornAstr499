"""
Photometry Pipeline Configuation File for VATT/VATT4k
2017-02-05, mommermiscience@gmail.com
"""

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


vatt4k_param = {'telescope_instrument': 'VATT/VATT4k',
                telescope/instrument name 'telescope_keyword': 'VATT4K',
                telescope/instrument keyword 'observatory_code': '290',  # MPC
                observatory code 'secpix': (0.1875, 0.1875),  # pixel size
                (arcsec)  # before binning # image orientation preferences 'flipx'
                : True, 'flipy': False, 'rotate': 0,

                # instrument-specific FITS header keywords
                'binning': ('CCDBIN1', 'CCDBIN2'),  # binning in x/y
                'extent': ('NAXIS1', 'NAXIS2'),   # N_pixels in x/y
                'ra': 'RA',  # telescope pointing, RA
                'dec': 'DEC',  # telescope pointin, Dec
                'radec_separator': ':',   # RA/Dec hms separator, use 'XXX'
                # if already in degrees
                'date_keyword': 'DATE-OBS|TIME-OBS',  # obs date/time
                # keyword; use
                # 'date|time' if
                # separate
                'object': 'OBJECT',  # object name keyword
                'filter': 'FILTER',  # filter keyword
                'filter_translations': {'TOP 2 BOT 1': 'V', 'TOP 3 BOT 1': 'R',
                                        'TOP 4 BOT 1': 'I', 'TOP 5 BOT 1': 'B'},
                # filtername translation dictionary
                'exptime': 'EXPTIME',  # exposure time keyword (s)
                'airmass': 'AIRMASS',  # airmass keyword

                # source extractor settings
                'sex': {'DETECT_MINAREA':  12,
                        'DETECT_THRESH':   3,
                        'ANALYSIS_THRESH': 3,
                        'CATALOG_NAME':    'VATT4K.ldac',
                        'aprad_default':   5,
                        # [minimum, maximum] aperture rad (px)
                        'aprad_range':     [2, 10],
                        'PHOT_APERTURES':  5,
                        'WEIGHT_TYPE':     'NONE',
                        'WEIGHT_IMAGE':    'NONE',
                        'mask_files': {},  # as a function of x,y binning
                        'BACKPHOTO_TYPE':  'GLOBAL',
                        'PARAMETERS_NAME': '$PHOTPIPEDIR/setup/singleaperture.sexparam',
                        'SATUR_LEVEL:':    50000,
                        'SATUR_KEY':       'NONE'},

                # scamp settings
                'scamp-config-file': rootpath+'/setup/vatt4k.scamp',
                'reg_max_mag': 19,
                'reg_search_radius': 0.5,  # deg

                # swarp settings
                'copy_keywords': ('TELESCOP,INSTRUME,FILTER,EXPTIME,OBJECT,' +
                                  'DATE-OBS,TIME-OBS,RA,DEC,SECPIX,AIRMASS,' +
                                  'TEL_KEYW'),
                #                         keywords to be copied in image
                #                         combination using swarp
                'swarp-config-file': rootpath+'/setup/vatt4k.swarp',

                # default catalog settings
                'astrometry_catalogs': ['GAIA'],
                'photometry_catalogs': ['SDSS-R9', 'APASS9']
                }

# add telescope parameters to according lists and dictionaries
implemented_telescopes.append('VATT4K')
instrument_identifiers['= "Vatt4k"'] = 'VATT4K'
telescope_parameters['VATT4K'] = vatt4k_param

"""
Photometry Pipeline Configuation File
2016-03-09, mommermiscience@gmail.com
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

# telescope/instrument configurations

# VATT, VATT4k
vatt4k_param = {
    'telescope_instrument': 'VATT/VATT4k',  # telescope/instrument name
    'telescope_keyword': 'VATT4K',      # telescope/instrument keyword
    'observatory_code': '290',         # MPC observatory code
    'secpix': (0.1875, 0.1875),  # pixel size (arcsec)
    # before binning
    'ext_coeff': 0.05,          # typical extinction coefficient


    # image orientation preferences
    'flipx': True,
    'flipy': False,
    'rotate': 0,

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
    'obsmidtime_jd': 'MIDTIMJD',  # obs midtime jd keyword
    # (usually provided by
    # pp_prepare
    'object': 'OBJECT',  # object name keyword
    'filter': 'FILTER',  # filter keyword
    'filter_translations': {'TOP 2 BOT 1': 'V', 'TOP 3 BOT 1': 'R',
                            'TOP 4 BOT 1': 'I', 'TOP 5 BOT 1': 'B',
                            'TOP 1 BOT 1': None},
    # filtername translation dictionary
    'exptime': 'EXPTIME',  # exposure time keyword (s)
    'airmass': 'AIRMASS',  # airmass keyword


    # source extractor settings
    'source_minarea': 12,  # default sextractor source minimum N_pixels
    'source_snr': 3,  # default sextractor source snr for registration
    'aprad_default': 5,  # default aperture radius in px
    'aprad_range': [2, 10],  # [minimum, maximum] aperture radius (px)
    'sex-config-file': rootpath+'/setup/vatt4k.sex',
    'mask_file': {},
    #                        mask files as a function of x,y binning

    # registration settings (Scamp)
    'scamp-config-file': rootpath+'/setup/vatt4k.scamp',
    'reg_max_mag': 19,
    'reg_search_radius': 0.5,  # deg
    'source_tolerance': 'high',

    # swarp settings
    'copy_keywords': ('TELESCOP,INSTRUME,FILTER,EXPTIME,OBJECT,' +
                      'DATE-OBS,TIME-OBS,RA,DEC,SECPIX,AIRMASS,' +
                      'TEL_KEYW,CCDBIN1,CCDBIN2,MIDTIMJD'),
    #                         keywords to be copied in image
    #                         combination using swarp
    'swarp-config-file': rootpath+'/setup/vatt4k.swarp',

    # default catalog settings
    'astrometry_catalogs': ['GAIA'],
    'photometry_catalogs': ['PANSTARRS', 'SDSS-R9', 'APASS9']
}


# DCT, LMI
dctlmi_param = {
    'telescope_instrument': 'DCT/LMI',  # telescope/instrument name
    'telescope_keyword': 'DCTLMI',  # telescope/instrument keyword
    'observatory_code': 'G37',         # MPC observatory code
    'secpix': (0.12, 0.12),  # pixel size (arcsec)
    # before binning
    'ext_coeff': 0.05,          # typical extinction coefficient


    # image orientation preferences
    'flipx': True,
    'flipy': False,
    'rotate': 0,

    # instrument-specific FITS header keywords
    'binning': ('CCDSUM#blank0', 'CCDSUM#blank1'),
    # binning in x/y, '_blankN' denotes that both axes
    # are listed in one keyword, sep. by blanks
    'extent': ('NAXIS1', 'NAXIS2'),   # N_pixels in x/y
    'ra': 'RA',  # telescope pointing, RA
    'dec': 'DEC',  # telescope pointin, Dec
    'radec_separator': ':',   # RA/Dec hms separator, use 'XXX'
    # if already in degrees
    'date_keyword': 'DATE-OBS',  # obs date/time
    # keyword; use
    # 'date|time' if
    # separate
    'obsmidtime_jd': 'MIDTIMJD',  # obs midtime jd keyword
    # (usually provided by
    # pp_prepare
    'object': 'OBJECT',  # object name keyword
    'filter': 'FILTERS',  # filter keyword
    'filter_translations': {'V': 'V', 'R': 'R', 'B': 'B', 'VR': None,
                              'I': 'I', 'SDSS-U': 'u', 'SDSS-G': 'g',
                              'SDSS-R': 'r', 'SDSS-I': 'i',
                              'SDSS-Z': 'z', 'OH': None, 'CN': None,
                              'UC': None, 'NH': None, 'BC': None,
                              'C2': None, 'C3': None, 'CO+': None,
                              'H2O+': None, 'GC': None, 'RC': None},
    # filtername translation dictionary
    'exptime': 'EXPTIME',  # exposure time keyword (s)
    'airmass': 'AIRMASS',  # airmass keyword


    # source extractor settings
    'source_minarea': 9,  # default sextractor source minimum N_pixels
    'source_snr': 3,  # default sextractor source snr for registration
    'aprad_default': 4,  # default aperture radius in px
    'aprad_range': [2, 10],  # [minimum, maximum] aperture radius (px)
    'sex-config-file': rootpath+'/setup/dctlmi.sex',
    'mask_file': {},
    #                        mask files as a function of x,y binning

    # registration settings (Scamp)
    'scamp-config-file': rootpath+'/setup/dctlmi.scamp',
    'reg_max_mag': 19,
    'reg_search_radius': 0.5,  # deg
    'source_tolerance': 'high',

    # swarp settings
    'copy_keywords': ('OBSERVAT,INSTRUME,EXPTIME,OBJECT,' +
                      'DATE-OBS,RA,DEC,AIRMASS,TEL_KEYW,CCDSUM,' +
                      'FILTERS,MIDTIMJD'),
    #                        keywords to be copied in image
    #                        combination using swarp
    'swarp-config-file': rootpath+'/setup/dctlmi.swarp',

    # default catalog settings
    'astrometry_catalogs': ['GAIA'],
    'photometry_catalogs': ['PANSTARRS', 'SDSS-R9', 'APASS9']
}


# Apache Point ARC 3.5m, ARCTIC
arc35arctic_param = {
    'telescope_instrument': 'ARC35m/ARCTIC',  # telescope/instrument name
    'telescope_keyword': 'ARC35ARCTIC',   # telescope/instrument keyword
    'observatory_code': '705',         # MPC observatory code
    'secpix': (0.115, 0.115),  # pixel size (arcsec)
    # before binning
    'ext_coeff': 0.05,          # typical extinction coefficient


    # image orientation preferences
    'flipx': True,
    'flipy': False,
    'rotate': 0,

    # instrument-specific FITS header keywords
    'binning': ('CCDBIN1', 'CCDBIN2'),
    # binning in x/y, '_blankN' denotes that both axes
    # are listed in one keyword, sep. by blanks
    'extent': ('NAXIS1', 'NAXIS2'),   # N_pixels in x/y
    'ra': 'RA',  # telescope pointing, RA
    'dec': 'DEC',  # telescope pointin, Dec
    'radec_separator': ':',   # RA/Dec hms separator, use 'XXX'
    # if already in degrees
    'date_keyword': 'DATE-OBS',  # obs date/time
    # keyword; use
    # 'date|time' if
    # separate
    'obsmidtime_jd': 'MIDTIMJD',  # obs midtime jd keyword
    # (usually provided by
    # pp_prepare
    'object': 'OBJNAME',  # object name keyword
    'filter': 'FILTER',  # filter keyword
    'filter_translations': {'SDSS U': 'u', 'SDSS G': 'g', 'SDSS R': 'r',
                              'SDSS I': 'i', 'SDSS Z': 'z', 'clear': None},
    # filtername translation dictionary
    'exptime': 'EXPTIME',  # exposure time keyword (s)
    'airmass': 'AIRMASS',  # airmass keyword


    # source extractor settings
    'source_minarea': 12,  # default sextractor source minimum N_pixels
    'source_snr': 3,  # default sextractor source snr for registration
    'aprad_default': 4,  # default aperture radius in px
    'aprad_range': [2, 10],  # [minimum, maximum] aperture radius (px)
    'sex-config-file': rootpath+'/setup/arc35arctic.sex',
    'mask_file': {},
    #                        mask files as a function of x,y binning

    # registration settings (Scamp)
    'scamp-config-file': rootpath+'/setup/arc35arctic.scamp',
    'reg_max_mag': 19,
    'reg_search_radius': 0.5,  # deg
    'source_tolerance': 'high',

    # swarp settings
    'copy_keywords': ('OBSERVAT,INSTRUME,FILTER,EXPTIME,OBJNAME,' +
                      'DATE-OBS,RA,DEC,AIRMASS,SECPIX,TEL_KEYW'),
    #                       keywords to be copied in image
    #                        combination using swarp
    'swarp-config-file': rootpath+'/setup/arc35arctic.swarp',

    # default catalog settings
    'astrometry_catalogs': ['GAIA'],
    'photometry_catalogs': ['SDSS-R9', 'APASS9', 'PANSTARRS', '2MASS']
}


# Apache Point ARC 3.5m, AGILE
arc35agile_param = {
    'telescope_instrument': 'ARC35m/AGILE',  # telescope/instrument name
    'telescope_keyword': 'ARC35AGILE',   # telescope/instrument keyword
    'observatory_code': '705',         # MPC observatory code
    'secpix': (0.13, 0.13),  # pixel size (arcsec)
    # before binning
    'ext_coeff': 0.05,          # typical extinction coefficient


    # image orientation preferences
    'flipx': True,
    'flipy': False,
    'rotate': 0,

    # instrument-specific FITS header keywords
    'binning': ('BINX', 'BINY'),
    # binning in x/y, '_blankN' denotes that both axes
    # are listed in one keyword, sep. by blanks
    'extent': ('NAXIS1', 'NAXIS2'),   # N_pixels in x/y
    'ra': 'RA',  # telescope pointing, RA
    'dec': 'DEC',  # telescope pointin, Dec
    'radec_separator': ':',   # RA/Dec hms separator, use 'XXX'
    # if already in degrees
    'date_keyword': 'DATE-OBS',  # obs date/time
    # keyword; use
    # 'date|time' if
    # separate
    'obsmidtime_jd': 'MIDTIMJD',  # obs midtime jd keyword
    # (usually provided by
    # pp_prepare
    'object': 'OBJNAME',  # object name keyword
    'filter': 'FILTER',  # filter keyword
    'filter_translations': {'SDSS u': 'u', 'SDSS g': 'g', 'SDSS r': 'r',
                              'SDSS i': 'i', 'SDSS z': 'z'},
    # filtername translation dictionary
    'exptime': 'EXPTIME',  # exposure time keyword (s)
    'airmass': 'AIRMASS',  # airmass keyword


    # source extractor settings
    'source_minarea': 7,  # default sextractor source minimum N_pixels
    'source_snr': 3,  # default sextractor source snr for registration
    'aprad_default': 4,  # default aperture radius in px
    'aprad_range': [2, 10],  # [minimum, maximum] aperture radius (px)
    'sex-config-file': rootpath+'/setup/arc35agile.sex',
    'mask_file': {},
    #                        mask files as a function of x,y binning

    # registration settings (Scamp)
    'scamp-config-file': rootpath+'/setup/arc35agile.scamp',
    'reg_max_mag': 19,
    'reg_search_radius': 0.5,  # deg
    'source_tolerance': 'high',

    # swarp settings
    'copy_keywords': ('OBSERVAT,INSTRUME,FILTER,EXPTIME,OBJNAME,' +
                      'DATE-OBS,RA,DEC,AIRMASS,SECPIX,TEL_KEYW'),
    #                        keywords to be copied in image
    #                        combination using swarp
    'swarp-config-file': rootpath+'/setup/arc35agile.swarp',

    # default catalog settings
    'astrometry_catalogs': ['GAIA'],
    'photometry_catalogs': ['SDSS-R9', 'APASS9', 'PANSTARRS', '2MASS']
}


# Magellan, IMACS long camera
magimacsl_param = {
    'telescope_instrument': 'Magellan/IMACS long',  # telescope/instrument name
    'telescope_keyword': 'MAGIMACSL',   # telescope/instrument keyword
    'observatory_code': '269',         # MPC observatory code
    'secpix': (0.11, 0.11),  # pixel size (arcsec)
    # before binning
    'ext_coeff': 0.05,          # typical extinction coefficient


    # image orientation preferences (for each chip)
    'chip_id': 'CHIP',        # chip identifier (remove,
    # if not existent)
    # the following keys are dictionaries if 'chip_id' exists, single
    # values otherwise
    'flipx': {1: True, 2: True, 3: True, 4: True, 5: True, 6: True,
                              7: True, 8: True},
    'flipy': {1: False, 2: False, 3: False, 4: False, 5: False,
              6: False, 7: False, 8: False},
    'rotate': {1: 270, 2: 270, 3: 270, 4: 270, 5: 90, 6: 90,
               7: 90, 8: 90},
    'chip_offset_fixed': {1: (-0.033, -0.099), 2: (-0.033, -0.033),
                          3: (-0.033, 0.033),  4: (-0.033, 0.099),
                          5: (0.033, -0.033),  6: (0.033, -0.099),
                          7: (0.033, 0.099),   8: (0.033, 0.033)},
    # chip offset (ra, dec in degress) [optional]

    # instrument-specific FITS header keywords
    'binning': ('BINNING#x1', 'BINNING#x2'),
    # binning in x/y, '_blankN' denotes that both axes
    # are listed in one keyword, sep. by blanks
    'extent': ('NAXIS1', 'NAXIS2'),   # N_pixels in x/y
    'ra': 'RA',  # telescope pointing, RA
    'dec': 'DEC',  # telescope pointin, Dec
    'radec_separator': ':',   # RA/Dec hms separator, use 'XXX'
    # if already in degrees
    'date_keyword': 'DATE-OBS|UT-TIME',  # obs date/time
    # keyword; use
    # 'date|time' if
                                         # separate
    'obsmidtime_jd': 'MIDTIMJD',  # obs midtime jd keyword
                                         # (usually provided by
                                         # pp_prepare
    'object': 'OBJECT',  # object name keyword
    'filter': 'FILTER',  # filter keyword
    'filter_translations': {'Sloan_u': 'u', 'Sloan_g': 'g', 'Sloan_r': 'r',
                            'Sloan_i': 'i', 'Sloan_z': 'z',
                            'Bessell_V1': 'V', 'WB4800-7800': None},
    # filtername translation dictionary
    'exptime': 'EXPTIME',  # exposure time keyword (s)
    'airmass': 'AIRMASS',  # airmass keyword


    # source extractor settings
    'source_minarea': 10,  # default sextractor source minimum N_pixels
    'source_snr': 3,  # default sextractor source snr for registration
    'aprad_default': 8,  # default aperture radius in px
    'aprad_range': [3, 15],  # [minimum, maximum] aperture radius (px)
    'sex-config-file': rootpath+'/setup/magimacs.sex',
    'mask_file': {},
    #                        mask files as a function of x,y binning

    # registration settings (Scamp)
    'scamp-config-file': rootpath+'/setup/magimacs.scamp',
    'reg_max_mag': 21,
    'reg_search_radius': 0.5,  # deg
    'source_tolerance': 'high',

    # swarp settings
    'copy_keywords': ('TELESCOP,INSTRUME,FILTER,EXPTIME,OBJECT,CHIP,EXPTYPE,' +
                      'DATE-OBS,UT-TIME,BINNING,RA,DEC,AIRMASS,' +
                      'SECPIX,TEL_KEYW'),
    #                        keywords to be copied in image
    #                        combination using swarp
    'swarp-config-file': rootpath+'/setup/magimacs.swarp',

    # default catalog settings
    'astrometry_catalogs': ['GAIA'],
    'photometry_catalogs': ['SDSS-R9', 'PANSTARRS', 'SkyMapper', 'APASS9', '2MASS']
}

# Magellan, IMACS short camera
magimacss_param = {
    'telescope_instrument': 'Magellan/IMACS short',  # telescope/instrument name
    'telescope_keyword': 'MAGIMACSS',   # telescope/instrument keyword
    'observatory_code': '269',         # MPC observatory code
    'secpix': (0.2, 0.2),  # pixel size (arcsec)
    # before binning
    'ext_coeff': 0.05,          # typical extinction coefficient


    # image orientation preferences (for each chip)
    'chip_id': 'CHIP',        # chip identifier (remove,
    # if not existent)
    # the following keys are dictionaries if 'chip_id' exists, single
    # values otherwise
    'flipx': {1: True, 2: False, 3: True, 4: True, 5: True, 6: True,
                              7: True, 8: True},
    'flipy': {1: False, 2: False, 3: False, 4: False, 5: False,
              6: False, 7: False, 8: False},
    'rotate': {1: 270, 2: 270, 3: 270, 4: 270, 5: 90, 6: 90,
               7: 90, 8: 90},
    'chip_offset_fixed': {1: (None, None), 2: (-0.11, 0.06),
                          3: (None, None),  4: (None, None),
                          5: (None, None),  6: (None, None),
                          7: (None, None),   8: (None, None)},
    # chip offset (ra, dec in degress) [optional]

    # instrument-specific FITS header keywords
    'binning': ('BINNING#x1', 'BINNING#x2'),
    # binning in x/y, '_blankN' denotes that both axes
    # are listed in one keyword, sep. by blanks
    'extent': ('NAXIS1', 'NAXIS2'),   # N_pixels in x/y
    'ra': 'RA',  # telescope pointing, RA
    'dec': 'DEC',  # telescope pointin, Dec
    'radec_separator': ':',   # RA/Dec hms separator, use 'XXX'
    # if already in degrees
    'date_keyword': 'DATE-OBS|UT-TIME',  # obs date/time
    # keyword; use
    # 'date|time' if
                                         # separate
    'obsmidtime_jd': 'MIDTIMJD',  # obs midtime jd keyword
                                         # (usually provided by
                                         # pp_prepare
    'object': 'OBJECT',  # object name keyword
    'filter': 'FILTER',  # filter keyword
    'filter_translations': {'Sloan_u': 'u', 'Sloan_g': 'g', 'Sloan_r': 'r',
                            'Sloan_i': 'i', 'Sloan_z': 'z',
                            'Bessell_V1': 'V', 'WB4800-7800': None},
    # filtername translation dictionary
    'exptime': 'EXPTIME',  # exposure time keyword (s)
    'airmass': 'AIRMASS',  # airmass keyword


    # source extractor settings
    'source_minarea': 10,  # default sextractor source minimum N_pixels
    'source_snr': 3,  # default sextractor source snr for registration
    'aprad_default': 8,  # default aperture radius in px
    'aprad_range': [3, 15],  # [minimum, maximum] aperture radius (px)
    'sex-config-file': rootpath+'/setup/magimacs.sex',
    'mask_file': {},
    #                        mask files as a function of x,y binning

    # registration settings (Scamp)
    'scamp-config-file': rootpath+'/setup/magimacs.scamp',
    'reg_max_mag': 21,
    'reg_search_radius': 0.5,  # deg
    'source_tolerance': 'high',

    # swarp settings
    'copy_keywords': ('TELESCOP,INSTRUME,FILTER,EXPTIME,OBJECT,CHIP,EXPTYPE,' +
                      'DATE-OBS,UT-TIME,BINNING,RA,DEC,AIRMASS,' +
                      'SECPIX,TEL_KEYW'),
    #                        keywords to be copied in image
    #                        combination using swarp
    'swarp-config-file': rootpath+'/setup/magimacs.swarp',

    # default catalog settings
    'astrometry_catalogs': ['GAIA'],
    'photometry_catalogs': ['SDSS-R9', 'PANSTARRS', 'SkyMapper', 'APASS9', '2MASS']
}


# Calar Alto 1.23m, DLR-MKIII
ca123dlrmkiii_param = {
    'telescope_instrument': 'Calar Alto 1.23m/DLR-MKIII',  # telescope/instrument
    'telescope_keyword': 'CA123DLRMKIII',   # telescope/instrument keyword
    'observatory_code': '493',         # MPC observatory code
    'secpix': (0.3132, 0.3132),  # pixel size (arcsec)
    # before binning
    'ext_coeff': 0.05,          # typical extinction coefficient


    # image orientation preferences
    'flipx': True,
    'flipy': True,
    'rotate': 0,

    # instrument-specific FITS header keywords
    'binning': ('CCDBINX', 'CCDBINY'),
    # binning in x/y, '_blankN' denotes that both axes
    # are listed in one keyword, sep. by blanks
    'extent': ('NAXIS1', 'NAXIS2'),   # N_pixels in x/y
    'ra': 'RA',  # telescope pointing, RA
    'dec': 'DEC',  # telescope pointin, Dec
    'radec_separator': 'XXX',   # RA/Dec hms separator, use 'XXX'
    # if already in degrees
    'date_keyword': 'DATE-OBS',  # obs date/time
    # keyword; use
    # 'date|time' if
    # separate
    'obsmidtime_jd': 'MIDTIMJD',  # obs midtime jd keyword
    # (usually provided by
    # pp_prepare
    'object': 'OBJECT',  # object name keyword
    'filter': 'FILTER',  # filter keyword
    'filter_translations': {'V_Johnson': 'V', 'R_Johnson': 'R', 'free': None},
    # filtername translation dictionary
    'exptime': 'EXPTIME',  # exposure time keyword (s)
    'airmass': 'AIRMASS',  # airmass keyword


    # source extractor settings
    'source_minarea': 15,  # default sextractor source minimum N_pixels
    'source_snr': 3,  # default sextractor source snr for registration
    'aprad_default': 4,  # default aperture radius in px
    'aprad_range': [2, 15],  # [minimum, maximum] aperture radius (px)
    'sex-config-file': rootpath+'/setup/ca123dlrmkiii.sex',
    'mask_file': {},
    #                        mask files as a function of x,y binning

    # registration settings (Scamp)
    'scamp-config-file': rootpath+'/setup/ca123dlrmkiii.scamp',
    'reg_max_mag': 19,
    'reg_search_radius': 0.5,  # deg
    'source_tolerance': 'high',

    # swarp settings
    'copy_keywords': ('TELESCOP,INSTRUME,FILTER,EXPTIME,OBJECT,' +
                      'DATE-OBS,RA,DEC,AIRMASS,' +
                      'SECPIX,TEL_KEYW'),
    #                        keywords to be copied in image
    #                        combination using swarp
    'swarp-config-file': rootpath+'/setup/ca123dlrmkiii.swarp',

    # default catalog settings
    'astrometry_catalogs': ['GAIA'],
    'photometry_catalogs': ['SDSS-R9', 'APASS9', 'PANSTARRS', '2MASS']
}


# Lowell31, NASACAM
lowell31_param = {
    'telescope_instrument': 'Lowell31/NASACAM',  # telescope/instrument name
    'telescope_keyword': 'LOWELL31',      # telescope/instrument keyword
    'observatory_code': '688',         # MPC observatory code
    'secpix': (0.456, 0.456),  # pixel size (arcsec)
    # before binning
    'ext_coeff': 0.05,          # typical extinction coefficient

    # image orientation preferences
    'flipx': True,
    'flipy': False,
    'rotate': 270,

    # instrument-specific FITS header keywords
    #    'binning'              : ('CCDSUM#blank1', 'CCDSUM#blank2'), # binning in x/y
    'binning': (1, 1),  # binning in x/y
    'extent': ('NAXIS1', 'NAXIS2'),   # N_pixels in x/y
    'ra': 'TELRA',  # telescope pointing, RA
    'dec': 'TELDEC',  # telescope pointin, Dec
    #'ra'                   : 'RA',  # telescope pointing, RA
    #'dec'                  : 'DEC', # telescope pointin, Dec
    'radec_separator': ':',   # RA/Dec hms separator, use 'XXX'
    # if already in degrees
    'date_keyword': 'DATE-OBS|UT',  # obs date/time
    # keyword; use
    # 'date|time' if
    # separate
    'obsmidtime_jd': 'MIDTIMJD',  # obs midtime jd keyword
    # (usually provided by
    # pp_prepare
    'object': 'OBJECT',  # object name keyword
    'filter': 'FILTER2',  # filter keyword
    'filter_translations': {'V': 'V', 'I': 'I', 'R': 'R',
                            'clear': None},
    # filtername translation dictionary
    'exptime': 'EXPTIME',  # exposure time keyword (s)
    'airmass': 'AIRMASS',  # airmass keyword


    # source extractor settings
    'source_minarea': 12,  # default sextractor source minimum N_pixels
    'source_snr': 3,  # default sextractor source snr for registration
    'aprad_default': 5,  # default aperture radius in px
    'aprad_range': [2, 10],  # [minimum, maximum] aperture radius (px)
    'sex-config-file': rootpath+'/setup/lowell31.sex',
    'mask_file': {},
    #                        mask files as a function of x,y binning

    # registration settings (Scamp)
    'scamp-config-file': rootpath+'/setup/lowell31.scamp',
    'reg_max_mag': 19,
    'reg_search_radius': 0.5,  # deg
    'source_tolerance': 'high',

    # swarp settings
    'copy_keywords': ('TELESCOP,INSTRUME,FILTER1,FILTER2,EXPTIME,OBJECT,' +
                      'DATE-OBS,UT,TELRA,TELDEC,SCALE,AIRMASS,' +
                      'TEL_KEYW'),
    #                         keywords to be copied in image
    #                         combination using swarp
    'swarp-config-file': rootpath+'/setup/lowell31.swarp',

    # default catalog settings
    'astrometry_catalogs': ['GAIA'],
    'photometry_catalogs': ['SDSS-R9', 'APASS9', 'PANSTARRS', '2MASS']
}


# Lowell42, NASA42
lowell42_param = {
    'telescope_instrument': 'Lowell42/NASA42',  # telescope/instrument name
    'telescope_keyword': 'LOWELL42',      # telescope/instrument keyword
    'observatory_code': '688',         # MPC observatory code
    'secpix': (0.327, 0.327),  # pixel size (arcsec)
    # before binning
    'ext_coeff': 0.05,          # typical extinction coefficient


    # image orientation preferences
    'flipx': True,
    'flipy': False,
    'rotate': 90,

    # instrument-specific FITS header keywords
    'binning': ('CCDSUM#blank1', 'CCDSUM#blank2'),
    # binning in x/y, '_blankN' denotes that both axes
    # are listed in one keyword, sep. by blanks
    'extent': ('NAXIS1', 'NAXIS2'),   # N_pixels in x/y
    'ra': 'TELRA',  # telescope pointing, RA
    'dec': 'TELDEC',  # telescope pointin, Dec
    'radec_separator': ':',   # RA/Dec hms separator, use 'XXX'
    # if already in degrees
    'date_keyword': 'DATE-OBS|UTC-OBS',  # obs date/time
    # keyword; use
    # 'date|time' if
    # separate
    'obsmidtime_jd': 'MIDTIMJD',  # obs midtime jd keyword
    # (usually provided by
                                         # pp_prepare
    'object': 'OBJECT',  # object name keyword
    'filter': 'FILTNAME',  # filter keyword
    'filter_translations': {'V': 'V', 'I': 'I', 'VR': None,
                            'R': 'R', 'clear': None, 'CN': None,
                            'open': None, 'Open': None},
    # filtername translation dictionary
    'exptime': 'EXPTIME',  # exposure time keyword (s)
    'airmass': 'AIRMASS',  # airmass keyword


    # source extractor settings
    'source_minarea': 12,  # default sextractor source minimum N_pixels
    'source_snr': 3,  # default sextractor source snr for registration
    'aprad_default': 5,  # default aperture radius in px
    'aprad_range': [2, 10],  # [minimum, maximum] aperture radius (px)
    'sex-config-file': rootpath+'/setup/lowell42.sex',
    'mask_file': {},
    #                        mask files as a function of x,y binning

    # registration settings (Scamp)
    'scamp-config-file': rootpath+'/setup/lowell42.scamp',
    'reg_max_mag': 19,
    'reg_search_radius': 0.5,  # deg
    'source_tolerance': 'high',

    # swarp settings
    'copy_keywords': ('TELESCOP,INSTRUME,CCDSUM,FILTNAME,EXPTIME,' +
                      'OBJECT,' +
                      'DATE-OBS,UTC-OBS,TELRA,TELDEC,PIXSCAL,AIRMASS,' +
                      'TEL_KEYW'),
    #                         keywords to be copied in image
    #                         combination using swarp
    'swarp-config-file': rootpath+'/setup/lowell42.swarp',

    # default catalog settings
    'astrometry_catalogs': ['GAIA'],
    'photometry_catalogs': ['PANSTARRS', 'SDSS-R9', 'APASS9', '2MASS']
}

# Lowell42, SITE
lowell42site_param = {
    'telescope_instrument': 'Lowell42/SITE',  # telescope/instrument name
    'telescope_keyword': 'LOWELL42SITE',      # telescope/instrument keyword
    'observatory_code': '688',         # MPC observatory code
    'secpix': (0.6, 0.6),  # pixel size (arcsec)
    # before binning
    'ext_coeff': 0.05,          # typical extinction coefficient


    # image orientation preferences
    'flipx': True,
    'flipy': True,
    'rotate': 0,

    # instrument-specific FITS header keywords
    'binning': ('CCDSUM#blank1', 'CCDSUM#blank2'),
    # binning in x/y, '_blankN' denotes that both axes
    # are listed in one keyword, sep. by blanks
    'extent': ('NAXIS1', 'NAXIS2'),   # N_pixels in x/y
    'ra': 'RA',  # telescope pointing, RA
    'dec': 'DEC',  # telescope pointin, Dec
    'radec_separator': ':',   # RA/Dec hms separator, use 'XXX'
    # if already in degrees
    'date_keyword': 'DATE-OBS|UT',  # obs date/time
    # keyword; use
    # 'date|time' if
    # separate
    'obsmidtime_jd': 'MIDTIMJD',  # obs midtime jd keyword
    # (usually provided by
    # pp_prepare
    'object': 'OBJECT',  # object name keyword
    'filter': 'FILTNAME',  # filter keyword
    'filter_translations': {'R': 'R', 'I': 'I', 'V': 'V'},
    # filtername translation dictionary
    'exptime': 'EXPTIME',  # exposure time keyword (s)
    'airmass': 'AIRMASS',  # airmass keyword


    # source extractor settings
    'source_minarea': 12,  # default sextractor source minimum N_pixels
    'source_snr': 3,  # default sextractor source snr for registration
    'aprad_default': 5,  # default aperture radius in px
    'aprad_range': [2, 10],  # [minimum, maximum] aperture radius (px)
    'sex-config-file': rootpath+'/setup/lowell42.sex',
    'mask_file': {},
    #                        mask files as a function of x,y binning

    # registration settings (Scamp)
    'scamp-config-file': rootpath+'/setup/lowell42.scamp',
    'reg_max_mag': 19,
    'reg_search_radius': 0.5,  # deg
    'source_tolerance': 'high',

    # swarp settings
    'copy_keywords': ('TELESCOP,INSTRUME,CCDSUM,FILTNAME,EXPTIME,' +
                      'OBJECT,' +
                      'DATE-OBS,UT,RA,DEC,PIXSCAL,AIRMASS,' +
                      'TEL_KEYW'),
    #                         keywords to be copied in image
    #                         combination using swarp
    'swarp-config-file': rootpath+'/setup/lowell42.swarp',

    # default catalog settings
    'astrometry_catalogs': ['GAIA'],
    'photometry_catalogs': ['SDSS-R9', 'APASS9', 'PANSTARRS', '2MASS']
}


# Lowell72 (Perkins), PRISM
lowell72_param = {
    'telescope_instrument': 'Lowell72/PRISM',  # telescope/instrument name
    'telescope_keyword': 'LOWELL72',      # telescope/instrument keyword
    'observatory_code': '688',         # MPC observatory code
    'secpix': (0.39, 0.39),  # pixel size (arcsec)
    # before binning
    'ext_coeff': 0.05,          # typical extinction coefficient


    # image orientation preferences
    'flipx': True,
    'flipy': False,
    'rotate': 0,

    # instrument-specific FITS header keywords
    'binning': ('CCDSUM#blank1', 'CCDSUM#blank2'),
    # binning in x/y, '_blankN' denotes that both axes
    # are listed in one keyword, sep. by blanks
    'extent': ('NAXIS1', 'NAXIS2'),   # N_pixels in x/y
    'ra': 'TELRA',  # telescope pointing, RA
    'dec': 'TELDEC',  # telescope pointin, Dec
    'radec_separator': ':',   # RA/Dec hms separator, use 'XXX'
    # if already in degrees
    'date_keyword': 'DATE-OBS',  # obs date/time
    # keyword; use
    # 'date|time' if
    # separate
    'obsmidtime_jd': 'MIDTIMJD',  # obs midtime jd keyword
    # (usually provided by
    # pp_prepare
    'object': 'OBJECT',  # object name keyword
    'filter': 'FILTNME3',  # filter keyword
    'filter_translations': {'B': 'B', 'V': 'V', 'I': 'I', 'R': 'R', 'VR': None},
    # filtername translation dictionary
    'exptime': 'EXPTIME',  # exposure time keyword (s)
    'airmass': 'AIRMASS',  # airmass keyword


    # source extractor settings
    'source_minarea': 5,  # default sextractor source minimum N_pixels
    'source_snr': 3,  # default sextractor source snr for registration
    'aprad_default': 2,  # default aperture radius in px
    'aprad_range': [1, 7],  # [minimum, maximum] aperture radius (px)
    'sex-config-file': rootpath+'/setup/lowell72.sex',
    'mask_file': {},
    #                        mask files as a function of x,y binning

    # registration settings (Scamp)
    'scamp-config-file': rootpath+'/setup/lowell72.scamp',
    'reg_max_mag': 19,
    'reg_search_radius': 0.5,  # deg
    'source_tolerance': 'high',

    # swarp settings
    'copy_keywords': ('TELESCOP,INSTRUME,FILTNME3,' +
                      'EXPTIME,OBJECT,' +
                      'DATE-OBS,TELRA,TELDEC,PIXSCAL,AIRMASS,' +
                      'TEL_KEYW'),
    #                         keywords to be copied in image
    #                         combination using swarp
    'swarp-config-file': rootpath+'/setup/lowell72.swarp',

    # default catalog settings
    'astrometry_catalogs': ['GAIA'],
    'photometry_catalogs': ['PANSTARRS', 'SDSS-R9', 'APASS9', '2MASS']
}

# CTIO 0.9m, CFCCD
ctio09_param = {
    'telescope_instrument': 'CTIO09/CFCCD',  # telescope/instrument name
    'telescope_keyword': 'CTIO09',      # telescope/instrument keyword
    'observatory_code': '807',         # MPC observatory code
    'secpix': (0.396, 0.396),  # pixel size (arcsec)
    # before binning
    'ext_coeff': 0.05,          # typical extinction coefficient


    # image orientation preferences
    'flipx': False,
    'flipy': False,
    'rotate': 180,

    # instrument-specific FITS header keywords
    'binning': ('CCDSUM#blank1', 'CCDSUM#blank2'),
    # binning in x/y
    'extent': ('NAXIS1', 'NAXIS2'),   # N_pixels in x/y
    'ra': 'RA',  # telescope pointing, RA
    'dec': 'DEC',  # telescope pointin, Dec
    'radec_separator': ':',   # RA/Dec hms separator, use 'XXX'
    # if already in degrees
    'date_keyword': 'DATE-OBS',  # obs date/time
    # keyword; use
    # 'date|time' if
    # separate
    'obsmidtime_jd': 'MIDTIMJD',  # obs midtime jd keyword
    # (usually provided by
    # pp_prepare
    'object': 'OBJECT',  # object name keyword
    'filter': 'FILTERS',  # filter keyword
    'filter_translations': {'dia v': 'V', 'dia ov': 'V', 'dia i': 'I',
                              'dia r': 'R'},
    # filtername translation dictionary
    'exptime': 'EXPTIME',  # exposure time keyword (s)
    'airmass': 'AIRMASS',  # airmass keyword


    # source extractor settings
    'source_minarea': 15,  # default sextractor source minimum N_pixels
    'source_snr': 3,  # default sextractor source snr for registration
    'aprad_default': 5,  # default aperture radius in px
    'aprad_range': [2, 10],  # [minimum, maximum] aperture radius (px)
    'sex-config-file': rootpath+'/setup/ctio09.sex',
    'mask_file': {},
    #                        mask files as a function of x,y binning

    # registration settings (Scamp)
    'scamp-config-file': rootpath+'/setup/ctio09.scamp',
    'reg_max_mag': 19,
    'reg_search_radius': 0.5,  # deg
    'source_tolerance': 'high',

    # swarp settings
    'copy_keywords': ('TELESCOP,INSTRUME,FILTERS,EXPTIME,OBJECT,' +
                      'DATE-OBS,RA,DEC,CCDSUM,AIRMASS,' +
                      'TEL_KEYW'),
    #                         keywords to be copied in image
    #                         combination using swarp
    'swarp-config-file': rootpath+'/setup/ctio09.swarp',

    # default catalog settings
    'astrometry_catalogs': ['GAIA'],
    'photometry_catalogs': ['SDSS-R9', 'APASS9', 'PANSTARRS', '2MASS']
}


# CTIO 1.0m, Y4KCAM
ctio10_param = {
    'telescope_instrument': 'CTIO10/Y4KCAM',  # telescope/instrument name
    'telescope_keyword': 'CTIO10',      # telescope/instrument keyword
    'observatory_code': '807',         # MPC observatory code
    'secpix': (0.289, 0.289),  # pixel size (arcsec)
    # before binning
    'ext_coeff': 0.05,          # typical extinction coefficient


    # image orientation preferences
    'flipx': True,
    'flipy': False,
    'rotate': 180,

    # instrument-specific FITS header keywords
    'binning': ('CCDSUM#blank1', 'CCDSUM#blank2'),
    # binning in x/y
    'extent': ('NAXIS1', 'NAXIS2'),   # N_pixels in x/y
    'ra': 'RA',  # telescope pointing, RA
    'dec': 'DEC',  # telescope pointin, Dec
    'radec_separator': ':',   # RA/Dec hms separator, use 'XXX'
    # if already in degrees
    'date_keyword': 'DATE-OBS|TIME-OBS',  # obs date/time
    # keyword; use
    # 'date|time' if
    # separate
    'obsmidtime_jd': 'MIDTIMJD',  # obs midtime jd keyword
    # (usually provided by
    # pp_prepare
    'object': 'OBJECT',  # object name keyword
    'filter': 'FILTERID',  # filter keyword
    'filter_translations': {'V': 'V', 'I': 'I'},
    # filtername translation dictionary
    'exptime': 'EXPTIME',  # exposure time keyword (s)
    'airmass': 'AIRMASS',  # airmass keyword


    # source extractor settings
    'source_minarea': 12,  # default sextractor source minimum N_pixels
    'source_snr': 3,  # default sextractor source snr for registration
    'aprad_default': 5,  # default aperture radius in px
    'aprad_range': [2, 10],  # [minimum, maximum] aperture radius (px)
    'sex-config-file': rootpath+'/setup/ctio10.sex',
    'mask_file': {},
    #                        mask files as a function of x,y binning

    # registration settings (Scamp)
    'scamp-config-file': rootpath+'/setup/ctio10.scamp',
    'reg_max_mag': 19,
    'reg_search_radius': 0.5,  # deg
    'source_tolerance': 'high',

    # swarp settings
    'copy_keywords': ('TELESCOP,INSTRUME,FILTER,EXPTIME,OBJECT,' +
                      'DATE-OBS,TIME-OBS, RA,DEC,CCDSUM,AIRMASS,' +
                      'TEL_KEYW'),
    #                         keywords to be copied in image
    #                         combination using swarp
    'swarp-config-file': rootpath+'/setup/ctio10.swarp',

    # default catalog settings
    'astrometry_catalogs': ['GAIA'],
    'photometry_catalogs': ['SDSS-R9', 'APASS9', 'PANSTARRS', '2MASS']
}


# CTIO 1.3m, ANDICAM (CCD)
ctio13ccd_param = {
    'telescope_instrument': 'CTIO/ANDICAM_CCD',  # telescope/instrument name
    'telescope_keyword': 'CTIO13CCD',        # telescope/instrument keyword
    'observatory_code': '807',         # MPC observatory code
    'secpix': (0.185, 0.185),  # pixel size (arcsec)
    # before binning
    'ext_coeff': 0.05,          # typical extinction coefficient


    # image orientation preferences
    'flipx': True,
    'flipy': False,
    'rotate': 0,

    # instrument-specific FITS header keywords
    'binning': ('CCDXBIN', 'CCDYBIN'),  # binning in x/y
    'extent': ('NAXIS1', 'NAXIS2'),   # N_pixels in x/y
    'ra': 'RA',  # telescope pointing, RA
    'dec': 'DEC',  # telescope pointin, Dec
    'radec_separator': ':',   # RA/Dec hms separator, use 'XXX'
    # if already in degrees
    'date_keyword': 'DATE-OBS|TIME-OBS',  # obs date/time
    # keyword; use
    # 'date|time' if
    # separate
    'obsmidtime_jd': 'MIDTIMJD',  # obs midtime jd keyword
    # (usually provided by
    # pp_prepare
    'object': 'OBJECT',  # object name keyword
    'filter': 'CCDFLTID',  # filter keyword
    'filter_translations': {'V': 'V', 'R': 'R', 'B': 'B'},
    # filtername translation dictionary
    'exptime': 'EXPTIME',  # exposure time keyword (s)
    'airmass': 'SECZ',  # airmass keyword


    # source extractor settings
    'source_minarea': 15,  # default sextractor source minimum N_pixels
    'source_snr': 3,  # default sextractor source snr for registration
    'aprad_default': 4,  # default aperture radius in px
    'aprad_range': [2, 10],  # [minimum, maximum] aperture radius (px)
    'sex-config-file': rootpath+'/setup/andicam.sex',
    'mask_file': {'2,2': rootpath+'/setup/mask_andicam_2x2.fits'},
    #                        mask files as a function of x,y binning

    # registration settings (Scamp)
    'scamp-config-file': rootpath+'/setup/andicam.scamp',
    'reg_max_mag': 19,
    'reg_search_radius': 0.5,  # deg
    'source_tolerance': 'high',

    # swarp settings
    'copy_keywords': ('OBSERVAT,INSTRUME,CCDFLTID,EXPTIME,OBJECT,' +
                      'DATE-OBS,TIME-OBS,RA,DEC,SECPIX,SECZ,' +
                      'TEL_KEYW'),
    #                        keywords to be copied in image
    #                        combination using swarp
    'swarp-config-file': rootpath+'/setup/andicam.swarp',

    # default catalog settings
    'astrometry_catalogs': ['GAIA'],
    'photometry_catalogs': ['SDSS-R9', 'APASS9', 'PANSTARRS', '2MASS']
}


# UH 88inch, SNIFS (imaging mode)
uh88snifs_param = {
    'telescope_instrument': 'UH88/SNIFS',  # telescope/instrument name
    'telescope_keyword': 'UH88SNIFS',        # telescope/instrument keyword
    'observatory_code': '568',         # MPC observatory code
    'secpix': (0.1365, 0.1365),  # pixel size (arcsec)
    # before binning
    'ext_coeff': 0.05,          # typical extinction coefficient


    # image orientation preferences
    'flipx': True,
    'flipy': True,
    'rotate': 0,

    # instrument-specific FITS header keywords
    'binning': ('CCDBIN1', 'CCDBIN2'),  # binning in x/y
    'extent': ('NAXIS1', 'NAXIS2'),   # N_pixels in x/y
    'ra': 'RA',  # telescope pointing, RA
    'dec': 'DEC',  # telescope pointin, Dec
    'radec_separator': ':',   # RA/Dec hms separator, use 'XXX'
    # if already in degrees
    'date_keyword': 'DATE-OBS',  # obs date/time
    # keyword; use
    # 'date|time' if
    # separate
    'obsmidtime_jd': 'MIDTIMJD',  # obs midtime jd keyword
    # (usually provided by
    # pp_prepare
    'object': 'OBJECT',  # object name keyword
    'filter': 'FILTER',  # filter keyword
    'filter_translations': {'SDSS u': 'u', 'SDSS g': 'g', 'SDSS r': 'r',
                              'SDSS i': 'i', 'SDSS z': 'z', 'Bessell B': 'B',
                              'Bessell V': 'V', 'Bessell R': 'R',
                              'Bessell I': 'I'},
    # filtername translation dictionary
    'exptime': 'EXPTIME',  # exposure time keyword (s)
    'airmass': 'AIRMASS',  # airmass keyword


    # source extractor settings
    'source_minarea': 12,  # default sextractor source minimum N_pixels
    'source_snr': 3,  # default sextractor source snr for registration
    'aprad_default': 3,  # default aperture radius in px
    'aprad_range': [2, 10],  # [minimum, maximum] aperture radius (px)
    'sex-config-file': rootpath+'/setup/uh88snifs.sex',
    'mask_file': {'2,2': rootpath+'/setup/mask_snifs_2x2.fits'},
    #                        mask files as a function of x,y binning

    # registration settings (Scamp)
    'scamp-config-file': rootpath+'/setup/uh88snifs.scamp',
    'source_tolerance': 'high',

    # swarp settings
    'copy_keywords': ('OBSERVAT,INSTRUME,FILTER,EXPTIME,OBJECT,' +
                      'DATE-OBS,RA,DEC,SECPIX,AIRMASS,' +
                      'TEL_KEYW'),
    #                        keywords to be copied in image
    #                        combination using swarp
    'swarp-config-file': rootpath+'/setup/uh88snifs.swarp',
    'reg_max_mag': 19,
    'reg_search_radius': 0.5,  # deg

    # default catalog settings
    'astrometry_catalogs': ['GAIA'],
    'photometry_catalogs': ['SDSS-R9', 'APASS9', 'PANSTARRS', '2MASS']
}

# WIYN 0.9m, Half Degree Imager (HDI)
wiyn09hdi_param = {
    'telescope_instrument': 'WIYN09/HDI',  # telescope/instrument name
    'telescope_keyword': 'WIYN09HDI',  # telescope/instrument keyword
    'observatory_code': '695',         # MPC observatory code
    'secpix': (0.43, 0.43),  # pixel size (arcsec)
    # before binning
    'ext_coeff': 0.05,          # typical extinction coefficient


    # image orientation preferences
    'flipx': False,
    'flipy': False,
    'rotate': 0,

    # instrument-specific FITS header keywords
    'binning': ('CCDBIN1', 'CCDBIN2'),  # binning in x/y
    'extent': ('NAXIS1', 'NAXIS2'),   # N_pixels in x/y
    'ra': 'RASTRNG',  # telescope pointing, RA
    'dec': 'DECSTRNG',  # telescope pointin, Dec
    'radec_separator': ':',   # RA/Dec hms separator, use 'XXX'
    # if already in degrees
    'date_keyword': 'DATE-OBS',  # obs date/time
    # keyword; use
    # 'date|time' if
    # separate
    'obsmidtime_jd': 'MIDTIMJD',  # obs midtime jd keyword
    # (usually provided by
    # pp_prepare
    'object': 'OBJECT',  # object name keyword
    'filter': 'FILTER1',  # filter keyword
    'filter_translations': {'u': 'u', 'g': 'g', 'r': 'r',
                              'i': 'i', 'z': 'z'},
    # filtername translation dictionary
    'exptime': 'EXPTIME',  # exposure time keyword (s)
    'airmass': 'AIRMASS',  # airmass keyword


    # source extractor settings
    'source_minarea': 9,  # default sextractor source minimum N_pixels
    'source_snr': 3,  # default sextractor source snr for registration
    'aprad_default': 3,  # default aperture radius in px
    'aprad_range': [2, 8],  # [minimum, maximum] aperture radius (px)
    'sex-config-file': rootpath+'/setup/wiyn09hdi.sex',
    'mask_file': {},
    #                        mask files as a function of x,y binning

    # registration settings (Scamp)
    'scamp-config-file': rootpath+'/setup/wiyn09hdi.scamp',
    'reg_max_mag': 19,
    'reg_search_radius': 0.5,  # deg
    'source_tolerance': 'high',

    # swarp settings
    'copy_keywords': ('OBSERVAT,INSTRUME,FILTER1,EXPTIME,OBJECT,' +
                      'DATE-OBS,RASTRNG,DECSTRNG,SECPIX,AIRMASS,' +
                      'TEL_KEYW'),
    #                        keywords to be copied in image
    #                        combination using swarp
    'swarp-config-file': rootpath+'/setup/wiyn09hdi.swarp',

    # default catalog settings
    'astrometry_catalogs': ['GAIA'],
    'photometry_catalogs': ['SDSS-R9', 'APASS9', 'PANSTARRS', '2MASS']
}


# Generic telescope (e.g., amateur telescope)
generic_param = {
    'telescope_instrument': 'Generic',  # telescope/instrument name
    'telescope_keyword': 'GENERIC',     # telescope/instrument keyword
    'observatory_code': '500',         # MPC observatory code
    'secpix': (None, None),  # pixel size (arcsec)
    # before binning
    'ext_coeff': 0.05,          # typical extinction coefficient


    # image orientation preferences
    'flipx': False,
    'flipy': False,
    'rotate': 0,

    # instrument-specific FITS header keywords
    # in this GENERIC setup, all keywords are suggested
    # and will be checked in pp_prepare
    'binning': ('BINX', 'BINY'),  # binning in x/y
    'extent': ('NAXIS1', 'NAXIS2'),   # N_pixels in x/y
    'ra': 'OBJCTRA',  # telescope pointing, RA
    'dec': 'OBJCTDEC',  # telescope pointin, Dec
    'radec_separator': ' ',   # RA/Dec hms separator, use 'XXX'
    # if already in degrees
    'date_keyword': 'DATE-OBS',  # obs date/time
    # keyword; use
    # 'date|time' if
    # separate
    'obsmidtime_jd': 'MIDTIMJD',  # obs midtime jd keyword
    # (usually provided by
    # pp_prepare
    'object': 'OBJECT',  # object name keyword
    'filter': 'FILTER',  # filter keyword
    'filter_translations': {'V': 'V', 'R': 'R', 'clear': None, '': None},
    # filtername translation dictionary
    'exptime': 'EXPTIME',  # exposure time keyword (s)
    'airmass': 'AIRMASS',  # airmass keyword


    # source extractor settings
    'source_minarea': 8,  # default sextractor source minimum N_pixels
    'source_snr': 3,  # default sextractor source snr for registration
    'aprad_default': 3,  # default aperture radius in px
    'aprad_range': [2, 10],  # [minimum, maximum] aperture radius (px)
    'sex-config-file': rootpath+'/setup/generic.sex',
    'mask_file': {},  # '2,2' : rootpath+'/setup/mask_snifs_2x2.fits'},
    #                        mask files as a function of x,y binning

    # registration settings (Scamp)
    'scamp-config-file': rootpath+'/setup/generic.scamp',
    'reg_max_mag': 19,
    'reg_search_radius': 0.5,  # deg
    'source_tolerance': 'high',

    # swarp settings
    'copy_keywords': ('OBSERVAT,INSTRUME,EXPTIME,OBJECT,' +
                      'DATE-OBS,TEL_KEYW'),
    #                        keywords to be copied in image
    #                        combination using swarp
    'swarp-config-file': rootpath+'/setup/generic.swarp',

    # default catalog settings
    'astrometry_catalogs': ['GAIA'],
    'photometry_catalogs': ['SDSS-R9', 'APASS9', 'PANSTARRS', '2MASS']
}

# SPM 1.3m, RATIR
ratir_param = {
    'telescope_instrument': 'RATIR',  # telescope/instrument name
    'telescope_keyword': 'RATIR',        # telescope/instrument keyword
    'observatory_code': '679',         # MPC observatory code
    'secpix': (0.15, 0.15),  # pixel size (arcsec)
    # before binning
    'ext_coeff': 0.05,          # typical extinction coefficient


    # image orientation preferences
    'flipx': True,
    'flipy': False,
    'rotate': 0,

    # instrument-specific FITS header keywords
    'binning': ('SC#CH#DTBN', 'SC#CH#DTBN'),  # binning in x/y
    # _CH_ gets replaced with Channel number
    'extent': ('NAXIS1', 'NAXIS2'),   # N_pixels in x/y
    'ra': 'STRSTRA',  # telescope pointing, RA
    'dec': 'STRSTDE',  # telescope pointin, Dec
    'radec_separator': 'XXX',   # RA/Dec hms separator, use 'XXX'
    # if already in degrees
    'date_keyword': 'SSHT',  # obs date/time
    # keyword; use
    # 'date|time' if
    # separate
    'obsmidtime_jd': 'MIDTIMJD',  # obs midtime jd keyword
    # (usually provided by
    # pp_prepare
    'object': 'OBJECT',  # object name keyword
    'filter': 'FILTER',  # filter keyword
    'filter_translations': {'g': 'g', 'r': 'r', 'i': 'i', 'z': 'z',
                            'J': 'J', 'H': 'H'},
    # filtername translation dictionary
    'exptime': 'EXPTIME',  # exposure time keyword (s)
    'airmass': 'STROBAM',  # airmass keyword


    # source extractor settings
    'source_minarea': 18,  # default sextractor source minimum N_pixels
    'source_snr': 3,  # default sextractor source snr for registration
    'aprad_default': 4,  # default aperture radius in px
    'aprad_range': [2, 10],  # [minimum, maximum] aperture radius (px)
    'sex-config-file': rootpath+'/setup/ratir.sex',
    'mask_file': {},
    #                        mask files as a function of x,y binning

    # registration settings (Scamp)
    'scamp-config-file': rootpath+'/setup/ratir.scamp',
    'reg_max_mag': 19,
    'reg_search_radius': 0.5,  # deg
    'source_tolerance': 'high',

    # swarp settings
    'copy_keywords': ('TELESCOP,INSTRUME,FILTER,EXPTIME,OBJECT' +
                      'SSHT,STRSTRA,STRSTDE,SECPIX,STROBAM,' +
                      'TEL_KEYW'),
    #                        keywords to be copied in image
    #                        combination using swarp
    'swarp-config-file': rootpath+'/setup/ratir.swarp',

    # default catalog settings
    'astrometry_catalogs': ['GAIA'],
    'photometry_catalogs': ['SDSS-R9', 'APASS9', 'PANSTARRS', '2MASS']
}

sl40in_param = {
    'telescope_instrument': 'Sutherland 40inch/SHOC',  # telescope/instrument name
    'telescope_keyword': 'SL40IN',  # telescope/instrument keyword
    'observatory_code': 'K94',  # MPC observatory code
    'secpix': (0.167, 0.167),  # pixel size (arcsec) before binning

    # image orientation preferences
    'flipx': False,
    'flipy': False,
    'rotate': 0,

    # instrument-specific FITS header keywords
    'binning': ('HBIN', 'VBIN'),  # binning in x/y
    'extent': ('NAXIS1', 'NAXIS2'),  # N_pixels in x/y
    'ra': 'TELRA',  # telescope pointing, RA
    'dec': 'TELDEC',  # telescope pointin, Dec
    'radec_separator': ':',  # RA/Dec hms separator, use 'XXX'
    # if already in degrees
    'date_keyword': 'FRAME',  # obs date/time
    # keyword; use
    # 'date|time' if
    # separate
    'obsmidtime_jd': 'MJD-OBS',  # obs midtime jd keyword
    # (usually provided by
    # pp_prepare
    'object': 'OBJECT',  # object name keyword
    # 'filter': 'FILTERA',  # filter keyword
    # 'filter_translations': {'V - Green': 'V','R - Red': 'R','I - Infrared': 'I'},
    'filter': 'FILTERB',  # filter keyword
    'filter_translations': {'u\'': 'u', 'g\'': 'g', 'r\'': 'r', 'i\'': 'i', 'z\'': 'z'},
    # filtername translation dictionary
    'exptime': 'EXPOSURE',  # exposure time keyword (s)
    'airmass': 'AIRMASS',  # airmass keyword

    # source extractor settings
    'source_minarea': 12,  # default sextractor source minimum N_pixels
    'source_snr': 3,  # default sextractor source snr for registration
    'aprad_default': 5,  # default aperture radius in px
    'aprad_range': [2, 10],  # [minimum, maximum] aperture radius (px)
    'sex-config-file': rootpath + '/setup/sl40in.sex',
    'mask_file': {},
    #                        mask files as a function of x,y binning

    # scamp settings
    'scamp-config-file': rootpath + '/setup/sl40in.scamp',
    'reg_max_mag': 19,
    'reg_search_radius': 0.5,  # deg
    'source_tolerance': 'high',

    # default catalog settings
    'astrometry_catalogs': ['GAIA'],
    'photometry_catalogs': ['SkyMapper', 'SDSS-R9', 'APASS9', 'PANSTARRS', '2MASS']
}


# SOAR, Goodman [old] (imager)
# running Goodman data requires the removal of header keywords
# PARAM0, PARAM61, PARAM62, PARAM63 (degree symbol is non-ASCII)
soargoodmanold_param = {
    'telescope_instrument': 'SOAR/GOODMANOLD',  # telescope/instrument name
    'telescope_keyword': 'SOARGOODMANOLD',  # telescope/instrument keyword
    'observatory_code': 'I33',         # MPC observatory code
    'secpix': (0.15, 0.15),  # pixel size (arcsec)
    # before binning
    'ext_coeff': 0.05,          # typical extinction coefficient


    # image orientation preferences
    'flipx': True,
    'flipy': True,
    'rotate': 0,

    # instrument-specific FITS header keywords
    'binning': ('PARAM18', 'PARAM22'),  # binning in x/y
    'extent': ('NAXIS1', 'NAXIS2'),   # N_pixels in x/y
    'ra': 'RA',  # telescope pointing, RA
    'dec': 'DEC',  # telescope pointin, Dec
    'radec_separator': ':',   # RA/Dec hms separator, use 'XXX'
    # if already in degrees
    'date_keyword': 'DATE-OBS',  # obs date/time
    # keyword; use
    # 'date|time' if
    # separate
    'obsmidtime_jd': 'MIDTIMJD',  # obs midtime jd keyword
    # (usually provided by
    # pp_prepare
    'object': 'OBJECT',  # object name keyword
    'filter': 'FILTER',  # filter keyword
    'filter_translations': {'Rc': 'R', 'R-Bessel': 'R',
                            'V': 'V', 'B': 'B', 'u': 'u',
                              'g-SDSS': 'g', 'r-SDSS': 'r', 'i-SDSS': 'i',
                              'z-SDSS': 'z', 'VR': None},
    # filtername translation dictionary
    'exptime': 'EXPTIME',  # exposure time keyword (s)
    'airmass': 'AIRMASS',  # airmass keyword


    # source extractor settings
    'source_minarea': 15,  # default sextractor source minimum N_pixels
    'source_snr': 3,  # default sextractor source snr for registration
    'aprad_default': 4,  # default aperture radius in px
    'aprad_range': [2, 10],  # [minimum, maximum] aperture radius (px)
    'sex-config-file': rootpath+'/setup/soargoodman.sex',
    'mask_file': {'1,1': rootpath+'/setup/mask_soargoodman_1x1.fits',
                  '2,2': rootpath+'/setup/mask_soargoodman_2x2.fits'},
    #                        mask files as a function of x,y binning

    # registration settings (Scamp)
    'scamp-config-file': rootpath+'/setup/soargoodman.scamp',
    'reg_max_mag': 19,
    'reg_search_radius': 0.5,  # deg
    'source_tolerance': 'high',

    # swarp settings
    'copy_keywords': ('TELESCOP,INSTRUME,FILTER,EXPTIME,OBJECT,' +
                      'DATE-OBS,RA,DEC,SECPIX,AIRMASS,' +
                      'TEL_KEYW'),
    #                        keywords to be copied in image
    #                        combination using swarp
    'swarp-config-file': rootpath+'/setup/soargoodman.swarp',

    # default catalog settings
    'astrometry_catalogs': ['GAIA'],
    'photometry_catalogs': ['PANSTARRS', 'SkyMapper', 'SDSS-R9']
}

# SOAR, Goodman [new] (imager)
# running Goodman data requires the removal of header keywords
# PARAM0, PARAM61, PARAM62, PARAM63 (degree symbol is non-ASCII)
soargoodman_param = {
    'telescope_instrument': 'SOAR/GOODMAN',  # telescope/instrument name
    'telescope_keyword': 'SOARGOODMAN',  # telescope/instrument keyword
    'observatory_code': 'I33',         # MPC observatory code
    'secpix': (0.15, 0.15),  # pixel size (arcsec)
    # before binning
    'ext_coeff': 0.05,          # typical extinction coefficient


    # image orientation preferences
    'flipx': True,
    'flipy': True,
    'rotate': 0,

    # instrument-specific FITS header keywords
    'binning': ('PG5_9', 'PG5_4'),  # binning in x/y
    'extent': ('NAXIS1', 'NAXIS2'),   # N_pixels in x/y
    'ra': 'RA',  # telescope pointing, RA
    'dec': 'DEC',  # telescope pointin, Dec
    'radec_separator': ':',   # RA/Dec hms separator, use 'XXX'
    # if already in degrees
    'date_keyword': 'DATE-OBS',  # obs date/time
    # keyword; use
    # 'date|time' if
    # separate
    'obsmidtime_jd': 'MIDTIMJD',  # obs midtime jd keyword
    # (usually provided by
    # pp_prepare
    'object': 'OBJECT',  # object name keyword
    'filter': 'FILTER',  # filter keyword
    'filter_translations': {'Rc': 'R', 'V': 'V', 'B': 'B', 'u': 'u',
                              'B-Bessel': 'B', 'V-Bessel': 'V',
                              'R-Bessel': 'R', 'I-Bessel': 'I',
                              'g-SDSS': 'g', 'r-SDSS': 'r', 'i-SDSS': 'i',
                              'z-SDSS': 'z', 'VR': None,
                              '<NO FILTER>': None},
    # filtername translation dictionary
    'exptime': 'EXPTIME',  # exposure time keyword (s)
    'airmass': 'AIRMASS',  # airmass keyword


    # source extractor settings
    'source_minarea': 15,  # default sextractor source minimum N_pixels
    'source_snr': 2,  # default sextractor source snr for registration
    'aprad_default': 4,  # default aperture radius in px
    'aprad_range': [2, 10],  # [minimum, maximum] aperture radius (px)
    'sex-config-file': rootpath+'/setup/soargoodman.sex',
    'mask_file': {'1,1': rootpath+'/setup/mask_soargoodman_1x1.fits',
                  '2,2': rootpath+'/setup/mask_soargoodman_2x2.fits'},
    #                        mask files as a function of x,y binning

    # registration settings (Scamp)
    'scamp-config-file': rootpath+'/setup/soargoodman.scamp',
    'reg_max_mag': 21,
    'reg_search_radius': 0.5,  # deg
    'source_tolerance': 'high',

    # swarp settings
    'copy_keywords': ('TELESCOP,INSTRUME,FILTER,EXPTIME,OBJECT,' +
                      'DATE-OBS,RA,DEC,SECPIX,AIRMASS,PG5_9,PG5_4,' +
                      'MIDTIMJD,TEL_KEYW'),
    #                        keywords to be copied in image
    #                        combination using swarp
    'swarp-config-file': rootpath+'/setup/soargoodman.swarp',

    # default catalog settings
    'astrometry_catalogs': ['GAIA'],
    'photometry_catalogs': ['PANSTARRS', 'SkyMapper', 'SDSS-R9']
}


# Observatoire de Haute-Provence, 120cm, CCD
ohp120_param = {
    'telescope_instrument': 'OHP120cm/CCD',  # telescope/instrument name
    'telescope_keyword': 'OHP120',      # telescope/instrument keyword
    'observatory_code': '511',         # MPC observatory code
    'secpix': (0.385, 0.385),  # pixel size (arcsec) before binning

    # image orientation preferences
    'flipx': True,
    'flipy': True,
    'rotate': 0,

    # instrument-specific FITS header keywords
    'binning': ('XBINNING', 'YBINNING'),  # binning in x/y
    'extent': ('NAXIS1', 'NAXIS2'),   # N_pixels in x/y
    'ra': 'RA',  # telescope pointing, RA
    'dec': 'DEC',  # telescope pointing, Dec
    'radec_separator': 'XXX',   # RA/Dec hms separator, use 'XXX'
    # if already in degrees
    'date_keyword': 'DATE-OBS',  # obs date/time
    # keyword; use
    # 'date|time' if
    # separate
    'obsmidtime_jd': 'MIDTIMJD',  # obs midtime jd keyword
    # (usually provided by
    # pp_prepare
    'object': 'OBJECT',  # object name keyword
    'filter': 'FILTER',  # filter keyword
    'filter_translations': {'R_Cousins': 'R', 'V_Cousins': 'V',
                            'B_Cousins': 'B', 'H-alpha': None,
                            'i_Gunn': 'i'},
    # filtername translation dictionary

    'exptime': 'EXPTIME',  # exposure time keyword (s)
    'airmass': 'AIRMASS',  # airmass keyword

    # source extractor settings
    'source_minarea': 15,  # default sextractor source minimum N_pixels
    'source_snr': 3,  # default sextractor source snr for registration
    'aprad_default': 8,  # default aperture radius in px
    'aprad_range': [2, 15],  # [minimum, maximum] aperture radius (px)
    'sex-config-file': rootpath+'/setup/ohp120.sex',
    'mask_file': {},
    #                        mask files as a function of x,y binning

    # registration settings (Scamp)
    'scamp-config-file': rootpath+'/setup/ohp120.scamp',
    'reg_max_mag': 19,
    'reg_search_radius': 0.5,  # deg
    'source_tolerance': 'high',

    # default catalog settings
    'astrometry_catalogs': ['GAIA'],
    'photometry_catalogs': ['SDSS-R9', 'APASS9', 'PANSTARRS', '2MASS']
}


# Telescopio Nazionale Galileo, DOLORES
tngdolores_param = {
    'telescope_instrument': 'TNG/DOLORES',  # telescope/instrument name
    'telescope_keyword': 'TNGDOLORES',      # telescope/instrument keyword
    'observatory_code': 'Z19',         # MPC observatory code
    'secpix': (0.252, 0.252),  # pixel size (arcsec) before binning

    # image orientation preferences
    'flipx': True,
    'flipy': True,
    'rotate': 0,

    # instrument-specific FITS header keywords
    'binning': ('CRDELT1', 'CRDELT2'),  # binning in x/y
    'extent': ('NAXIS1', 'NAXIS2'),   # N_pixels in x/y
    'ra': 'RA',  # telescope pointing, RA
    'dec': 'DEC',  # telescope pointing, Dec
    'radec_separator': 'XXX',   # RA/Dec hms separator, use 'XXX'
    # if already in degrees
    'date_keyword': 'DATE-OBS',  # obs date/time
    # keyword; use
    # 'date|time' if
    # separate
    'obsmidtime_jd': 'MIDTIMJD',  # obs midtime jd keyword
    # (usually provided by
    # pp_prepare
    'object': 'OBJCAT',  # object name keyword
    'filter': 'FLT_ID',  # filter keyword
    'filter_translations': {'B_JOHN_10': 'B', 'V_JOHN_11': 'V',
                            'R_JOHN_12': 'R', 'I_JOHN_13': 'I'},
    # filtername translation dictionary
    'exptime': 'EXPTIME',  # exposure time keyword (s)
    'airmass': 'AIRMASS',  # airmass keyword

    # source extractor settings
    'source_minarea': 15,  # default sextractor source minimum N_pixels
    'source_snr': 3,  # default sextractor source snr for registration
    'aprad_default': 8,  # default aperture radius in px
    'aprad_range': [2, 15],  # [minimum, maximum] aperture radius (px)
    'sex-config-file': rootpath+'/setup/tngdolores.sex',
    'mask_file': {},
    #                        mask files as a function of x,y binning

    # registration settings (Scamp)
    'scamp-config-file': rootpath+'/setup/tngdolores.scamp',
    'reg_max_mag': 17,
    'reg_search_radius': 0.5,  # deg
    'source_tolerance': 'high',

    # default catalog settings
    'astrometry_catalogs': ['GAIA'],
    'photometry_catalogs': ['SDSS-R9', 'APASS9', 'PANSTARRS', '2MASS']
}

# KPNO 4m Mayall, MOSAIC-1, 1.1 and 3
kpno4mos1_param = {
    'telescope_instrument': 'KPNO4m/MOSAIC',  # telescope/instrument name
    'telescope_keyword': 'KPNO4MOS1',  # telescope/instrument keyword
    'observatory_code': '695',  # MPC observatory code
    'secpix': (0.2666, 0.2666),  # pixel size (arcsec) before binning
    #'secpix': (0.25, 0.25),  # pixel size (arcsec) before binning

    # image orientation preferences
    'flipx': True,
    'flipy': False,
    'rotate': 0,

    # instrument-specific FITS header keywords
    'binning': (1, 1),  # binning in x/y
    'extent': ('NAXIS1', 'NAXIS2'),  # N_pixels in x/y
    'ra': 'RA',  # telescope pointing, RA
    'dec': 'DEC',  # telescope pointin, Dec
    'radec_separator': ':',  # RA/Dec hms separator, use 'XXX'
    # if already in degrees
    'date_keyword': 'DATE-OBS',  # obs date/time
    # keyword; use
    # 'date|time' if
    # separate
    'obsmidtime_jd': 'MJD-OBS',  # obs midtime jd keyword
    # (usually provided by
    # pp_prepare
    'object': 'OBJECT',  # object name keyword
    'filter': 'FILTER',  # filter keyword
    'filter_translations': {'g SDSS k1017': 'g',
                            'r SDSS k1018': 'r',
                            'i SDSS k1019': 'i',
                            'z SDSS c6020': 'z',  # mosaic1.1/3
                            'z SDSS k1020': 'z',  # broken filter
                            'KXs': 'K',
                            'none': None},
    # filtername translation dictionary
    'exptime': 'EXPTIME',  # exposure time keyword (s)
    'airmass': 'AIRMASS',  # airmass keyword

    # source extractor settings
    'source_minarea': 5,  # default sextractor source minimum N_pixels
    'source_snr': 1.5,  # default sextractor snr for registration
    'aprad_default': 5,  # default aperture radius in px
    'aprad_range': [2, 10],  # [minimum, maximum] aperture radius (px)
    'sex-config-file': rootpath + '/setup/kpno4mos1.sex',
    'mask_file': {},
    #                        mask files as a function of x,y binning

    # scamp settings
    'scamp-config-file': rootpath + '/setup/kpno4mos1.scamp',
    'reg_max_mag': 19,
    'reg_search_radius': 0.5,  # deg
    'source_tolerance': 'low',

    # default catalog settings
    'astrometry_catalogs': ['GAIA'],
    'photometry_catalogs': ['SDSS-R13', 'PANSTARRS', 'APASS9', '2MASS']
}

# KPNO 4m Mayall, MOSAIC3
kpno4mos3_param = {
    'telescope_instrument': 'KPNO4m/MOSAIC',  # telescope/instrument name
    'telescope_keyword': 'KPNO4MOS3',  # telescope/instrument keyword
    'observatory_code': '695',  # MPC observatory code
    'secpix': (0.25, 0.25),  # pixel size (arcsec) before binning
    #'secpix': (0.2666, 0.2666),  # pixel size (arcsec) before binning

    # image orientation preferences
    'flipx': True,
    'flipy': False,
    'rotate': 0,

    # instrument-specific FITS header keywords
    'binning': (1, 1),  # binning in x/y
    'extent': ('NAXIS1', 'NAXIS2'),  # N_pixels in x/y
    'ra': 'RA',  # telescope pointing, RA
    'dec': 'DEC',  # telescope pointin, Dec
    'radec_separator': ':',  # RA/Dec hms separator, use 'XXX'
    # if already in degrees
    'date_keyword': 'DATE-OBS',  # obs date/time
    # keyword; use
    # 'date|time' if
    # separate
    'obsmidtime_jd': 'MJD-OBS',  # obs midtime jd keyword
    # (usually provided by
    # pp_prepare
    'object': 'OBJECT',  # object name keyword
    'filter': 'FILTER',  # filter keyword
    'filter_translations': {'g SDSS k1017': 'g',
                            'r SDSS k1018': 'r',
                            'i SDSS k1019': 'i',
                            'z SDSS c6020': 'z',
                            'B Harris k1002': 'B',
                            'KXs': 'K',
                            'none': None},
    # filtername translation dictionary
    'exptime': 'EXPTIME',  # exposure time keyword (s)
    'airmass': 'AIRMASS',  # airmass keyword

    # source extractor settings
    'source_minarea': 9,  # default sextractor source minimum N_pixels
    'source_snr': 3,  # default sextractor snr for registration
    'aprad_default': 5,  # default aperture radius in px
    'aprad_range': [3, 15],  # [minimum, maximum] aperture radius (px)
    'sex-config-file': rootpath + '/setup/kpno4mos1.sex',
    'mask_file': {},
    #                        mask files as a function of x,y binning

    # scamp settings
    'scamp-config-file': rootpath + '/setup/kpno4mos1.scamp',
    'reg_max_mag': 21,
    'reg_search_radius': 0.5,  # deg
    'source_tolerance': 'low',

    # default catalog settings
    'astrometry_catalogs': ['GAIA'],
    'photometry_catalogs': ['SDSS-R13', 'PANSTARRS', 'APASS9', '2MASS']
}

# KPNO 4m Mayall, NEWFIRM
kpno4newf_param = {
    'telescope_instrument': 'KPNO4m/NEWFIRM',  # telescope/instrument name
    'telescope_keyword': 'KPNO4NEWF',  # telescope/instrument keyword
    'observatory_code': '695',  # MPC observatory code
    'secpix': (0.4, 0.4),  # pixel size (arcsec) before binning
    #'secpix': (0.2666, 0.2666),  # pixel size (arcsec) before binning

    # image orientation preferences
    'flipx': True,
    'flipy': False,
    'rotate': 0,

    # instrument-specific FITS header keywords
    'binning': (1, 1),  # binning in x/y
    'extent': ('NAXIS1', 'NAXIS2'),  # N_pixels in x/y
    'ra': 'RA',  # telescope pointing, RA
    'dec': 'DEC',  # telescope pointin, Dec
    'radec_separator': ':',  # RA/Dec hms separator, use 'XXX'
    # if already in degrees
    'date_keyword': 'DATE-OBS',  # obs date/time
    # keyword; use
    # 'date|time' if
    # separate
    'obsmidtime_jd': 'MJD-OBS',  # obs midtime jd keyword
    # (usually provided by
    # pp_prepare
    'object': 'OBJECT',  # object name keyword
    'filter': 'FILTER',  # filter keyword
    'filter_translations': {'KXs': 'K',
                            'JX': 'J',
                            'none': None},
    # filtername translation dictionary
    'exptime': 'EXPTIME',  # exposure time keyword (s)
    'airmass': 'AIRMASS',  # airmass keyword

    # source extractor settings
    'source_minarea': 5,  # default sextractor source minimum N_pixels
    'source_snr': 1.5,  # default sextractor snr for registration
    'aprad_default': 10,  # default aperture radius in px
    'aprad_range': [5, 20],  # [minimum, maximum] aperture radius (px)
    'sex-config-file': rootpath + '/setup/kpno4mos1.sex',
    'mask_file': {},
    #                        mask files as a function of x,y binning

    # scamp settings
    'scamp-config-file': rootpath + '/setup/kpno4mos1.scamp',
    'reg_max_mag': 19,
    'reg_search_radius': 0.5,  # deg
    'source_tolerance': 'low',

    # default catalog settings
    'astrometry_catalogs': ['GAIA'],
    'photometry_catalogs': ['2MASS']
}

# KMTNET-S
kmtnets_param = {
    'telescope_instrument': 'KMTNET-S',  # telescope/instrument name
    'telescope_keyword': 'KMTNETS',  # telescope/instrument keyword
    'observatory_code': 'K94',  # MPC observatory code
    'secpix': (0.398, 0.398),  # pixel size (arcsec) before binning

    # image orientation preferences
    'flipx': True,
    'flipy': False,
    'rotate': 0,

    # instrument-specific FITS header keywords
    'binning': ('CCDXBIN', 'CCDYBIN'),  # binning in x/y
    'extent': ('NAXIS1', 'NAXIS2'),  # N_pixels in x/y
    # use for crop fields
    # 'ra': 'CRVAL1',  # telescope pointing, RA
    # 'dec': 'CRVAL2',  # telescope pointin, Dec
    # use for single CCDs
    'ra': 'CCD_RA',  # telescope pointing, RA
    'dec': 'CCD_DEC',  # telescope pointin, Dec
    'radec_separator': 'XXX',  # RA/Dec hms separator, use 'XXX'
    # use full CCD mosaic
    # 'ra': 'RA',  # telescope pointing, RA
    # 'dec': 'DEC',  # telescope pointin, Dec
    # 'radec_separator': ':',  # RA/Dec hms separator, use 'XXX'

    # default distortion parameters per CCD
    # using individual CCD centers for CRPIX1/2, CRVAL1/2
    # using 1x1 binning
    'distort': {'functionof': 'CCD_NAME',
                'K': {'PV1_0': -0.000711762135419,
                      'PV1_1': 1.00171981124,
                      'PV1_2': -0.000228317298275,
                      'PV1_4': 0.0160316513101,
                      'PV1_5': -0.0113636730644,
                      'PV1_6': 0.00549290221274,
                      'PV1_7': -0.010384662223,
                      'PV1_8': 0.00081137447258,
                      'PV1_9': -0.0104005033165,
                      'PV1_10': 0.000396116525231,
                      'PV2_0': 0.00319659825814,
                      'PV2_1': 1.00057217492,
                      'PV2_2': 0.0116970617468,
                      'PV2_4': -0.0165165992945,
                      'PV2_5': 0.0108419617955,
                      'PV2_6': -0.00548712635243,
                      'PV2_7': -0.0103126090939,
                      'PV2_8': 0.000686169735533,
                      'PV2_9': -0.0103739930263,
                      'PV2_10': 0.000240139308823},
                'M': {'PV1_0': 0.000460033317474,
                      'PV1_1': 1.00141889978,
                      'PV1_2': -0.000430380594516,
                      'PV1_4': -0.015400655054,
                      'PV1_5': -0.0115995667827,
                      'PV1_6': -0.00518535937805,
                      'PV1_7': -0.0101118044677,
                      'PV1_8': -3.19250493138e-05,
                      'PV1_9': -0.0106737708283,
                      'PV1_10': -0.000431356736006,
                      'PV2_0': 0.00198339122013,
                      'PV2_1': 0.999670425747,
                      'PV2_2': -0.011005193782,
                      'PV2_4': -0.0167779694087,
                      'PV2_5': -0.0106335253045,
                      'PV2_6': -0.00526313446543,
                      'PV2_7': -0.0101955642118,
                      'PV2_8': -0.000255088245494,
                      'PV2_9': -0.0094035107269,
                      'PV2_10': -0.000292075883415},
                'T': {'PV1_0': -0.00127421419519,
                      'PV1_1': 1.00104160823,
                      'PV1_2': -0.000660886473555,
                      'PV1_4': 0.0158990914667,
                      'PV1_5': 0.01169760742,
                      'PV1_6': 0.00543678807381,
                      'PV1_7': -0.0103261215423,
                      'PV1_8': -0.000794914727406,
                      'PV1_9': -0.0103649052751,
                      'PV1_10': -0.000279241301327,
                      'PV2_0': -0.00392674586144,
                      'PV2_1': 0.999343102486,
                      'PV2_2': -0.0111411751205,
                      'PV2_4': 0.017084775899,
                      'PV2_5': 0.010790771213,
                      'PV2_6': 0.00566154555136,
                      'PV2_7': -0.0102149292801,
                      'PV2_8': -0.000623538149787,
                      'PV2_9': -0.0102808588946,
                      'PV2_10': -0.000162027267646},
                'N': {'PV1_0':  0.000594262478073,
                      'PV1_1': 0.998957324393,
                      'PV1_2': 0.00141244617639,
                      'PV1_4': -0.015924297871,
                      'PV1_5': 0.0114088931271,
                      'PV1_6': -0.00544474906971,
                      'PV1_7': -0.010256940745,
                      'PV1_8': 0.000789499342505,
                      'PV1_9': -0.0101626175866,
                      'PV1_10': 0.000179932486564,
                      'PV2_0': -0.00166315601111,
                      'PV2_1': 0.997680853652,
                      'PV2_2': 0.0102233336556,
                      'PV2_4': 0.0166335481995,
                      'PV2_5': -0.0107970324189,
                      'PV2_6': 0.00556273542573,
                      'PV2_7': -0.0103304464752,
                      'PV2_8': 0.00070600868855,
                      'PV2_9': -0.0102590932248,
                      'PV2_10': 0.000205537467251}},


    # if already in degrees
    'date_keyword': 'DATE-OBS',  # obs date/time
    # keyword; use
    # 'date|time' if
    # separate
    'obsmidtime_jd': 'MJD-OBS',  # obs midtime jd keyword
    # (usually provided by
    # pp_prepare
    'object': 'OBJECT',  # object name keyword
    'filter': 'FILTER',  # filter keyword
    'filter_translations': {'V': 'V',
                            'R': 'R',
                            'I': 'I'},
    # filtername translation dictionary
    'exptime': 'EXPTIME',  # exposure time keyword (s)
    'airmass': 'SECZ',  # airmass keyword

    # source extractor settings
    'source_minarea': 10,  # default sextractor source minimum N_pixels
    'source_snr': 20,  # default sextractor source snr for registration
    'aprad_default': 5,  # default aperture radius in px
    'aprad_range': [2, 10],  # [minimum, maximum] aperture radius (px)
    'sex-config-file': rootpath + '/setup/kmtnets.sex',
    #'mask_file': {'1,1' : rootpath+'/setup/mask_kmtnets_1x1.fits'},
    'mask_file': {},
    #                        mask files as a function of x,y binning

    # registration settings (Scamp)
    'scamp-config-file': rootpath + '/setup/kmtnets.scamp',
    'reg_max_mag': 18,
    'reg_search_radius': 0.5,  # deg
    'source_tolerance': 'none',
    'scamp': {'ASTRINSTRU_KEY': 'FILTER,CCD_NAME'},

    # default catalog settings
    'astrometry_catalogs': ['GAIA', 'GAIA'],  # run registration twice
    # due to large field distortions
    'photometry_catalogs': ['SDSS-R9', 'APASS9', 'PANSTARRS', '2MASS'],

    # list of header keywords that should not be removed
    'dont_remove': 'CCD_NAME'
}

# Flagstaff Robotic survey telescopes
frost_param = {
    'telescope_instrument': 'FRoST',  # telescope/instrument name
    'telescope_keyword': 'FROST',      # telescope/instrument keyword
    'observatory_code': 'V04',         # MPC observatory code
    'secpix': (2.81, 2.81),  # pixel size (arcsec)
    # before binning
    'ext_coeff': 0.05,          # typical extinction coefficient


    # image orientation preferences
    'flipx': True,
    'flipy': False,
    'rotate': 0,

    # instrument-specific FITS header keywords
    'binning': (1, 1),  # binning in x/y
    'extent': ('NAXIS1', 'NAXIS2'),   # N_pixels in x/y
    'ra': 'RA',  # telescope pointing, RA
    'dec': 'DEC',  # telescope pointin, Dec
    'radec_separator': 'XXX',   # RA/Dec hms separator, use 'XXX'
    # if already in degrees
    'date_keyword': 'DATE-OBS',  # obs date/time
    # keyword; use
    # 'date|time' if
    # separate
    'obsmidtime_jd': 'MIDTIMJD',  # obs midtime jd keyword
    # (usually provided by
    # pp_prepare
    'object': 'OBJECT',  # object name keyword
    'filter': 'FILTER',  # filter keyword
    'filter_translations': {'clear': 'V', 'Clear': 'V', 'CLEAR': 'V'},
    # filtername translation dictionary
    'exptime': 'EXPTIME',  # exposure time keyword (s)
    'airmass': 'AIRMASS',  # airmass keyword


    # source extractor settings
    'source_minarea': 5,  # default sextractor source minimum N_pixels
    'source_snr': 3,  # default sextractor source snr for registration
    'aprad_default': 3,  # default aperture radius in px
    'aprad_range': [1, 8],  # [minimum, maximum] aperture radius (px)
    'sex-config-file': rootpath+'/setup/frost.sex',
    'mask_file': {},
    #                        mask files as a function of x,y binning

    # registration settings (Scamp)
    'scamp-config-file': rootpath+'/setup/frost.scamp',
    'reg_max_mag': 17,
    'reg_search_radius': 0.5,  # deg
    'source_tolerance': 'high',

    # swarp settings
    'copy_keywords': ('TELESCOP,INSTRUME,FILTER,EXPTIME,OBJECT,' +
                      'DATE-OBS,TIME-OBS,RA,DEC,SECPIX,AIRMASS,' +
                      'TEL_KEYW'),
    #                         keywords to be copied in image
    #                         combination using swarp
    'swarp-config-file': rootpath+'/setup/frost.swarp',

    # default catalog settings
    'astrometry_catalogs': ['GAIA'],
    'photometry_catalogs': ['SDSS-R9', 'APASS9', 'PANSTARRS', '2MASS']
}


# SPM 84cm, Mexman
mexman_param = {
    'telescope_instrument': '0.84cm/MEXMAN',  # telescope/instrument name
    'telescope_keyword': 'MEXMAN',  # telescope/instrument keyword
    'observatory_code': '679',  # MPC observatory code
    'secpix': (0.22, 0.22),  # pixel size (arcsec) before binning


    # image orientation preferences
    'flipx': False,
    'flipy': True,
    'rotate': 90,


    # instrument-specific FITS header keywords
    'binning': ('CCDXBIN', 'CCDYBIN'),  # binning in x/y
    'extent': ('NAXIS1', 'NAXIS2'),  # N_pixels in x/y
    'ra': 'RA',  # telescope pointing, RA
    'dec': 'DEC',  # telescope pointin, Dec
    'radec_separator': ':',  # RA/Dec hms separator, use 'XXX'
    # if already in degrees
    'date_keyword': 'DATE-OBS|UT',  # obs date/time
    # keyword; use
    # 'date|time' if
    # separate
    'obsmidtime_jd': 'JD',  # obs midtime jd keyword
    # (usually provided by
    # pp_prepare
    'object': 'OBJECT',  # object name keyword
    'filter': 'FILTER',  # filter keyword
    'filter_translations': {'U': 'U', 'B': 'B', 'V': 'V', 'R': 'R', 'I': 'I'},
    # filtername translation dictionary
    'exptime': 'EXPTIME',  # exposure time keyword (s)
    'airmass': 'AIRMASS',  # airmass keyword


    # source extractor settings
    'source_minarea': 12,  # default sextractor source minimum N_pixels
    'source_snr': 3,  # default sextractor source snr for registration
    'aprad_default': 5,  # default aperture radius in px
    'aprad_range': [2, 10],  # [minimum, maximum] aperture radius (px)
    'sex-config-file': rootpath + '/setup/mexman.sex',
    'mask_file': {},
    #                        mask files as a function of x,y binning


    # scamp settings
    'scamp-config-file': rootpath + '/setup/mexman.scamp',
    'reg_max_mag': 19,
    'reg_search_radius': 0.5,  # deg
    'source_tolerance': 'high',


    # swarp settings
    'copy_keywords': ('TELESCOP,INSTRUME,FILTER,EXPTIME,OBJECT,' +
                      'DATE-OBS,RA,DEC,AIRMASS'),
    #                         keywords to be copied in image
    #                         combination using swarp
    'swarp-config-file': rootpath+'/setup/mexman.swarp',


    # default catalog settings
    'astrometry_catalogs': ['GAIA'],
    'photometry_catalogs': ['SDSS-R9', 'APASS9', '2MASS']
}

# UKIRT, WFCAM
ukirtwfcam_param = {
    'telescope_instrument': 'UKIRT/WFCAM',  # telescope/instrument name
    'telescope_keyword': 'UKIRTWFCAM',      # telescope/instrument keyword
    'observatory_code': '568',         # MPC observatory code
    'secpix': (0.402, 0.402),  # pixel size (arcsec)
    # before binning
    'ext_coeff': 0.05,          # typical extinction coefficient


    # # image orientation preferences
    # 'flipx'                : False,
    # 'flipy'                : False,
    # 'rotate'               : 0,

    # image orientation preferences (for each chip)
    'chip_id': 'CAMNUM',        # chip identifier (remove,
    # if not existent)
    # the following keys are dictionaries if 'chip_id' exists, single
    # values otherwise
    'flipx': {1: True, 2: True, 3: True, 4: True},
    'flipy': {1: True, 2: True, 3: True, 4: True},
    'rotate': {1: 270, 2: 0, 3: 90, 4: 180},
    'chip_offset_fixed': {1: (-0.44, -0.44), 2: (0, -0.44),
                          3: (0, 0), 4: (-0.440, 0)},
    # chip offset (ra, dec in degress) [optional]


    # instrument-specific FITS header keywords
    'binning': (1, 1),  # binning in x/y
    'extent': ('NAXIS1', 'NAXIS2'),   # N_pixels in x/y
    'ra': 'TELRA',  # telescope pointing, RA
    'dec': 'TELDEC',  # telescope pointin, Dec
    'radec_separator': 'XXX',   # RA/Dec hms separator, use 'XXX'
    # if already in degrees
    'date_keyword': 'DATE-OBS',  # obs date/time
    # keyword; use
    # 'date|time' if
    # separate
    'obsmidtime_jd': 'MIDTIMJD',  # obs midtime jd keyword
    # (usually provided by
    # pp_prepare
    'object': 'OBJECT',  # object name keyword
    'filter': 'FILTER',  # filter keyword
    'filter_translations': {'J': 'J', 'Z': 'Z',
                            'H': 'H', 'K': 'K'},
    # filtername translation dictionary
    'exptime': 'EXP_TIME',  # exposure time keyword (s)
    'airmass': 'AMSTART',  # airmass keyword

    # source extractor settings
    'source_minarea': 12,  # default sextractor source minimum N_pixels
    'source_snr': 3,  # default sextractor source snr for registration
    'aprad_default': 5,  # default aperture radius in px
    'aprad_range': [1, 8],  # [minimum, maximum] aperture radius (px)
    'sex-config-file': rootpath+'/setup/ukirtwfcam.sex',
    'mask_file': {},
    #                        mask files as a function of x,y binning

    # registration settings (Scamp)
    'scamp-config-file': rootpath+'/setup/ukirtwfcam.scamp',
    'reg_max_mag': 19,
    'reg_search_radius': 0.05,  # deg
    'source_tolerance': 'high',

    # swarp settings
    'copy_keywords': ('TELESCOP,INSTRUME,FILTER,EXP_TIME,PROJECT,' +
                      'AIRMASS,OBJECT,DATE-OBS,MJD-OBS,GAIN,' +
                      'READNOIS,TELRA,TELDEC,SECPIX1,SECPIX2,' +
                      'CAMNUM,MIDTIMJD,TEL_KEYW,AMSTART'),
    #                         keywords to be copied in image
    #                         combination using swarp
    'swarp-config-file': rootpath+'/setup/ukirtwfcam.swarp',

    # default catalog settings
    'astrometry_catalogs': ['GAIA'],
    'photometry_catalogs': ['2MASS']
}

# IRSF 1.4m, SIRIUS
irsfsirius_param = {
    'telescope_instrument': 'IRSF/SIRIUS',  # telescope/instrument name
    'telescope_keyword': 'IRSFSIRIUS',      # telescope/instrument keyword
    'observatory_code': 'K94',         # MPC observatory code
    'secpix': (0.453, 0.453),  # pixel size (arcsec)
    # before binning
    'ext_coeff': 0.05,          # typical extinction coefficient


    # image orientation preferences
    'flipx': True,
    'flipy': False,
    'rotate': 0,

    # instrument-specific FITS header keywords
    'binning': (1, 1),  # binning in x/y
    'extent': ('NAXIS1', 'NAXIS2'),   # N_pixels in x/y
    'ra': 'RA',  # telescope pointing, RA
    'dec': 'DEC',  # telescope pointin, Dec
    'radec_separator': ':',   # RA/Dec hms separator, use 'XXX'
    # if already in degrees
    'date_keyword': 'DATE_UTC|TIME_UTC',  # obs date/time
    # keyword; use
    # 'date|time' if
    # separate
    'obsmidtime_jd': 'MIDTIMJD',  # obs midtime jd keyword
    # (usually provided by
    # pp_prepare
    'object': 'OBJECT',  # object name keyword
    'filter': 'FILTER',  # filter keyword
    'filter_translations': {'J': 'J', 'H': 'H', 'Ks': 'Ks'},
    # filtername translation dictionary
    'exptime': 'EXPOS',  # exposure time keyword (s)
    'airmass': 'AIRMASS',  # airmass keyword


    # source extractor settings
    'source_minarea': 12,  # default sextractor source minimum N_pixels
    'source_snr': 3,  # default sextractor source snr for registration
    'aprad_default': 5,  # default aperture radius in px
    'aprad_range': [2, 10],  # [minimum, maximum] aperture radius (px)
    'sex-config-file': rootpath+'/setup/irsfsirius.sex',
    'mask_file': {},
    #                        mask files as a function of x,y binning

    # registration settings (Scamp)
    'scamp-config-file': rootpath+'/setup/irsfsirius.scamp',
    'reg_max_mag': 19,
    'reg_search_radius': 0.5,  # deg
    'source_tolerance': 'high',

    # swarp settings
    'copy_keywords': ('TELESCOP,INSTRUME,FILTER,EXPTIME,OBJECT,' +
                      'DATE-OBS,TIME-OBS,RA,DEC,SECPIX,AIRMASS,' +
                      'TEL_KEYW'),
    #                         keywords to be copied in image
    #                         combination using swarp
    'swarp-config-file': rootpath+'/setup/vatt4k.swarp',

    # default catalog settings
    'astrometry_catalogs': ['GAIA'],
    'photometry_catalogs': ['2MASS']
}

# VLT, FORS2
vltfors2_param = {
    'telescope_instrument': 'VLT/FORS2',  # telescope/instrument name
    'telescope_keyword': 'VLTFORS2',  # telescope/instrument keyword
    'observatory_code': '309',         # MPC observatory code
    'secpix': (0.126, 0.126),  # pixel size (arcsec)
    # before binning
    'ext_coeff': 0.05,          # typical extinction coefficient


    # image orientation preferences
    'flipx': True,
    'flipy': False,
    'rotate': 0,

    # instrument-specific FITS header keywords
    'binning': ('ESO DET WIN1 BINX', 'ESO DET WIN1 BINY'),
    # binning in x/y, '_blankN' denotes that both axes
    # are listed in one keyword, sep. by blanks
    'extent': ('NAXIS1', 'NAXIS2'),   # N_pixels in x/y
    'ra': 'RA',  # telescope pointing, RA
    'dec': 'DEC',  # telescope pointin, Dec
    'radec_separator': 'XXX',   # RA/Dec hms separator, use 'XXX'
    # if already in degrees
    'date_keyword': 'DATE-OBS',  # obs date/time
    # keyword; use
    # 'date|time' if
    # separate
    'obsmidtime_jd': 'MIDTIMJD',  # obs midtime jd keyword
    # (usually provided by
    # pp_prepare
    'object': 'OBJECT',  # object name keyword
    'filter': 'ESO INS FILT1 NAME',  # filter keyword
    'filter_translations': {'R_SPECIAL': 'r'},
    # filtername translation dictionary
    'exptime': 'EXPTIME',  # exposure time keyword (s)
    'airmass': 'AIRMASS',  # airmass keyword


    # source extractor settings
    'source_minarea': 10,  # default sextractor source minimum N_pixels
    'source_snr': 5,  # default sextractor source snr for registration
    'aprad_default': 4,  # default aperture radius in px
    'aprad_range': [2, 20],  # [minimum, maximum] aperture radius (px)
    'sex-config-file': rootpath+'/setup/vltfors2.sex',
    'mask_file': {},
    #                        mask files as a function of x,y binning

    # registration settings (Scamp)
    'scamp-config-file': rootpath+'/setup/vltfors2.scamp',
    'reg_max_mag': 22,
    'reg_search_radius': 0.1,  # deg
    'source_tolerance': 'high',

    # swarp settings
    # swarp does not work for FORS2 due to hierarchical header keywords
    # requires a workaround
    'copy_keywords': ('OBSERVAT,INSTRUME,'
                      'EXPTIME,'
                      'OBJECT,DATE-OBS,RA,DEC,SCALE,AIRMASS,'
                      'TEL_KEYW'),
    #                        keywords to be copied in image
    #                        combination using swarp
    'swarp-config-file': rootpath+'/setup/vltfors2.swarp',

    # default catalog settings
    'astrometry_catalogs': ['GAIA'],
    'photometry_catalogs': ['SDSS-R9', 'PANSTARRS', 'APASS9']
}

# Pluto plate
plutoplate_param = {
    'telescope_instrument': 'Pluto/plate',  # telescope/instrument name
    'telescope_keyword': 'PLUTOPLATE',      # telescope/instrument keyword
    'observatory_code': '690',         # MPC observatory code
    'secpix': (0.9, 0.9),  # pixel size (arcsec)
    # before binning
    'ext_coeff': 0.05,          # typical extinction coefficient


    # image orientation preferences
    'flipx': False,
    'flipy': False,
    'rotate': 0,

    # instrument-specific FITS header keywords
    'binning': (1, 1),  # binning in x/y
    'extent': ('NAXIS1', 'NAXIS2'),   # N_pixels in x/y
    'ra': 'RA',  # telescope pointing, RA
    'dec': 'DEC',  # telescope pointin, Dec
    'radec_separator': 'XXX',   # RA/Dec hms separator, use 'XXX'
    # if already in degrees
    'date_keyword': 'DATE-OBS',  # obs date/time
    # keyword; use
    # 'date|time' if
    # separate
    'obsmidtime_jd': 'MIDTIMJD',  # obs midtime jd keyword
    # (usually provided by
    # pp_prepare
    'object': 'OBJECT',  # object name keyword
    'filter': 'FILTER',  # filter keyword
    'filter_translations': {'plate': None},
    # filtername translation dictionary
    'exptime': 'EXPTIME',  # exposure time keyword (s)
    'airmass': 'AIRMASS',  # airmass keyword


    # source extractor settings
    'source_minarea': 12,  # default sextractor source minimum N_pixels
    'source_snr': 3,  # default sextractor source snr for registration
    'aprad_default': 5,  # default aperture radius in px
    'aprad_range': [2, 10],  # [minimum, maximum] aperture radius (px)
    'sex-config-file': rootpath+'/setup/plutoplate.sex',
    'mask_file': {},
    #                        mask files as a function of x,y binning

    # registration settings (Scamp)
    'scamp-config-file': rootpath+'/setup/plutoplate.scamp',
    'reg_max_mag': 19,
    'reg_search_radius': 0.5,  # deg
    'source_tolerance': 'high',

    # swarp settings
    'copy_keywords': ('TELESCOP,INSTRUME,FILTER,EXPTIME,OBJECT,' +
                      'DATE-OBS,TIME-OBS,RA,DEC,SECPIX,AIRMASS,' +
                      'TEL_KEYW,CCDBIN1,CCDBIN2,MIDTIMJD'),
    #                         keywords to be copied in image
    #                         combination using swarp
    'swarp-config-file': rootpath+'/setup/vatt4k.swarp',

    # default catalog settings
    'astrometry_catalogs': ['GAIA'],
    'photometry_catalogs': ['PANSTARRS', 'SDSS-R9', 'APASS9']
}

# TCS1.5m/MUSCAT2
tcs15muscat2_param = {
    'telescope_instrument': 'TCS1.5m/MUSCAT2',  # telescope/instrument name
    'telescope_keyword': 'TCS15MUSCAT2',      # telescope/instrument keyword
    'observatory_code': 'J04',         # MPC observatory code
    'secpix': (0.43, 0.43),  # pixel size (arcsec)
    # before binning
    'ext_coeff': 0.05,          # typical extinction coefficient


    # image orientation preferences
    'flipx': True,
    'flipy': False,
    'rotate': 0,

    # instrument-specific FITS header keywords
    'binning': ('BINX', 'BINY'),  # binning in x/y
    'extent': ('NAXIS1', 'NAXIS2'),   # N_pixels in x/y
    'ra': 'RA',  # telescope pointing, RA
    'dec': 'DEC',  # telescope pointin, Dec
    'radec_separator': ':',   # RA/Dec hms separator, use 'XXX'
    # if already in degrees
    'date_keyword': 'DATE-OBS|EXP-STRT',  # obs date/time
    # keyword; use
    # 'date|time' if
    # separate
    'obsmidtime_jd': 'MIDTIMJD',  # obs midtime jd keyword
    # (usually provided by
    # pp_prepare
    'object': 'OBJECT',  # object name keyword
    'filter': 'FILTER',  # filter keyword
    'filter_translations': {'g': 'g', 'r': 'r', 'i': 'i',
                            'z': 'z', 'clear': None},
    # filtername translation dictionary
    'exptime': 'EXPTIME',  # exposure time keyword (s)
    'airmass': 'AIRMASS',  # airmass keyword


    # source extractor settings
    'source_minarea': 12,  # default sextractor source minimum N_pixels
    'source_snr': 3,  # default sextractor source snr for registration
    'aprad_default': 5,  # default aperture radius in px
    'aprad_range': [2, 10],  # [minimum, maximum] aperture radius (px)
    'sex-config-file': rootpath+'/setup/tcs15muscat2.sex',
    'mask_file': {},
    #                        mask files as a function of x,y binning

    # registration settings (Scamp)
    'scamp-config-file': rootpath+'/setup/tcs15muscat2.scamp',
    'reg_max_mag': 19,
    'reg_search_radius': 0.5,  # deg
    'source_tolerance': 'high',

    # swarp settings
    'copy_keywords': ('TELESCOP,INSTRUME,FILTER,EXPTIME,OBJECT,' +
                      'DATE-OBS,EXP-STRT,RA,DEC,SECPIX,AIRMASS,' +
                      'TEL_KEYW,BINX,BINY,MIDTIMJD'),
    #                         keywords to be copied in image
    #                         combination using swarp
    'swarp-config-file': rootpath+'/setup/vatt4k.swarp',

    # default catalog settings
    'astrometry_catalogs': ['GAIA'],
    'photometry_catalogs': ['PANSTARRS', 'SDSS-R9', 'APASS9']
}

# LCOGT, SBIG camera (LSC, KB78)
lcosbigkb78_param = {
    'telescope_instrument': 'LCOGT(LSC)/SBIG',  # telescope/instrument name
    'telescope_keyword': 'LCOSBIGKB78',      # telescope/instrument keyword
    'observatory_code': 'W85',         # MPC observatory code
    'secpix': (0.234, 0.234),  # pixel size (arcsec)
    # before binning
    'ext_coeff': 0.05,          # typical extinction coefficient


    # image orientation preferences
    'flipx': True,
    'flipy': True,
    'rotate': 0,

    # instrument-specific FITS header keywords
    'binning': ('CCDSUM#blank0', 'CCDSUM#blank1'),  # binning in x/y
    'extent': ('NAXIS1', 'NAXIS2'),   # N_pixels in x/y
    'ra': 'RA',  # telescope pointing, RA
    'dec': 'DEC',  # telescope pointin, Dec
    'radec_separator': ':',   # RA/Dec hms separator, use 'XXX'
    # if already in degrees
    'date_keyword': 'DATE-OBS',  # obs date/time
    # keyword; use
    # 'date|time' if
    # separate
    'obsmidtime_jd': 'MIDTIMJD',  # obs midtime jd keyword
    # (usually provided by
    # pp_prepare
    'object': 'OBJECT',  # object name keyword
    'filter': 'FILTER',  # filter keyword
    'filter_translations': {'gp': 'g', 'rp': 'r',
                            'ip': 'i', 'zp': 'z', 'w': None},
    # filtername translation dictionary
    'exptime': 'EXPTIME',  # exposure time keyword (s)
    'airmass': 'AIRMASS',  # airmass keyword


    # source extractor settings
    'source_minarea': 9,  # default sextractor source minimum N_pixels
    'source_snr': 3,  # default sextractor source snr for registration
    'aprad_default': 5,  # default aperture radius in px
    'aprad_range': [2, 10],  # [minimum, maximum] aperture radius (px)
    'sex-config-file': rootpath+'/setup/lcosbig.sex',
    'mask_file': {},
    #                        mask files as a function of x,y binning

    # registration settings (Scamp)
    'scamp-config-file': rootpath+'/setup/lcosbig.scamp',
    'reg_max_mag': 18,
    'reg_search_radius': 0.5,  # deg
    'source_tolerance': 'high',

    # swarp settings
    'copy_keywords': ('TELESCOP,INSTRUME,FILTER,EXPTIME,OBJECT,' +
                      'DATE-OBS,RA,DEC,SECPIX,AIRMASS,' +
                      'TEL_KEYW,MIDTIMJD'),
    #                         keywords to be copied in image
    #                         combination using swarp
    'swarp-config-file': rootpath+'/setup/vatt4k.swarp',

    # default catalog settings
    'astrometry_catalogs': ['GAIA'],
    'photometry_catalogs': ['PANSTARRS', 'SDSS-R9', 'APASS9']
}

# LCOGT, SINISTRO camera (CTIO, FA03)
lcosinfa03_param = {
    # telescope/instrument name
    'telescope_instrument': 'LCOGT(CTIO)/SINISTRO',
    'telescope_keyword': 'LCOSINFA03',      # telescope/instrument keyword
    'observatory_code': 'W86',         # MPC observatory code
    'secpix': (0.389, 0.389),  # pixel size (arcsec)
    # before binning
    'ext_coeff': 0.05,          # typical extinction coefficient


    # image orientation preferences
    'flipx': False,
    'flipy': True,
    'rotate': 0,

    # instrument-specific FITS header keywords
    'binning': ('CCDSUM#blank0', 'CCDSUM#blank1'),  # binning in x/y
    'extent': ('NAXIS1', 'NAXIS2'),   # N_pixels in x/y
    'ra': 'RA',  # telescope pointing, RA
    'dec': 'DEC',  # telescope pointin, Dec
    'radec_separator': ':',   # RA/Dec hms separator, use 'XXX'
    # if already in degrees
    'date_keyword': 'DATE-OBS',  # obs date/time
    # keyword; use
    # 'date|time' if
    # separate
    'obsmidtime_jd': 'MJD-OBS',  # obs midtime jd keyword
    # (usually provided by
    # pp_prepare
    'object': 'OBJECT',  # object name keyword
    'filter': 'FILTER',  # filter keyword
    'filter_translations': {'gp': 'g', 'rp': 'r',
                            'ip': 'i', 'zp': 'z', 'w': None},
    # filtername translation dictionary
    'exptime': 'EXPTIME',  # exposure time keyword (s)
    'airmass': 'AIRMASS',  # airmass keyword


    # source extractor settings
    'source_minarea': 9,  # default sextractor source minimum N_pixels
    'source_snr': 3,  # default sextractor source snr for registration
    'aprad_default': 5,  # default aperture radius in px
    'aprad_range': [2, 10],  # [minimum, maximum] aperture radius (px)
    'sex-config-file': rootpath+'/setup/lcosin.sex',
    'mask_file': {},
    #                        mask files as a function of x,y binning

    # registration settings (Scamp)
    'scamp-config-file': rootpath+'/setup/lcosin.scamp',
    'reg_max_mag': 18,
    'reg_search_radius': 0.5,  # deg
    'source_tolerance': 'high',

    # swarp settings
    'copy_keywords': ('TELESCOP,INSTRUME,FILTER,EXPTIME,OBJECT,' +
                      'DATE-OBS,RA,DEC,SECPIX,AIRMASS,' +
                      'TEL_KEYW,MIDTIMJD'),
    #                         keywords to be copied in image
    #                         combination using swarp
    'swarp-config-file': rootpath+'/setup/vatt4k.swarp',

    # default catalog settings
    'astrometry_catalogs': ['GAIA'],
    'photometry_catalogs': ['PANSTARRS', 'SDSS-R9', 'APASS9']
}


# LCOGT, SINISTRO camera (CTIO, FA15)
lcosinfa15_param = {
    # telescope/instrument name
    'telescope_instrument': 'LCOGT(CTIO)/SINISTRO',
    'telescope_keyword': 'LCOSINFA15',      # telescope/instrument keyword
    'observatory_code': 'W86',         # MPC observatory code
    'secpix': (0.389, 0.389),  # pixel size (arcsec)
    # before binning
    'ext_coeff': 0.05,          # typical extinction coefficient


    # image orientation preferences
    'flipx': False,
    'flipy': True,
    'rotate': 0,

    # instrument-specific FITS header keywords
    'binning': ('CCDSUM#blank0', 'CCDSUM#blank1'),  # binning in x/y
    'extent': ('NAXIS1', 'NAXIS2'),   # N_pixels in x/y
    'ra': 'RA',  # telescope pointing, RA
    'dec': 'DEC',  # telescope pointin, Dec
    'radec_separator': ':',   # RA/Dec hms separator, use 'XXX'
    # if already in degrees
    'date_keyword': 'DATE-OBS',  # obs date/time
    # keyword; use
    # 'date|time' if
    # separate
    'obsmidtime_jd': 'MJD-OBS',  # obs midtime jd keyword
    # (usually provided by
    # pp_prepare
    'object': 'OBJECT',  # object name keyword
    'filter': 'FILTER',  # filter keyword
    'filter_translations': {'gp': 'g', 'rp': 'r',
                            'ip': 'i', 'zp': 'z', 'w': None},
    # filtername translation dictionary
    'exptime': 'EXPTIME',  # exposure time keyword (s)
    'airmass': 'AIRMASS',  # airmass keyword


    # source extractor settings
    'source_minarea': 9,  # default sextractor source minimum N_pixels
    'source_snr': 3,  # default sextractor source snr for registration
    'aprad_default': 5,  # default aperture radius in px
    'aprad_range': [2, 10],  # [minimum, maximum] aperture radius (px)
    'sex-config-file': rootpath+'/setup/lcosin.sex',
    'mask_file': {},
    #                        mask files as a function of x,y binning

    # registration settings (Scamp)
    'scamp-config-file': rootpath+'/setup/lcosin.scamp',
    'reg_max_mag': 18,
    'reg_search_radius': 0.5,  # deg
    'source_tolerance': 'high',

    # swarp settings
    'copy_keywords': ('TELESCOP,INSTRUME,FILTER,EXPTIME,OBJECT,' +
                      'DATE-OBS,RA,DEC,SECPIX,AIRMASS,' +
                      'TEL_KEYW,MIDTIMJD'),
    #                         keywords to be copied in image
    #                         combination using swarp
    'swarp-config-file': rootpath+'/setup/vatt4k.swarp',

    # default catalog settings
    'astrometry_catalogs': ['GAIA'],
    'photometry_catalogs': ['PANSTARRS', 'SDSS-R9', 'APASS9']
}

# LCOGT, SINISTRO camera (SSO, FL11)
lcosinfl11_param = {
    'telescope_instrument': 'LCOGT(SSO)/SINISTRO',  # telescope/instrument name
    'telescope_keyword': 'LCOSINFL11',      # telescope/instrument keyword
    'observatory_code': 'Q63',         # MPC observatory code
    'secpix': (0.389, 0.389),  # pixel size (arcsec)
    # before binning
    'ext_coeff': 0.05,          # typical extinction coefficient


    # image orientation preferences
    'flipx': False,
    'flipy': True,
    'rotate': 0,

    # instrument-specific FITS header keywords
    'binning': ('CCDSUM#blank0', 'CCDSUM#blank1'),  # binning in x/y
    'extent': ('NAXIS1', 'NAXIS2'),   # N_pixels in x/y
    'ra': 'RA',  # telescope pointing, RA
    'dec': 'DEC',  # telescope pointin, Dec
    'radec_separator': ':',   # RA/Dec hms separator, use 'XXX'
    # if already in degrees
    'date_keyword': 'DATE-OBS',  # obs date/time
    # keyword; use
    # 'date|time' if
    # separate
    'obsmidtime_jd': 'MJD-OBS',  # obs midtime jd keyword
    # (usually provided by
    # pp_prepare
    'object': 'OBJECT',  # object name keyword
    'filter': 'FILTER',  # filter keyword
    'filter_translations': {'gp': 'g', 'rp': 'r',
                            'ip': 'i', 'zp': 'z', 'w': None},
    # filtername translation dictionary
    'exptime': 'EXPTIME',  # exposure time keyword (s)
    'airmass': 'AIRMASS',  # airmass keyword


    # source extractor settings
    'source_minarea': 9,  # default sextractor source minimum N_pixels
    'source_snr': 3,  # default sextractor source snr for registration
    'aprad_default': 5,  # default aperture radius in px
    'aprad_range': [2, 10],  # [minimum, maximum] aperture radius (px)
    'sex-config-file': rootpath+'/setup/lcosin.sex',
    'mask_file': {},
    #                        mask files as a function of x,y binning

    # registration settings (Scamp)
    'scamp-config-file': rootpath+'/setup/lcosin.scamp',
    'reg_max_mag': 18,
    'reg_search_radius': 0.5,  # deg
    'source_tolerance': 'high',

    # swarp settings
    'copy_keywords': ('TELESCOP,INSTRUME,FILTER,EXPTIME,OBJECT,' +
                      'DATE-OBS,RA,DEC,SECPIX,AIRMASS,' +
                      'TEL_KEYW,MIDTIMJD'),
    #                         keywords to be copied in image
    #                         combination using swarp
    'swarp-config-file': rootpath+'/setup/vatt4k.swarp',

    # default catalog settings
    'astrometry_catalogs': ['GAIA'],
    'photometry_catalogs': ['PANSTARRS', 'SDSS-R9', 'APASS9']
}


# LCOGT, SINISTRO camera (SAAO, FL16)
lcosinfl16_param = {
    # telescope/instrument name
    'telescope_instrument': 'LCOGT(SAAO)/SINISTRO',
    'telescope_keyword': 'LCOSINFL16',      # telescope/instrument keyword
    'observatory_code': 'K92',         # MPC observatory code
    'secpix': (0.389, 0.389),  # pixel size (arcsec)
    # before binning
    'ext_coeff': 0.05,          # typical extinction coefficient


    # image orientation preferences
    'flipx': False,
    'flipy': True,
    'rotate': 0,

    # instrument-specific FITS header keywords
    'binning': ('CCDSUM#blank0', 'CCDSUM#blank1'),  # binning in x/y
    'extent': ('NAXIS1', 'NAXIS2'),   # N_pixels in x/y
    'ra': 'RA',  # telescope pointing, RA
    'dec': 'DEC',  # telescope pointin, Dec
    'radec_separator': ':',   # RA/Dec hms separator, use 'XXX'
    # if already in degrees
    'date_keyword': 'DATE-OBS',  # obs date/time
    # keyword; use
    # 'date|time' if
    # separate
    'obsmidtime_jd': 'MJD-OBS',  # obs midtime jd keyword
    # (usually provided by
    # pp_prepare
    'object': 'OBJECT',  # object name keyword
    'filter': 'FILTER',  # filter keyword
    'filter_translations': {'gp': 'g', 'rp': 'r',
                            'ip': 'i', 'zp': 'z', 'w': None},
    # filtername translation dictionary
    'exptime': 'EXPTIME',  # exposure time keyword (s)
    'airmass': 'AIRMASS',  # airmass keyword


    # source extractor settings
    'source_minarea': 9,  # default sextractor source minimum N_pixels
    'source_snr': 3,  # default sextractor source snr for registration
    'aprad_default': 5,  # default aperture radius in px
    'aprad_range': [2, 10],  # [minimum, maximum] aperture radius (px)
    'sex-config-file': rootpath+'/setup/lcosin.sex',
    'mask_file': {},
    #                        mask files as a function of x,y binning

    # registration settings (Scamp)
    'scamp-config-file': rootpath+'/setup/lcosin.scamp',
    'reg_max_mag': 18,
    'reg_search_radius': 0.5,  # deg
    'source_tolerance': 'high',

    # swarp settings
    'copy_keywords': ('TELESCOP,INSTRUME,FILTER,EXPTIME,OBJECT,' +
                      'DATE-OBS,RA,DEC,SECPIX,AIRMASS,' +
                      'TEL_KEYW,MIDTIMJD'),
    #                         keywords to be copied in image
    #                         combination using swarp
    'swarp-config-file': rootpath+'/setup/vatt4k.swarp',

    # default catalog settings
    'astrometry_catalogs': ['GAIA'],
    'photometry_catalogs': ['PANSTARRS', 'SDSS-R9', 'APASS9']
}


# LCOGT, SINISTRO camera (LSC, FL03)
lcosinfl03_param = {
    'telescope_instrument': 'LCOGT(LSC)/SINISTRO',  # telescope/instrument name
    'telescope_keyword': 'LCOSINFL03',      # telescope/instrument keyword
    'observatory_code': 'W85',         # MPC observatory code
    'secpix': (0.387, 0.387),  # pixel size (arcsec)
    # before binning
    'ext_coeff': 0.05,          # typical extinction coefficient


    # image orientation preferences
    'flipx': True,
    'flipy': True,
    'rotate': 0,

    # instrument-specific FITS header keywords
    'binning': ('CCDSUM#blank0', 'CCDSUM#blank1'),  # binning in x/y
    'extent': ('NAXIS1', 'NAXIS2'),   # N_pixels in x/y
    'ra': 'RA',  # telescope pointing, RA
    'dec': 'DEC',  # telescope pointin, Dec
    'radec_separator': ':',   # RA/Dec hms separator, use 'XXX'
    # if already in degrees
    'date_keyword': 'DATE-OBS',  # obs date/time
    # keyword; use
    # 'date|time' if
    # separate
    'obsmidtime_jd': 'MIDTIMJD',  # obs midtime jd keyword
    # (usually provided by
    # pp_prepare
    'object': 'OBJECT',  # object name keyword
    'filter': 'FILTER',  # filter keyword
    'filter_translations': {'gp': 'g', 'rp': 'r',
                            'ip': 'i', 'zp': 'z', 'w': None},
    # filtername translation dictionary
    'exptime': 'EXPTIME',  # exposure time keyword (s)
    'airmass': 'AIRMASS',  # airmass keyword


    # source extractor settings
    'source_minarea': 9,  # default sextractor source minimum N_pixels
    'source_snr': 3,  # default sextractor source snr for registration
    'aprad_default': 5,  # default aperture radius in px
    'aprad_range': [2, 10],  # [minimum, maximum] aperture radius (px)
    'sex-config-file': rootpath+'/setup/lcosin.sex',
    'mask_file': {},
    #                        mask files as a function of x,y binning

    # registration settings (Scamp)
    'scamp-config-file': rootpath+'/setup/lcosin.scamp',
    'reg_max_mag': 18,
    'reg_search_radius': 0.5,  # deg
    'source_tolerance': 'high',

    # swarp settings
    'copy_keywords': ('TELESCOP,INSTRUME,FILTER,EXPTIME,OBJECT,' +
                      'DATE-OBS,RA,DEC,SECPIX,AIRMASS,' +
                      'TEL_KEYW,MIDTIMJD'),
    #                         keywords to be copied in image
    #                         combination using swarp
    'swarp-config-file': rootpath+'/setup/vatt4k.swarp',

    # default catalog settings
    'astrometry_catalogs': ['GAIA'],
    'photometry_catalogs': ['PANSTARRS', 'SDSS-R9', 'APASS9']
}

# LCOGT, Spectral camera (COJ, FS01)
lcospecfs01_param = {
    'telescope_instrument': 'LCOGT(COJ)/SPECTRAL',  # telescope/instrument name
    'telescope_keyword': 'LCOSPECFS01',      # telescope/instrument keyword
    'observatory_code': '413',         # MPC observatory code
    'secpix': (0.15, 0.15),  # pixel size (arcsec)
    # before binning
    'ext_coeff': 0.05,          # typical extinction coefficient


    # image orientation preferences
    'flipx': True,
    'flipy': False,
    'rotate': 0,

    # instrument-specific FITS header keywords
    'binning': ('CCDSUM#blank0', 'CCDSUM#blank1'),  # binning in x/y
    'extent': ('NAXIS1', 'NAXIS2'),   # N_pixels in x/y
    'ra': 'RA',  # telescope pointing, RA
    'dec': 'DEC',  # telescope pointin, Dec
    'radec_separator': ':',   # RA/Dec hms separator, use 'XXX'
    # if already in degrees
    'date_keyword': 'DATE-OBS',  # obs date/time
    # keyword; use
    # 'date|time' if
    # separate
    'obsmidtime_jd': 'MIDTIMJD',  # obs midtime jd keyword
    # (usually provided by
    # pp_prepare
    'object': 'OBJECT',  # object name keyword
    'filter': 'FILTER',  # filter keyword
    'filter_translations': {'gp': 'g', 'rp': 'r',
                            'ip': 'i', 'zp': 'z',
                            'clear': None,
                            'V': 'V'},
    # filtername translation dictionary
    'exptime': 'EXPTIME',  # exposure time keyword (s)
    'airmass': 'AIRMASS',  # airmass keyword


    # source extractor settings
    'source_minarea': 9,  # default sextractor source minimum N_pixels
    'source_snr': 3,  # default sextractor source snr for registration
    'aprad_default': 5,  # default aperture radius in px
    'aprad_range': [2, 10],  # [minimum, maximum] aperture radius (px)
    'sex-config-file': rootpath+'/setup/lcospec.sex',
    'mask_file': {},
    #                        mask files as a function of x,y binning

    # registration settings (Scamp)
    'scamp-config-file': rootpath+'/setup/lcospec.scamp',
    'reg_max_mag': 18,
    'reg_search_radius': 0.5,  # deg
    'source_tolerance': 'high',

    # swarp settings
    'copy_keywords': ('TELESCOP,INSTRUME,FILTER,EXPTIME,OBJECT,' +
                      'DATE-OBS,RA,DEC,SECPIX,AIRMASS,' +
                      'TEL_KEYW,MIDTIMJD'),
    #                         keywords to be copied in image
    #                         combination using swarp
    'swarp-config-file': rootpath+'/setup/vatt4k.swarp',

    # default catalog settings
    'astrometry_catalogs': ['GAIA'],
    'photometry_catalogs': ['PANSTARRS', 'SDSS-R9', 'APASS9']
}


# LCOGT, SINISTRO camera (SAAO, FL06)
lcosinfl06_param = {
    # telescope/instrument name
    'telescope_instrument': 'LCOGT(SAAO)/SINISTRO',
    'telescope_keyword': 'LCOSINFL06',      # telescope/instrument keyword
    'observatory_code': 'K92',         # MPC observatory code
    'secpix': (0.389, 0.389),  # pixel size (arcsec)
    # before binning
    'ext_coeff': 0.05,          # typical extinction coefficient


    # image orientation preferences
    'flipx': False,
    'flipy': True,
    'rotate': 0,

    # instrument-specific FITS header keywords
    'binning': ('CCDSUM#blank0', 'CCDSUM#blank1'),  # binning in x/y
    'extent': ('NAXIS1', 'NAXIS2'),   # N_pixels in x/y
    'ra': 'RA',  # telescope pointing, RA
    'dec': 'DEC',  # telescope pointin, Dec
    'radec_separator': ':',   # RA/Dec hms separator, use 'XXX'
    # if already in degrees
    'date_keyword': 'DATE-OBS',  # obs date/time
    # keyword; use
    # 'date|time' if
    # separate
    'obsmidtime_jd': 'MJD-OBS',  # obs midtime jd keyword
    # (usually provided by
    # pp_prepare
    'object': 'OBJECT',  # object name keyword
    'filter': 'FILTER',  # filter keyword
    'filter_translations': {'gp': 'g', 'rp': 'r',
                            'ip': 'i', 'zp': 'z', 'w': None},
    # filtername translation dictionary
    'exptime': 'EXPTIME',  # exposure time keyword (s)
    'airmass': 'AIRMASS',  # airmass keyword

    # source extractor settings
    'source_minarea': 9,  # default sextractor source minimum N_pixels
    'source_snr': 3,  # default sextractor source snr for registration
    'aprad_default': 5,  # default aperture radius in px
    'aprad_range': [2, 10],  # [minimum, maximum] aperture radius (px)
    'sex-config-file': rootpath+'/setup/lcosin.sex',
    'mask_file': {},
    #                        mask files as a function of x,y binning

    # registration settings (Scamp)
    'scamp-config-file': rootpath+'/setup/lcosin.scamp',
    'reg_max_mag': 18,
    'reg_search_radius': 0.5,  # deg
    'source_tolerance': 'high',

    # swarp settings
    'copy_keywords': ('TELESCOP,INSTRUME,FILTER,EXPTIME,OBJECT,' +
                      'DATE-OBS,RA,DEC,SECPIX,AIRMASS,' +
                      'TEL_KEYW,MIDTIMJD'),
    #                         keywords to be copied in image
    #                         combination using swarp
    'swarp-config-file': rootpath+'/setup/vatt4k.swarp',

    # default catalog settings
    'astrometry_catalogs': ['GAIA'],
    'photometry_catalogs': ['PANSTARRS', 'SDSS-R9', 'APASS9']
}


# Apache Point Observatory 3.5m, SPICAM
arc35spicam_param = {
    'telescope_instrument': 'ARC35/SPICAM',  # telescope/instrument name
    'telescope_keyword': 'ARC35SPICAM',      # telescope/instrument keyword
    'observatory_code': '705',         # MPC observatory code
    'secpix': (0.141, 0.141),  # pixel size (arcsec)
    # before binning
    'ext_coeff': 0.05,          # typical extinction coefficient


    # image orientation preferences
    'flipx': True,
    'flipy': False,
    'rotate': 0,

    # instrument-specific FITS header keywords
    'binning': ('CCDSUM#blank0', 'CCDSUM#blank1'),  # binning in x/y
    'extent': ('NAXIS1', 'NAXIS2'),   # N_pixels in x/y
    'ra': 'RA',  # telescope pointing, RA
    'dec': 'DEC',  # telescope pointin, Dec
    'radec_separator': ':',   # RA/Dec hms separator, use 'XXX'
    # if already in degrees
    'date_keyword': 'DATE-OBS',  # obs date/time
    # keyword; use
    # 'date|time' if
    # separate
    'obsmidtime_jd': 'MIDTIMJD',  # obs midtime jd keyword
    # (usually provided by
    # pp_prepare
    'object': 'OBJECT',  # object name keyword
    'filter': 'FILTER',  # filter keyword
    'filter_translations': {'SDSS g': 'g', 'SDSS r': 'r',
                            'SDSS i': 'i', 'SDSS z': 'z'},
    # filtername translation dictionary
    'exptime': 'EXPTIME',  # exposure time keyword (s)
    'airmass': 'AIRMASS',  # airmass keyword


    # source extractor settings
    'source_minarea': 15,  # default sextractor source minimum N_pixels
    'source_snr': 3,  # default sextractor source snr for registration
    'aprad_default': 5,  # default aperture radius in px
    'aprad_range': [3, 15],  # [minimum, maximum] aperture radius (px)
    'sex-config-file': rootpath+'/setup/arc35spicam.sex',
    'mask_file': {},
    #                        mask files as a function of x,y binning

    # registration settings (Scamp)
    'scamp-config-file': rootpath+'/setup/arc35spicam.scamp',
    'reg_max_mag': 18,
    'reg_search_radius': 0.5,  # deg
    'source_tolerance': 'high',

    # swarp settings
    'copy_keywords': ('TELESCOP,INSTRUME,FILTER,EXPTIME,OBJECT,' +
                      'DATE-OBS,RA,DEC,SECPIX,AIRMASS,' +
                      'TEL_KEYW,MIDTIMJD'),
    #                         keywords to be copied in image
    #                         combination using swarp
    'swarp-config-file': rootpath+'/setup/vatt4k.swarp',

    # default catalog settings
    'astrometry_catalogs': ['GAIA'],
    'photometry_catalogs': ['PANSTARRS', 'SDSS-R9', 'APASS9']
}

# LCOGT, SINISTRO camera (LSC, FL03)
lcosinfl03_param = {
    'telescope_instrument': 'LCOGT(LSC)/SINISTRO',  # telescope/instrument name
    'telescope_keyword': 'LCOSINFL03',      # telescope/instrument keyword
    'observatory_code': 'W85',         # MPC observatory code
    'secpix': (0.387, 0.387),  # pixel size (arcsec)
    # before binning
    'ext_coeff': 0.05,          # typical extinction coefficient


    # image orientation preferences
    'flipx': True,
    'flipy': True,
    'rotate': 0,

    # instrument-specific FITS header keywords
    'binning': ('CCDSUM#blank0', 'CCDSUM#blank1'),  # binning in x/y
    'extent': ('NAXIS1', 'NAXIS2'),   # N_pixels in x/y
    'ra': 'RA',  # telescope pointing, RA
    'dec': 'DEC',  # telescope pointin, Dec
    'radec_separator': ':',   # RA/Dec hms separator, use 'XXX'
    # if already in degrees
    'date_keyword': 'DATE-OBS',  # obs date/time
    # keyword; use
    # 'date|time' if
    # separate
    'obsmidtime_jd': 'MIDTIMJD',  # obs midtime jd keyword
    # (usually provided by
    # pp_prepare
    'object': 'OBJECT',  # object name keyword
    'filter': 'FILTER',  # filter keyword
    'filter_translations': {'gp': 'g', 'rp': 'r',
                            'ip': 'i', 'zp': 'z'},
    # filtername translation dictionary
    'exptime': 'EXPTIME',  # exposure time keyword (s)
    'airmass': 'AIRMASS',  # airmass keyword


    # source extractor settings
    'source_minarea': 9,  # default sextractor source minimum N_pixels
    'source_snr': 3,  # default sextractor source snr for registration
    'aprad_default': 5,  # default aperture radius in px
    'aprad_range': [2, 10],  # [minimum, maximum] aperture radius (px)
    'sex-config-file': rootpath+'/setup/lcosin.sex',
    'mask_file': {},
    #                        mask files as a function of x,y binning

    # registration settings (Scamp)
    'scamp-config-file': rootpath+'/setup/lcosin.scamp',
    'reg_max_mag': 18,
    'reg_search_radius': 0.5,  # deg
    'source_tolerance': 'high',

    # swarp settings
    'copy_keywords': ('TELESCOP,INSTRUME,FILTER,EXPTIME,OBJECT,' +
                      'DATE-OBS,RA,DEC,SECPIX,AIRMASS,' +
                      'TEL_KEYW,MIDTIMJD'),
    #                         keywords to be copied in image
    #                         combination using swarp
    'swarp-config-file': rootpath+'/setup/vatt4k.swarp',

    # default catalog settings
    'astrometry_catalogs': ['GAIA'],
    'photometry_catalogs': ['PANSTARRS', 'SDSS-R9', 'APASS9']
}

# LCOGT, Spectral camera (COJ, FS01)
lcospecfs01_param = {
    'telescope_instrument': 'LCOGT(COJ)/SPECTRAL',  # telescope/instrument name
    'telescope_keyword': 'LCOSPECFS01',      # telescope/instrument keyword
    'observatory_code': '413',         # MPC observatory code
    'secpix': (0.15, 0.15),  # pixel size (arcsec)
    # before binning
    'ext_coeff': 0.05,          # typical extinction coefficient


    # image orientation preferences
    'flipx': True,
    'flipy': False,
    'rotate': 0,

    # instrument-specific FITS header keywords
    'binning': ('CCDSUM#blank0', 'CCDSUM#blank1'),  # binning in x/y
    'extent': ('NAXIS1', 'NAXIS2'),   # N_pixels in x/y
    'ra': 'RA',  # telescope pointing, RA
    'dec': 'DEC',  # telescope pointin, Dec
    'radec_separator': ':',   # RA/Dec hms separator, use 'XXX'
    # if already in degrees
    'date_keyword': 'DATE-OBS',  # obs date/time
    # keyword; use
    # 'date|time' if
    # separate
    'obsmidtime_jd': 'MIDTIMJD',  # obs midtime jd keyword
    # (usually provided by
    # pp_prepare
    'object': 'OBJECT',  # object name keyword
    'filter': 'FILTER',  # filter keyword
    'filter_translations': {'gp': 'g', 'rp': 'r',
                            'ip': 'i', 'zp': 'z',
                            'clear': None,
                            'V': 'V'},
    # filtername translation dictionary
    'exptime': 'EXPTIME',  # exposure time keyword (s)
    'airmass': 'AIRMASS',  # airmass keyword


    # source extractor settings
    'source_minarea': 9,  # default sextractor source minimum N_pixels
    'source_snr': 3,  # default sextractor source snr for registration
    'aprad_default': 5,  # default aperture radius in px
    'aprad_range': [2, 10],  # [minimum, maximum] aperture radius (px)
    'sex-config-file': rootpath+'/setup/lcospec.sex',
    'mask_file': {},
    #                        mask files as a function of x,y binning

    # registration settings (Scamp)
    'scamp-config-file': rootpath+'/setup/lcospec.scamp',
    'reg_max_mag': 18,
    'reg_search_radius': 0.5,  # deg
    'source_tolerance': 'high',

    # swarp settings
    'copy_keywords': ('TELESCOP,INSTRUME,FILTER,EXPTIME,OBJECT,' +
                      'DATE-OBS,RA,DEC,SECPIX,AIRMASS,' +
                      'TEL_KEYW,MIDTIMJD'),
    #                         keywords to be copied in image
    #                         combination using swarp
    'swarp-config-file': rootpath+'/setup/vatt4k.swarp',

    # default catalog settings
    'astrometry_catalogs': ['GAIA'],
    'photometry_catalogs': ['PANSTARRS', 'SDSS-R9', 'APASS9']
}

# Palomar 60-inch, optical facility camera
p60opt_param = {
    'telescope_instrument': 'Palomar60-inch/opticalfacility',  # telescope/instrument name
    'telescope_keyword': 'P60OPT',      # telescope/instrument keyword
    'observatory_code': '675',         # MPC observatory code
    'secpix': (0.378, 0.378),  # pixel size (arcsec)
    # before binning
    'ext_coeff': 0.05,          # typical extinction coefficient


    # image orientation preferences
    'flipx': True,
    'flipy': False,
    'rotate': 270,

    # instrument-specific FITS header keywords
    'binning': ('CCDSUM#blank0', 'CCDSUM#blank1'),  # binning in x/y
    'extent': ('NAXIS1', 'NAXIS2'),   # N_pixels in x/y
    'ra': 'RA',  # telescope pointing, RA
    'dec': 'DEC',  # telescope pointin, Dec
    'radec_separator': ':',   # RA/Dec hms separator, use 'XXX'
    # if already in degrees
    'date_keyword': 'UTSHUT',  # obs date/time
    # keyword; use
    # 'date|time' if
    # separate
    'obsmidtime_jd': 'MIDTIMJD',  # obs midtime jd keyword
    # (usually provided by
    # pp_prepare
    'object': 'OBJECT',  # object name keyword
    'filter': 'FILTER',  # filter keyword
    'filter_translations': {'V': 'V', 'R': 'R'},
    # filtername translation dictionary
    'exptime': 'EXPTIME',  # exposure time keyword (s)
    'airmass': 'AIRMASS',  # airmass keyword


    # source extractor settings
    'source_minarea': 9,  # default sextractor source minimum N_pixels
    'source_snr': 3,  # default sextractor source snr for registration
    'aprad_default': 5,  # default aperture radius in px
    'aprad_range': [2, 10],  # [minimum, maximum] aperture radius (px)
    'sex-config-file': rootpath+'/setup/p60opt.sex',
    'mask_file': {'1,1': rootpath+'/setup/mask_p60opt_1x1.fits'},
    #                        mask files as a function of x,y binning

    # registration settings (Scamp)
    'scamp-config-file': rootpath+'/setup/p60opt.scamp',
    'reg_max_mag': 18,
    'reg_search_radius': 0.5,  # deg
    'source_tolerance': 'high',

    # swarp settings
    'copy_keywords': ('TELESCOP,INSTRUME,FILTER,EXPTIME,OBJECT,' +
                      'UTSHUT,RA,DEC,SECPIX,AIRMASS,' +
                      'TEL_KEYW,MIDTIMJD'),
    #                         keywords to be copied in image
    #                         combination using swarp
    'swarp-config-file': rootpath+'/setup/vatt4k.swarp',

    # default catalog settings
    'astrometry_catalogs': ['GAIA'],
    'photometry_catalogs': ['PANSTARRS', 'SDSS-R9', 'APASS9']
}

# Palomar 60-inch, SED Machine
p60sedm_param = {
    'telescope_instrument': 'Palomar60-inch/SEDmachine',  # telescope/instrument name
    'telescope_keyword': 'P60SEDM',      # telescope/instrument keyword
    'observatory_code': '675',         # MPC observatory code
    'secpix': (0.38, 0.38),  # pixel size (arcsec)
    # before binning
    'ext_coeff': 0.05,          # typical extinction coefficient


    # image orientation preferences
    'flipx': True,
    'flipy': True,
    'rotate': 0,

    # instrument-specific FITS header keywords
    'binning': (1, 1),  # binning in x/y
    'extent': ('NAXIS1', 'NAXIS2'),   # N_pixels in x/y
    'ra': 'RA',  # telescope pointing, RA
    'dec': 'DEC',  # telescope pointin, Dec
    'radec_separator': ':',   # RA/Dec hms separator, use 'XXX'
    # if already in degrees
    'date_keyword': 'OBSDATE|OBSTIME',  # obs date/time
    # keyword; use
    # 'date|time' if
    # separate
    'obsmidtime_jd': 'MIDTIMJD',  # obs midtime jd keyword
    # (usually provided by
    # pp_prepare
    'object': 'OBJECT',  # object name keyword
    'filter': 'FILTER',  # filter keyword
    'filter_translations': {'g': 'g', 'r': 'r',
                            'i': 'i', 'z': 'z'},
    # filtername translation dictionary
    'exptime': 'EXPTIME',  # exposure time keyword (s)
    'airmass': 'AIRMASS',  # airmass keyword


    # source extractor settings
    'source_minarea': 7,  # default sextractor source minimum N_pixels
    'source_snr': 2,  # default sextractor source snr for registration
    'aprad_default': 5,  # default aperture radius in px
    'aprad_range': [2, 10],  # [minimum, maximum] aperture radius (px)
    'sex-config-file': rootpath+'/setup/p60sedm.sex',
    'mask_file': {},
    #                        mask files as a function of x,y binning

    # registration settings (Scamp)
    'scamp-config-file': rootpath+'/setup/p60sedm.scamp',
    'reg_max_mag': 20,
    'reg_search_radius': 0.3,  # deg
    'source_tolerance': 'high',

    # swarp settings
    'copy_keywords': ('TELESCOP,INSTRUME,FILTER,EXPTIME,OBJECT,' +
                      'OBSDATE,OBSTIME,RA,DEC,SECPIX,AIRMASS,' +
                      'TEL_KEYW,MIDTIMJD'),
    #                         keywords to be copied in image
    #                         combination using swarp
    'swarp-config-file': rootpath+'/setup/vatt4k.swarp',

    # default catalog settings
    'astrometry_catalogs': ['GAIA'],
    'photometry_catalogs': ['PANSTARRS', 'SDSS-R9', 'APASS9']
}


# Gemini North, GMOS
gmosn_param = {
    'telescope_instrument': 'Gemini-N/GMOS',  # telescope/instrument name
    'telescope_keyword': 'GMOSN',  # telescope/instrument keyword
    'observatory_code': '568',         # MPC observatory code
    'secpix': (0.081, 0.081),  # pixel size (arcsec)
    # before binning
    'ext_coeff': 0.05,          # typical extinction coefficient


    # image orientation preferences
    'flipx': False,
    'flipy': False,
    'rotate': 90,

    # instrument-specific FITS header keywords
    'binning': ('CCDSUM#blank0', 'CCDSUM#blank1'),
    # binning in x/y, '_blankN' denotes that both axes
    # are listed in one keyword, sep. by blanks
    'extent': ('NAXIS1', 'NAXIS2'),   # N_pixels in x/y
    'ra': 'CRVAL1',  # telescope pointing, RA
    'dec': 'CRVAL2',  # telescope pointing, Dec
    'radec_separator': 'XXX',   # RA/Dec hms separator, use 'XXX'
    # if already in degrees
    'date_keyword': 'DATE-OBS|UTSTART',  # obs date/time
    # keyword; guse
    # 'date|time' if
    # separate
    'obsmidtime_jd': 'MIDTIMJD',  # obs midtime jd keyword
    # (usually provided by
    # pp_prepare
    'object': 'OBJECT',  # object name keyword
    'filter': 'FILTER2',  # filter keyword
    'filter_translations': {'r_G0303': 'r', 'i_G0302': 'i', 'clear': None},
    # filtername translation dictionary
    'exptime': 'EXPTIME',  # exposure time keyword (s)
    'airmass': 'AIRMASS',  # airmass keyword


    # source extractor settings
    'source_minarea': 9,  # default sextractor source minimum N_pixels
    'source_snr': 10,  # default sextractor source snr for registration
    'aprad_default': 6,  # default aperture radius in px
    'aprad_range': [2, 15],  # [minimum, maximum] aperture radius (px)
    'sex-config-file': rootpath+'/setup/gmosn.sex',
    'mask_file': {},  # '2,2': rootpath+'/setup/gmosn_mask_2x2.fits'},
    #                #        mask files as a function of x,y binning

    # registration settings (Scamp)
    'scamp-config-file': rootpath+'/setup/gmosn.scamp',
    'reg_max_mag': 23,
    'reg_search_radius': 0.3,  # deg
    'source_tolerance': 'high',

    # swarp settings
    'copy_keywords': ('OBSERVAT,INSTRUME,EXPTIME,OBJECT,' +
                      'DATE-OBS,UTSTART,RA,DEC,AIRMASS,TEL_KEYW,CCDSUM,' +
                      'FILTER2,MIDTIMJD'),
    #                        keywords to be copied in image
    #                        combination using swarp
    'swarp-config-file': rootpath+'/setup/gmosn.swarp',

    # default catalog settings
    'astrometry_catalogs': ['GAIA'],
    'photometry_catalogs': ['SDSS-R9', 'PANSTARRS', 'APASS9']
}

# Danish 1.54m, DFOSC
dfosc_param = {
    'telescope_instrument': 'Danish/DFOSC',  # telescope/instrument name
    'telescope_keyword': 'DFOSC',      # telescope/instrument keyword
    'observatory_code': '809',         # MPC observatory code
    'secpix': (0.4, 0.4),  # pixel size (arcsec)
    # before binning
    'ext_coeff': 0.05,          # typical extinction coefficient


    # image orientation preferences
    'flipx': True,
    'flipy': False,
    'rotate': 180,

    # instrument-specific FITS header keywords
    'binning': ('BINX', 'BINY'),  # binning in x/y
    'extent': ('NAXIS1', 'NAXIS2'),   # N_pixels in x/y
    'ra': 'OBJRA',  # telescope pointing, RA
    'dec': 'OBJDEC',  # telescope pointin, Dec
    'radec_separator': ':',   # RA/Dec hms separator, use 'XXX'
    # if already in degrees
    'date_keyword': 'DATE-OBS',  # obs date/time
    # keyword; use
    # 'date|time' if
    # separate
    'obsmidtime_jd': 'MIDTIMJD',  # obs midtime jd keyword
    # (usually provided by
    # pp_prepare
    'object': 'OBJECT',  # object name keyword
    'filter': 'FILTB',  # filter keyword
    'filter_translations': {'V': 'V', 'R': 'R',
                            'I': 'I', 'B': 'B'},
    # filtername translation dictionary
    'exptime': 'EXPTIME',  # exposure time keyword (s)
    'airmass': 'AIRMASS',  # airmass keyword


    # source extractor settings
    'source_minarea': 9,  # default sextractor source minimum N_pixels
    'source_snr': 3,  # default sextractor source snr for registration
    'aprad_default': 5,  # default aperture radius in px
    'aprad_range': [2, 10],  # [minimum, maximum] aperture radius (px)
    'sex-config-file': rootpath+'/setup/dfosc.sex',
    'mask_file': {},
    #                        mask files as a function of x,y binning

    # registration settings (Scamp)
    'scamp-config-file': rootpath+'/setup/dfosc.scamp',
    'reg_max_mag': 19,
    'reg_search_radius': 0.5,  # deg
    'source_tolerance': 'high',

    # swarp settings
    'copy_keywords': ('TELESCOP,INSTRUME,FILTB,EXPTIME,OBJECT,' +
                      'DATE-OBS,OBJRA,OBJDEC,SECPIX,AIRMASS,' +
                      'TEL_KEYW,BINX,BINY,MIDTIMJD'),
    #                         keywords to be copied in image
    #                         combination using swarp
    'swarp-config-file': rootpath+'/setup/dfosc.swarp',

    # default catalog settings
    'astrometry_catalogs': ['GAIA'],
    'photometry_catalogs': ['PANSTARRS', 'SDSS-R9', 'APASS9']
}

# Lowell Near-Earth Object Survey (LONEOS)
loneos_param = {
    'telescope_instrument': 'LONEOS',  # telescope/instrument name
    'telescope_keyword': 'LONEOS',      # telescope/instrument keyword
    'observatory_code': '699',         # MPC observatory code
    'secpix': (2.53, 2.53),  # pixel size (arcsec)
    # before binning
    'ext_coeff': 0.05,          # typical extinction coefficient


    # image orientation preferences
    'flipx': False,
    'flipy': False,
    'rotate': 0,

    # instrument-specific FITS header keywords
    'binning': (1, 1),  # binning in x/y
    'extent': ('NAXIS1', 'NAXIS2'),   # N_pixels in x/y
    'ra': 'TELRA',  # telescope pointing, RA
    'dec': 'TELDEC',  # telescope pointin, Dec
    'radec_separator': ':',   # RA/Dec hms separator, use 'XXX'
    # if already in degrees
    'date_keyword': 'DATE-OBS',  # obs date/time
    # keyword; use
    # 'date|time' if
    # separate
    'obsmidtime_jd': 'MIDTIMJD',  # obs midtime jd keyword
    # (usually provided by
    # pp_prepare
    'object': 'OBJECT',  # object name keyword
    'filter': 'FILTER',  # filter keyword
    'filter_translations': {'open': None},
    # filtername translation dictionary
    'exptime': 'EXPTIME',  # exposure time keyword (s)
    'airmass': 'AIRMASS',  # airmass keyword


    # source extractor settings
    'source_minarea': 5,  # default sextractor source minimum N_pixels
    'source_snr': 3,  # default sextractor source snr for registration
    'aprad_default': 3,  # default aperture radius in px
    'aprad_range': [1, 8],  # [minimum, maximum] aperture radius (px)
    'sex-config-file': rootpath+'/setup/frost.sex',
    'mask_file': {'1,1': rootpath+'/setup/mask_loneos.fits'},
    #                        mask files as a function of x,y binning

    # registration settings (Scamp)
    'scamp-config-file': rootpath+'/setup/frost.scamp',
    'reg_max_mag': 18,
    'reg_search_radius': 0.5,  # deg
    'source_tolerance': 'high',

    # swarp settings
    'copy_keywords': ('TELESCOP,INSTRUME,FILTER,EXPTIME,OBJECT,' +
                      'DATE-OBS,TIME-OBS,TELRA,TELDEC,SECPIX,AIRMASS,' +
                      'TEL_KEYW'),
    #                         keywords to be copied in image
    #                         combination using swarp
    'swarp-config-file': rootpath+'/setup/frost.swarp',

    # default catalog settings
    'astrometry_catalogs': ['GAIA'],
    'photometry_catalogs': ['SDSS-R9', 'APASS9', 'PANSTARRS', '2MASS']
}

# Palmer Divide Observatory, 25cm f/6.3, SBIG ST-8
pdo25cmf63st8_param = {
    'telescope_instrument': 'PDO 25cm f6.3/SBIG ST-8',  # telescope/instrument name
    'telescope_keyword': 'PDO25CMF63ST8',      # telescope/instrument keyword
    'observatory_code': '716',         # MPC observatory code
    'secpix': (0.1875, 0.1875),  # pixel size (arcsec) before binning

    # image orientation preferences
    'flipx': True,
    'flipy': False,
    'rotate': 0,

    # instrument-specific FITS header keywords
    'binning': ('BINNING', 'BINNING'),  # binning in x/y
    'extent': ('NAXIS1', 'NAXIS2'),   # N_pixels in x/y
    'ra': 'OBJCTRA',  # telescope pointing, RA
    'dec': 'OBJCTDEC',  # telescope pointin, Dec
    'radec_separator': ' ',   # RA/Dec hms separator, use 'XXX'
    # if already in degrees
    'date_keyword': 'DATETIME',  # obs date/time
    # keyword; use 'date|time' if separate
    'obsmidtime_jd': 'MIDTIMJD',  # obs midtime jd keyword
    # (usually provided by pp_prepare
    'object': 'OBJECT',  # object name keyword
    'filter': 'FILTER',  # filter keyword
    'filter_translations': {'clear': None},
    # filtername translation dictionary
    'exptime': 'EXPOSURE',  # exposure time keyword (s)
    'airmass': 'AIRMASS',  # airmass keyword


    # source extractor settings
    'source_minarea': 10,  # default sextractor source minimum N_pixels
    'source_snr': 3,  # default sextractor source snr for registration
    'aprad_default': 5,  # default aperture radius in px
    'aprad_range': [2, 10],  # [minimum, maximum] aperture radius (px)
    'sex-config-file': rootpath+'/setup/pdo.sex',
    'mask_file': {},  # mask files as a function of x,y binning

    # registration settings (Scamp)
    'scamp-config-file': rootpath+'/setup/pdo.scamp',
    'reg_max_mag': 16,
    'reg_search_radius': 0.2,  # deg
    'source_tolerance': 'high',

    # swarp settings
    'copy_keywords': ('TELESCOP,INSTRUME,FILTER,EXPOSURE,OBJECT,' +
                      'DATETIME,OBJCTRA,OBJCTDEC,SECPIX,' +
                      'TEL_KEYW,BINNING,MIDTIMJD'),
    #                         keywords to be copied in image
    #                         combination using swarp
    'swarp-config-file': rootpath+'/setup/pdo.swarp',

    # default catalog settings
    'astrometry_catalogs': ['TGAS'],
    'photometry_catalogs': ['PANSTARRS', 'SDSS-R9', 'APASS9']
}

# Palmer Divide Observatory, 0.5m f/8.1 Ritchey-Chretien, FLI KAF1001E -- misspelling
pdo05mf81kaf1001e_param = {
    'telescope_instrument': 'PDO 0.5m f8.1/FLI KAF1001E',  # telescope/instrument name
    'telescope_keyword': 'POD05MF81KAF1001E',      # telescope/instrument keyword
    'observatory_code': '716',         # MPC observatory code
    'secpix': (1.2, 1.2),  # pixel size (arcsec) before binning
    'ext_coeff': 0.05,          # typical extinction coefficient

    # image orientation preferences
    'flipx': True,
    'flipy': True,
    'rotate': 0,

    # instrument-specific FITS header keywords
    'binning': ('BINNING', 'BINNING'),  # binning in x/y
    'extent': ('NAXIS1', 'NAXIS2'),   # N_pixels in x/y
    'ra': 'OBJCTRA',  # telescope pointing, RA
    'dec': 'OBJCTDEC',  # telescope pointin, Dec
    'radec_separator': ':',   # RA/Dec hms separator, use 'XXX'
    # if already in degrees
    'date_keyword': 'DATE-OBS',  # obs date/time
    # keyword; use
    # 'date|time' if
    # separate
    'obsmidtime_jd': 'MIDTIMJD',  # obs midtime jd keyword
    # (usually provided by
    # pp_prepare
    'object': 'OBJECT',  # object name keyword
    'filter': 'FILTER',  # filter keyword
    'filter_translations': {'Clear': None},
    # filtername translation dictionary
    'exptime': 'EXPOSURE',  # exposure time keyword (s)
    'airmass': 'AIRMASS',  # airmass keyword


    # source extractor settings
    'source_minarea': 12,  # default sextractor source minimum N_pixels
    'source_snr': 3,  # default sextractor source snr for registration
    'aprad_default': 5,  # default aperture radius in px
    'aprad_range': [2, 10],  # [minimum, maximum] aperture radius (px)
    'sex-config-file': rootpath+'/setup/pdo.sex',
    'mask_file': {},
    #                        mask files as a function of x,y binning

    # registration settings (Scamp)
    'scamp-config-file': rootpath+'/setup/pdo.scamp',
    'reg_max_mag': 19,
    'reg_search_radius': 0.5,  # deg
    'source_tolerance': 'high',

    # swarp settings
    'copy_keywords': ('TELESCOP,INSTRUME,FILTER,EXPOSURE,OBJECT,' +
                      'DATE-OBS,OBJCTRA,OBJCTDEC,SECPIX,AIRMASS,' +
                      'TEL_KEYW,BINNING,MIDTIMJD'),
    #                         keywords to be copied in image
    #                         combination using swarp
    'swarp-config-file': rootpath+'/setup/pdo.swarp',

    # default catalog settings
    'astrometry_catalogs': ['GAIA'],
    'photometry_catalogs': ['PANSTARRS', 'SDSS-R9', 'APASS9']
}


# CS3-PDS-2-14N, 35cm SCT, STL-1001E
pds35cmstl1001e_param = {
    # telescope/instrument name
    'telescope_instrument': 'CS3-PDS-2-14N 35cm SCT/STL-1001E',
    'telescope_keyword': 'PDS35CMSTL1001E',      # telescope/instrument keyword
    'observatory_code': 'U82',         # MPC observatory code
    'secpix': (1.48, 1.48),  # pixel size (arcsec) before binning

    # image orientation preferences
    'flipx': False,
    'flipy': False,
    'rotate': 0,

    # instrument-specific FITS header keywords
    'binning': ('XBINNING', 'YBINNING'),  # binning in x/y
    'extent': ('NAXIS1', 'NAXIS2'),   # N_pixels in x/y
    'ra': 'OBJCTRA',  # telescope pointing, RA
    'dec': 'OBJCTDEC',  # telescope pointin, Dec
    'radec_separator': ' ',   # RA/Dec hms separator, use 'XXX'
    # if already in degrees
    'date_keyword': 'DATE-OBS',  # obs date/time
    # keyword; use 'date|time' if separate
    'obsmidtime_jd': 'MIDTIMJD',  # obs midtime jd keyword
    # (usually provided by pp_prepare
    'object': 'OBJECT',  # object name keyword
    'filter': 'FILTER',  # filter keyword
    'filter_translations': {'NF': None},
    # filtername translation dictionary
    'exptime': 'EXPOSURE',  # exposure time keyword (s)
    'airmass': 'AIRMASS',  # airmass keyword


    # source extractor settings
    'source_minarea': 10,  # default sextractor source minimum N_pixels
    'source_snr': 3,  # default sextractor source snr for registration
    'aprad_default': 5,  # default aperture radius in px
    'aprad_range': [2, 10],  # [minimum, maximum] aperture radius (px)
    'sex-config-file': rootpath+'/setup/pds.sex',
    'mask_file': {},  # mask files as a function of x,y binning

    # registration settings (Scamp)
    'scamp-config-file': rootpath+'/setup/pds.scamp',
    'reg_max_mag': 16,
    'reg_search_radius': 0.2,  # deg
    'source_tolerance': 'high',

    # swarp settings
    'copy_keywords': ('TELESCOP,INSTRUME,FILTER,EXPOSURE,OBJECT,' +
                      'DATETIME,OBJCTRA,OBJCTDEC,SECPIX,' +
                      'TEL_KEYW,XBINNING,YBINNING,MIDTIMJD'),
    #                         keywords to be copied in image
    #                         combination using swarp
    'swarp-config-file': rootpath+'/setup/pds.swarp',

    # default catalog settings
    'astrometry_catalogs': ['GAIA'],
    'photometry_catalogs': ['PANSTARRS', 'SDSS-R9', 'APASS9']
}

# MMT, MMTCam
mmtcam_param = {
    'telescope_instrument': 'MMT/MMTCam',  # telescope/instrument name
    'telescope_keyword': 'MMTCAM',      # telescope/instrument keyword
    'observatory_code': '696',         # MPC observatory code
    'secpix': (0.082, 0.082),  # pixel size (arcsec)
    # before binning
    'ext_coeff': 0.05,          # typical extinction coefficient


    # image orientation preferences
    'flipx': True,
    'flipy': True,
    'rotate': 0,

    # instrument-specific FITS header keywords
    'binning': ('XBINNING', 'YBINNING'),  # binning in x/y
    'extent': ('NAXIS1', 'NAXIS2'),   # N_pixels in x/y
    'ra': 'OBJCTRA',  # telescope pointing, RA
    'dec': 'OBJCTDEC',  # telescope pointin, Dec
    'radec_separator': ' ',   # RA/Dec hms separator, use 'XXX'
    # if already in degrees
    'date_keyword': 'DATE-OBS',  # obs date/time
    # keyword; use
    # 'date|time' if
    # separate
    'obsmidtime_jd': 'MIDTIMJD',  # obs midtime jd keyword
    # (usually provided by
    # pp_prepare
    'object': 'OBJECT',  # object name keyword
    'filter': 'FILTER',  # filter keyword
    'filter_translations': {'g': 'g', 'r': 'r',
                            'i': 'i'},
    # filtername translation dictionary
    'exptime': 'EXPTIME',  # exposure time keyword (s)
    'airmass': 'AIRMASS',  # airmass keyword


    # source extractor settings
    'source_minarea': 25,  # default sextractor source minimum N_pixels
    'source_snr': 3,  # default sextractor source snr for registration
    'aprad_default': 7,  # default aperture radius in px
    'aprad_range': [3, 20],  # [minimum, maximum] aperture radius (px)
    'sex-config-file': rootpath+'/setup/mmtcam.sex',
    'mask_file': {},
    #                        mask files as a function of x,y binning

    # registration settings (Scamp)
    'scamp-config-file': rootpath+'/setup/mmtcam.scamp',
    'reg_max_mag': 25,
    'reg_search_radius': 0.1,  # deg
    'source_tolerance': 'high',

    # swarp settings
    'copy_keywords': ('TELESCOP,INSTRUME,FILTER,EXPTIME,OBJECT,' +
                      'DATE-OBS,OBJCTRA,OBJCTDEC,SECPIX,AIRMASS,' +
                      'TEL_KEYW,XBINNING,YBINNING,MIDTIMJD'),
    #                         keywords to be copied in image
    #                         combination using swarp
    'swarp-config-file': rootpath+'/setup/mmtcam.swarp',

    # default catalog settings
    'astrometry_catalogs': ['GAIA'],
    'photometry_catalogs': ['PANSTARRS', 'SDSS-R9', 'APASS9']
}

# Magellan, LDSS3-C
magldss3_param = {
    'telescope_instrument': 'Magellan/LDSS3',  # telescope/instrument name
    'telescope_keyword': 'MAGLDSS3',      # telescope/instrument keyword
    'observatory_code': '268',         # MPC observatory code
    'secpix': (0.189, 0.189),  # pixel size (arcsec)
    # before binning
    'ext_coeff': 0.05,          # typical extinction coefficient


    # image orientation preferences
    'flipx': True,
    'flipy': True,
    'rotate': 0,

    # instrument-specific FITS header keywords
    'binning': ('BINNING#x1', 'BINNING#x2'),  # binning in x/y
    'extent': ('NAXIS1', 'NAXIS2'),   # N_pixels in x/y
    'ra': 'RA',  # telescope pointing, RA
    'dec': 'DEC',  # telescope pointin, Dec
    'radec_separator': ':',   # RA/Dec hms separator, use 'XXX'
    # if already in degrees
    'date_keyword': 'DATE-OBS|TIME-OBS',  # obs date/time
    # keyword; use
    # 'date|time' if
    # separate
    'obsmidtime_jd': 'MIDTIMJD',  # obs midtime jd keyword
    # (usually provided by
    # pp_prepare
    'object': 'OBJECT',  # object name keyword
    'filter': 'FILTER',  # filter keyword
    'filter_translations': {'W4800_7800': 'VR'},
    # filtername translation dictionary
    'exptime': 'EXPTIME',  # exposure time keyword (s)
    'airmass': 'AIRMASS',  # airmass keyword


    # source extractor settings
    'source_minarea': 10,  # default sextractor source minimum N_pixels
    'source_snr': 3,  # default sextractor source snr for registration
    'aprad_default': 3,  # default aperture radius in px
    'aprad_range': [3, 15],  # [minimum, maximum] aperture radius (px)
    'sex-config-file': rootpath+'/setup/magldss3.sex',
    'mask_file': {},
    #                        mask files as a function of x,y binning

    # registration settings (Scamp)
    'scamp-config-file': rootpath+'/setup/magldss3.scamp',
    'reg_max_mag': 25,
    'reg_search_radius': 0.1,  # deg
    'source_tolerance': 'high',

    # swarp settings
    'copy_keywords': ('TELESCOP,INSTRUME,FILTER,EXPTIME,OBJECT,' +
                      'DATE-OBS,TIME-OBS,RA,DEC,SECPIX,AIRMASS,' +
                      'TEL_KEYW,BINNING,MIDTIMJD'),
    #                         keywords to be copied in image
    #                         combination using swarp
    'swarp-config-file': rootpath+'/setup/mmtcam.swarp',

    # default catalog settings
    'astrometry_catalogs': ['GAIA'],
    'photometry_catalogs': ['PANSTARRS', 'SDSS-R9', 'APASS9']
}

# Steward 90", SCC
steward90scc_param = {
    'telescope_instrument': 'Steward90/SCC',  # telescope/instrument name
    'telescope_keyword': 'STEWARD90SCC',  # telescope/instrument keyword
    'observatory_code': '695',         # MPC observatory code
    'secpix': (0.1455, 0.1455),  # pixel size (arcsec)
    # before binning
    'ext_coeff': 0.05,          # typical extinction coefficient


    # image orientation preferences
    'flipx': True,
    'flipy': True,
    'rotate': 0,

    # instrument-specific FITS header keywords
    'binning': ('XBINNING', 'YBINNING'),
    # binning in x/y, '_blankN' denotes that both axes
    # are listed in one keyword, sep. by blanks
    'extent': ('NAXIS1', 'NAXIS2'),   # N_pixels in x/y
    'ra': 'RA',  # telescope pointing, RA
    'dec': 'DEC',  # telescope pointin, Dec
    'radec_separator': ':',   # RA/Dec hms separator, use 'XXX'
    # if already in degrees
    'date_keyword': 'DATE-OBS',  # obs date/time
    # keyword; use
    # 'date|time' if
    # separate
    'obsmidtime_jd': 'MIDTIMJD',  # obs midtime jd keyword
    # (usually provided by
    # pp_prepare
    'object': 'OBJECT',  # object name keyword
    'filter': 'FILTER',  # filter keyword
    'filter_translations': {'r': 'r'},
    # filtername translation dictionary
    'exptime': 'EXPTIME',  # exposure time keyword (s)
    'airmass': 'AIRMASS',  # airmass keyword


    # source extractor settings
    'source_minarea': 9,  # default sextractor source minimum N_pixels
    'source_snr': 3,  # default sextractor source snr for registration
    'aprad_default': 4,  # default aperture radius in px
    'aprad_range': [2, 10],  # [minimum, maximum] aperture radius (px)
    'sex-config-file': rootpath+'/setup/steward90scc.sex',
    'mask_file': {},
    #                        mask files as a function of x,y binning

    # registration settings (Scamp)
    'scamp-config-file': rootpath+'/setup/steward90scc.scamp',
    'reg_max_mag': 19,
    'reg_search_radius': 0.5,  # deg
    'source_tolerance': 'high',

    # swarp settings
    'copy_keywords': ('OBSERVAT,INSTRUME,EXPTIME,OBJECT,' +
                      'DATE-OBS,RA,DEC,AIRMASS,TEL_KEYW,CCDSUM,' +
                      'FILTERS,MIDTIMJD'),
    #                        keywords to be copied in image
    #                        combination using swarp
    'swarp-config-file': rootpath+'/setup/steward90scc.swarp',

    # default catalog settings
    'astrometry_catalogs': ['GAIA'],
    'photometry_catalogs': ['PANSTARRS', 'SDSS-R9', 'APASS9']
}

# McDonald 2.1m Struve, CQUEAN
struvecquean_param = {
    'telescope_instrument': 'Struve/CQUEAN',  # telescope/instrument name
    'telescope_keyword': 'stuvezquean',      # telescope/instrument keyword
    'observatory_code': '711',         # MPC observatory code
    'secpix': (0.276, 0.276),  # pixel size (arcsec)
    # before binning
    'ext_coeff': 0.05,          # typical extinction coefficient


    # image orientation preferences
    'flipx': True,
    'flipy': False,
    'rotate': 0,

    # instrument-specific FITS header keywords
    'binning': (1, 1),  # binning in x/y
    'extent': ('NAXIS1', 'NAXIS2'),   # N_pixels in x/y
    'ra': 'TEL_RA',  # telescope pointing, RA
    'dec': 'TEL_DEC',  # telescope pointin, Dec
    'radec_separator': ':',   # RA/Dec hms separator, use 'XXX'
    # if already in degrees
    'date_keyword': 'DATE-OBS|TIME-OBS',  # obs date/time
    # keyword; use
    # 'date|time' if
    # separate
    'obsmidtime_jd': 'MIDTIMJD',  # obs midtime jd keyword
    # (usually provided by
    # pp_prepare
    'object': 'OBJECT',  # object name keyword
    'filter': 'FILTER',  # filter keyword
    'filter_translations': {'B': 'B', 'V': 'V',
                            'R': 'R', 'I': 'I'},
    # filtername translation dictionary
    'exptime': 'EXPTIME',  # exposure time keyword (s)
    'airmass': 'AIRMASS',  # airmass keyword


    # source extractor settings
    'source_minarea': 12,  # default sextractor source minimum N_pixels
    'source_snr': 3,  # default sextractor source snr for registration
    'aprad_default': 5,  # default aperture radius in px
    'aprad_range': [2, 10],  # [minimum, maximum] aperture radius (px)
    'sex-config-file': rootpath+'/setup/vatt4k.sex',
    'mask_file': {},
    #                        mask files as a function of x,y binning

    # registration settings (Scamp)
    'scamp-config-file': rootpath+'/setup/vatt4k.scamp',
    'reg_max_mag': 19,
    'reg_search_radius': 0.5,  # deg
    'source_tolerance': 'high',

    # swarp settings
    'copy_keywords': ('TELESCOP,FILTER,EXPTIME,OBJECT,' +
                      'DATE-OBS,TIME-OBS,TEL_RA,TEL_DEC,SECPIX,AIRMASS,' +
                      'TEL_KEYW,MIDTIMJD'),
    #                         keywords to be copied in image
    #                         combination using swarp
    'swarp-config-file': rootpath+'/setup/vatt4k.swarp',

    # default catalog settings
    'astrometry_catalogs': ['GAIA'],
    'photometry_catalogs': ['PANSTARRS', 'SDSS-R9', 'APASS9']
}

ztfmosaic_param = {
    'telescope_instrument': 'ZTF/MOSAIC',  # telescope/instrument name
    'telescope_keyword': 'ZTFMOSAIC',  # telescope/instrument keyword
    'observatory_code': 'I41',  # MPC observatory code
    'secpix': (1.012, 1.012),  # pixel size (arcsec) before binning

    # image orientation preferences
    'flipx': True,
    'flipy': False,
    'rotate': 0,

    # instrument-specific FITS header keywords
    'binning': ('CCDSUM#blank0', 'CCDSUM#blank1'),  # binning in x/y
    'extent': ('NAXIS1', 'NAXIS2'),  # N_pixels in x/y
    'ra': 'OBJRAD',  # telescope pointing, RA
    'dec': 'OBJDECD',  # telescope pointin, Dec
    'radec_separator': 'XXX',  # RA/Dec hms separator, use 'XXX'
    # if already in degrees
    'date_keyword': 'OBSJD',  # obs date/time
    # keyword; use
    # 'date|time' if
    # separate
    'obsmidtime_jd': 'MIDTIMJD',  # obs midtime jd keyword
    # (usually provided by
    # pp_prepare
    'object': 'OBJECT',  # object name keyword
    'filter': 'FILTER',  # filter keyword
    'filter_translations': {'ZTF_g': 'g', 'ZTF_r': 'r'},
    # filtername translation dictionary
    'exptime': 'EXPTIME',  # exposure time keyword (s)
    'airmass': 'AIRMASS',  # airmass keyword

    # source extractor settings
    'source_minarea': 12,  # default sextractor source minimum N_pixels
    'source_snr': 3,  # default sextractor source snr for registration
    'aprad_default': 5,  # default aperture radius in px
    'aprad_range': [2, 10],  # [minimum, maximum] aperture radius (px)
    'sex-config-file': rootpath + '/setup/ztfmosaic.sex',
    'mask_file': {},
    #                        mask files as a function of x,y binning

    # scamp settings
    'scamp-config-file': rootpath + '/setup/ztfmosaic.scamp',
    'reg_max_mag': 19,
    'reg_search_radius': 0.1,  # deg
    'source_tolerance': 'high',

    # default catalog settings
    'astrometry_catalogs': ['GAIA'],
    'photometry_catalogs': ['PANSTARRS', 'SDSS-R9']
}

# INT, WFC
intwfc_param = {
    'telescope_instrument': 'INT/WFC',  # telescope/instrument name
    'telescope_keyword': 'INTWFC',      # telescope/instrument keyword
    'observatory_code': '950',         # MPC observatory code
    'secpix': (0.333, 0.333),  # pixel size (arcsec)
    # before binning
    'ext_coeff': 0.05,          # typical extinction coefficient


    # image orientation preferences
    'flipx': True,
    'flipy': False,
    'rotate': 90,

    # instrument-specific FITS header keywords
    'binning': ('CCDXBIN', 'CCDYBIN'),  # binning in x/y
    'extent': ('NAXIS1', 'NAXIS2'),   # N_pixels in x/y
    'ra': 'RA',  # telescope pointing, RA
    'dec': 'DEC',  # telescope pointin, Dec
    'radec_separator': ':',   # RA/Dec hms separator, use 'XXX'
    # if already in degrees
    'date_keyword': 'DATE-OBS|UT',  # obs date/time
    # keyword; use
    # 'date|time' if
    # separate
    'obsmidtime_jd': 'MIDTIMJD',  # obs midtime jd keyword
    # (usually provided by
    # pp_prepare
    'object': 'OBJECT',  # object name keyword
    'filter': 'WFFBAND',  # filter keyword
    'filter_translations': {'V': 'V'},
    # filtername translation dictionary
    'exptime': 'EXPTIME',  # exposure time keyword (s)
    'airmass': 'AIRMASS',  # airmass keyword


    # source extractor settings
    'source_minarea': 12,  # default sextractor source minimum N_pixels
    'source_snr': 3,  # default sextractor source snr for registration
    'aprad_default': 5,  # default aperture radius in px
    'aprad_range': [2, 12],  # [minimum, maximum] aperture radius (px)
    'sex-config-file': rootpath+'/setup/intwfc.sex',
    'mask_file': {},
    #                        mask files as a function of x,y binning

    # registration settings (Scamp)
    'scamp-config-file': rootpath+'/setup/intwfc.scamp',
    'reg_max_mag': 19,
    'reg_search_radius': 0.5,  # deg
    'source_tolerance': 'high',

    # swarp settings
    'copy_keywords': ('TELESCOP,INSTRUME,FILTER,EXPTIME,OBJECT,' +
                      'DATE-OBS,UT,RA,DEC,SECPIX,AIRMASS,' +
                      'TEL_KEYW,CCDXBIN,CCDYBIN,MIDTIMJD'),
    #                         keywords to be copied in image
    #                         combination using swarp
    'swarp-config-file': rootpath+'/setup/intwfc.swarp',

    # default catalog settings
    'astrometry_catalogs': ['GAIA'],
    'photometry_catalogs': ['PANSTARRS', 'SDSS-R9', 'APASS9']
}

# NOT, ALFOSC
notalfosc_param = {
    'telescope_instrument': 'NOT/ALFOSC',  # telescope/instrument name
    'telescope_keyword': 'NOTALFOSC',      # telescope/instrument keyword
    'observatory_code': 'Z18',         # MPC observatory code
    'secpix': (0.2138, 0.2138),  # pixel size (arcsec)
    # before binning
    'ext_coeff': 0.05,          # typical extinction coefficient


    # image orientation preferences
    'flipx': True,
    'flipy': False,
    'rotate': 0,

    # instrument-specific FITS header keywords
    'binning': ('DETXBIN', 'DETYBIN'),  # binning in x/y
    'extent': ('NAXIS1', 'NAXIS2'),   # N_pixels in x/y
    'ra': 'RA',  # telescope pointing, RA
    'dec': 'DEC',  # telescope pointin, Dec
    'radec_separator': 'XXX',   # RA/Dec hms separator, use 'XXX'
    # if already in degrees
    'date_keyword': 'DATE-OBS',  # obs date/time
    # keyword; use
    # 'date|time' if
    # separate
    'obsmidtime_jd': 'MIDTIMJD',  # obs midtime jd keyword
    # (usually provided by
    # pp_prepare
    'object': 'OBJECT',  # object name keyword
    'filter': 'ALFLTNM',  # filter keyword
    'filter_translations': {'Open': None},
    # filtername translation dictionary
    'exptime': 'EXPTIME',  # exposure time keyword (s)
    'airmass': 'AIRMASS',  # airmass keyword


    # source extractor settings
    'source_minarea': 12,  # default sextractor source minimum N_pixels
    'source_snr': 3,  # default sextractor source snr for registration
    'aprad_default': 5,  # default aperture radius in px
    'aprad_range': [2, 10],  # [minimum, maximum] aperture radius (px)
    'sex-config-file': rootpath+'/setup/notalfosc.sex',
    'mask_file': {},
    #                        mask files as a function of x,y binning

    # registration settings (Scamp)
    'scamp-config-file': rootpath+'/setup/notalfosc.scamp',
    'reg_max_mag': 19,
    'reg_search_radius': 0.5,  # deg
    'source_tolerance': 'high',

    # swarp settings
    'copy_keywords': ('TELESCOP,INSTRUME,ALFTLNM,EXPTIME,OBJECT,' +
                      'DATE-OBS,RA,DEC,SECPIX,AIRMASS,' +
                      'TEL_KEYW,DETXBIN,DETYBIN,MIDTIMJD'),
    #                         keywords to be copied in image
    #                         combination using swarp
    'swarp-config-file': rootpath+'/setup/notalfosc.swarp',

    # default catalog settings
    'astrometry_catalogs': ['GAIA'],
    'photometry_catalogs': ['PANSTARRS', 'SDSS-R9', 'APASS9']
}

# NEXT
nextfli_param = {
    'telescope_instrument': 'NEXT',  # telescope/instrument name
    'telescope_keyword': 'NEXT',  # telescope/instrument keyword
    'observatory_code': 'C42',  # MPC observatory code
    'secpix': (0.64, 0.64),  # pixel size (arcsec) before binning

    # image orientation preferences
    'flipx': False,
    'flipy': False,
    'rotate': 0,

    # instrument-specific FITS header keywords
    'binning': ('XBINNING', 'YBINNING'),  # binning in x/y
    'extent': ('NAXIS1', 'NAXIS2'),  # N_pixels in x/y
    'ra': 'RA',  # telescope pointing, RA
    'dec': 'DEC',  # telescope pointin, Dec
    'radec_separator': ' ',  # RA/Dec hms separator, use 'XXX'
    # if already in degrees
    'date_keyword': 'DATE-OBS',  # obs date/time
    # keyword; use
    # 'date|time' if
    # separate
    'obsmidtime_jd': 'MJD-OBS',  # obs midtime jd keyword
    # (usually provided by
    # pp_prepare
    'object': 'OBJECT',  # object name keyword
    'filter': 'FILTER',  # filter keyword
    'filter_translations': {'V': 'V', 'R': 'R', 'B': 'B', 'g2': 'g', 'r2': 'r', 'i2': 'i'},
    # filtername translation dictionary
    'exptime': 'EXPTIME',  # exposure time keyword (s)
    'airmass': 'AIRMASS',  # airmass keyword

    # source extractor settings
    'source_minarea': 12,  # default sextractor source minimum N_pixels
    'source_snr': 3,  # default sextractor source snr for registration
    'aprad_default': 5,  # default aperture radius in px
    'aprad_range': [2, 10],  # [minimum, maximum] aperture radius (px)
    'sex-config-file': rootpath + '/setup/next.sex',
    'mask_file': {},
    #                        mask files as a function of x,y binning

    # scamp settings
    'scamp-config-file': rootpath + '/setup/next.scamp',
    'reg_max_mag': 19,
    'reg_search_radius': 0.5,  # deg
    'source_tolerance': 'high',
    
    # swarp settings
    'copy_keywords': ('OBSERVAT,INSTRUME,EXPTIME,OBJECT,' +
                      'DATE-OBS,TEL_KEYW,XBINNING,YBINNING,' +
                      'FILTER,RA,DEC'),
    #                        keywords to be copied in image
    #                        combination using swarp
    'swarp-config-file': rootpath+'/setup/next.swarp',

    # default catalog settings
    'astrometry_catalogs': ['GAIA'],
    'photometry_catalogs': ['PANSTARRS']
}


# access functions for telescope configurations


implemented_telescopes = ['VATT4K', 'DCTLMI', 'ARC35ARCTIC',
                          'ARC35AGILE', 'MAGIMACSL', 'MAGIMACSS',
                          'LOWELL31', 'LOWELL42',
                          'LOWELL72',
                          'CTIO09', 'CTIO10', 'CTIO13CCD', 'UH88SNIFS',
                          'WIYN09HDI', 'RATIR', 'SOARGOODMANold', 'SOARGOODMAN',
                          'OHP120',
                          'TNGDOLORES', 'GENERIC', 'KPNO4MOS1', 'FROST',
                          'MEXMAN', 'KPNO4MOS1', 'KPNOMOS3',
                          'KPNO4NEWF', 'UKIRTWFCAM', 'VLTFORS2',
                          'LOWELL42SITE', 'PLUTOPLATE', 'TCS15MUSCAT2',
                          'LCOSBIGKB78', 'ARC35SPICAM', 'LCOSINFL03',
                          'LCOSINFL06',
                          'LCOSINFL16', 'LCOSINFL11',
                          'LCOSINFA03', 'LCOSINFA15',
                          'LCOSPECFS01', 'P60OPT', 'P60SEDM', 'GMOSN',
                          'DFOSC', 'LONEOS', 'PDO25CMF63ST8', 'PDO05F81KAF1001E',
                          'PDS35CMSTL1001E', 'MMTCAM', 'MAGLDSS3',
                          'SL40IN', 'STEWARD90SCC', 'STRUVECQUEAN',
                          'ZTFMOSAIC', 'NOTALFOSC', 'NEXT']

# translate INSTRUME (or others, see _pp_conf.py) header keyword into
# PP telescope keyword
instrument_identifiers = {'= "Vatt4k"':        'VATT4K',
                          'LMI':               'DCTLMI',
                          'lmi': 'DCTLMI',
                          'arctic':            'ARC35ARCTIC',
                          'agile':             'ARC35AGILE',
                          'IMACS Long-Camera': 'MAGIMACSL',
                          'IMACS Short-Camera': 'MAGIMACSS',
                          'DLR-MKIII':         'CA123DLRMKIII',
                          'NASAcam':           'LOWELL31',
                          'nasa42':            'LOWELL42',
                          'E2V CCD-231 4096x4112': 'LOWELL42',
                          'PRISM Instrument':  'LOWELL72',
                          'Prism 2048x2048 CCD': 'LOWELL72',
                          'prism': 'LOWELL72',
                          'cfccd':             'CTIO09',
                          'Y4KCam':            'CTIO10',
                          'ANDICAM-CCD':       'CTIO13CCD',
                          'SNIFS':             'UH88SNIFS',
                          'hdi':               'WIYN09HDI',
                          'ArtemisHSC':        'GENERIC',
                          'GENERIC':           'GENERIC',
                          'C0':                'RATIR',
                          'C1':                'RATIR',
                          'C2':                'RATIR',
                          'C3':                'RATIR',
                          'C4':                'RATIR',
                          'SHA':               'SL40IN',
                          'Goodman Spectrograph': 'SOARGOODMANOLD',
                          'Andor Tech':        'OHP120',
                          'LRS':               'TNGDOLORES',
                          'mosaic_1_1':        'KPNO4MOS1',
                          'mosaic_1':          'KPNO4MOS1',
                          'KMTS':              'KMTNETS',
                          'SI 600-277': 'FROST',
                          'Mexman': 'MEXMAN',
                          'mosaic_1_1':        'KPNO4MOS1',
                          'mosaic_1':          'KPNO4MOS1',
                          'KMTS':              'KMTNETS',
                          'newfirm': 'KPNO4NEWF',
                          'Mosaic3': 'KPNO4MOS3',
                          'WFCAM': 'UKIRTWFCAM',
                          'SIRIUS': 'IRSFSIRIUS',
                          'Goodman Spectro': 'SOARGOODMAN',
                          'FORS2': 'VLTFORS2',
                          '2:1 f/17 direct': 'LOWELL42SITE',
                          'Pluto plate': 'PLUTOPLATE',
                          'MuSCAT2': 'TCS15MUSCAT2',
                          'kb78': 'LCOSBIGKB78',
                          'spicam': 'ARC35SPICAM',
                          'fl03': 'LCOSINFL03',
                          'fl06': 'LCOSINFL06',
                          'fl16': 'LCOSINFL16',
                          'fl11': 'LCOSINFL11',
                          'fa03': 'LCOSINFA03',
                          'fa15': 'LCOSINFA15',
                          'fs01': 'LCOSPECFS01',
                          'P60': 'P60OPT',
                          'Rainbow Cam': 'P60SEDM',
                          'GMOS-N': 'GMOSN',
                          'DFOSC_FASU': 'DFOSC',
                          'loneos': 'LONEOS',
                          '25cm f/6.3 SCT_SBIG ST-8': 'PDO25CMF63ST8',
                          '0.5m f/8.1 Ritchey-Chretien_FLI KAF1001E':
                          'PDO05F81KAF1001E',
                          '0.35-m SCT_STL-1001E': 'PDS35CMSTL1001E',
                          'MMT Rapid Imager': 'MMTCAM',
                          'LDSS3-C': 'MAGLDSS3',
                          'Finger Lakes Instr. ProLine Model PL23042, S/N PL0101015':
                              'STEWARD90SCC',
                          '2.1m Otto Struve': 'STRUVECQUEAN',
                          'ZTF/MOSAIC': 'ZTFMOSAIC',
                          'WFC': 'INTWFC',
                          'ALFOSC_FASU': 'NOTALFOSC',
                          'FLI': 'NEXT'}

# translate telescope keyword into parameter set defined here
telescope_parameters = {'VATT4K':       vatt4k_param,
                        'DCTLMI':        dctlmi_param,
                        'ARC35ARCTIC':   arc35arctic_param,
                        'ARC35AGILE':    arc35agile_param,
                        'MAGIMACSL':      magimacsl_param,
                        'MAGIMACSS':      magimacss_param,
                        'CA123DLRMKIII': ca123dlrmkiii_param,
                        'LOWELL31':      lowell31_param,
                        'LOWELL42':      lowell42_param,
                        'LOWELL72':      lowell72_param,
                        'CTIO09':        ctio09_param,
                        'CTIO10':        ctio10_param,
                        'CTIO13CCD':     ctio13ccd_param,
                        'UH88SNIFS':     uh88snifs_param,
                        'WIYN09HDI':     wiyn09hdi_param,
                        'GENERIC':       generic_param,
                        'RATIR':         ratir_param,
                        'SOARGOODMANOLD': soargoodmanold_param,
                        'SOARGOODMAN': soargoodman_param,
                        'OHP120':        ohp120_param,
                        'TNGDOLORES':    tngdolores_param,
                        'KPNO4MOS1':     kpno4mos1_param,
                        'KMTNETS':       kmtnets_param,
                        'FROST':         frost_param,
                        'MEXMAN':        mexman_param,
                        'KPNO4MOS1':     kpno4mos1_param,
                        'KMTNETS':       kmtnets_param,
                        'FROST':         frost_param,
                        'KPNO4MOS3': kpno4mos3_param,
                        'KPNO4NEWF': kpno4newf_param,
                        'UKIRTWFCAM': ukirtwfcam_param,
                        'IRSFSIRIUS': irsfsirius_param,
                        'VLTFORS2': vltfors2_param,
                        'LOWELL42SITE': lowell42site_param,
                        'PLUTOPLATE': plutoplate_param,
                        'TCS15MUSCAT2': tcs15muscat2_param,
                        'LCOSBIGKB78': lcosbigkb78_param,
                        'ARC35SPICAM': arc35spicam_param,
                        'LCOSINFL03': lcosinfl03_param,
                        'LCOSINFL06': lcosinfl06_param,
                        'LCOSINFL16': lcosinfl16_param,
                        'LCOSINFL11': lcosinfl11_param,
                        'LCOSINFA03': lcosinfa03_param,
                        'LCOSINFA15': lcosinfa15_param,
                        'LCOSPECFS01': lcospecfs01_param,
                        'P60OPT': p60opt_param,
                        'P60SEDM': p60sedm_param,
                        'GMOSN': gmosn_param,
                        'DFOSC': dfosc_param,
                        'LONEOS': loneos_param,
                        'PDO25CMF63ST8': pdo25cmf63st8_param,
                        'PDO05F81KAF1001E': pdo05mf81kaf1001e_param,
                        'PDS35CMSTL1001E': pds35cmstl1001e_param,
                        'MMTCAM': mmtcam_param,
                        'MAGLDSS3': magldss3_param,
                        'SL40IN': sl40in_param,
                        'STEWARD90SCC': steward90scc_param,
                        'STRUVECQUEAN': struvecquean_param,
                        'ZTFMOSAIC': ztfmosaic_param,
                        'INTWFC': intwfc_param,
                        'NOTALFOSC': notalfosc_param,
                        'NEXT': nextfli_param}


# append mytelescopes.py, if available
#
# mytelescopes.py allows you to setup your own telescope; that file is
# not part of the github repository, hence it will not be affected by
# pulls and pushes
#
# an example mytelescopes.py file is available here:
# http://134.114.60.45/photometrypipeline/mytelescopes.py
# more information are available on the PP documentation website:
# http://mommermi.github.io/pp/install.html

try:
    execfile(rootpath+'/setup/mytelescopes.py')
except IOError:
    pass

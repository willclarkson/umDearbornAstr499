from numpy import sqrt
""" PP_SETUP - photometry pipeline setup file
    v1.0: 2018-11-15, mommermiscience@gmail.com
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


"""This configuration file is meant to customize your pipeline
according to your needs. This is not an executable file. Instead, each
pp function will import the classes and use the definitions made in
here - parameters provided to the function through the command line
will override the definitions in this class. Definitions below are
group into different class based on their pp function association.
"""


class Conf():
    """General pipeline configurations"""

    diagnostics = True  # produce diagnostic files and website?


class ConfPrepare(Conf):
    """configuration setup for pp_prepare"""
    pass


class ConfRegister(Conf):
    """configuration setup for pp_register"""
    pass


class ConfPhotometry(Conf):
    """configuration setup for pp_photometry"""
    pass


class ConfCalibrate(Conf):
    """Configuration setup for pp_calibrate"""

    # write photometric calibration raw data into file
    save_caldata = True  # save data on calibration process
    save_caldata_format = 'ascii.basic'  # Table.write formats (ascii, csv...)
    save_caldata_suffix = '_cal.dat'  # file suffix ('.fits' will be clipped)
    save_caldata_usedonly = False  # only output stars used in calibration?

    # add photometric calibration raw data to frame database
    caldata_in_db = True  # add calibration data to database file?


class ConfDistill(Conf):
    """configuration setup for pp_distill"""

    # target rejection dictionary using `dat` as defined in pp_distill.distill
    # key refers to rejection scheme identifier
    # read as: schema `key` rejects sources with...
    rejection = {
        # geometric positional uncertainties > 10"
        'pos': lambda dat: (sqrt((dat[1]-dat[3])**2 +
                                 (dat[2]-dat[4])**2)*3600 > 10),
        'none': lambda dat: False,
    }


class ConfDiagnostics(Conf):
    """configuration setup for diagnostics"""

    # general settings
    diagnostics_path = '.'
    # '.' puts diagnostics.html into the data directory
    # ancillary files go into '.diagnostics' relative to diagnostics_path
    # can be any other absolute path; can be changed during runtime
    main_html = 'diagnostics.html'  # main html document
    image_file_format = 'png'  # output format for images and plots
    individual_frame_pages = True  # produce frame-specific pages

    # image settings
    image_stretch = 'linear'  # could be 'linear', 'log', 'power'
    image_size_lg_px = 1000  # cutout size for large image in px
    image_size_lg_in = 5     # image size for large image in inches
    image_size_thumb_px = 200  # cutout size for thumbnail in px
    image_size_thumb_in = 2.5     # image size for thumbnail in inches
    image_dpi = 150  # image output dpi
    overlay_lg_linewidth = 0.5  # linewidth in overlays for large images
    show_quickview_image = True  # show quickview image
    # (the quickview image is required for registration_star_map and
    #  calibration_star_map)

    # plot settings
    plot_dpi = 100  # plot output dpi

    # presentation of registration results
    show_registration_table = True  # show table with registration results
    show_registration_star_map = True  # registration catalog on images

    # presentation of calibration data
    show_phot_calibration_plots = True  # present calibration results plot
    show_calibration_star_map = True  # present map of calibration stars
    show_calibration_star_table = True  # present table of calibration stars

    # distill settings
    show_target_animations = True  # build and show target gif animations

    # target thumbnail overlay properties
    thumb_scalelength = 10  # length of pixelscale indicators in arcsec
    thumb_fontsize = 8  # fontsize for text in thumbnail overlay
    thumb_linewidth = 0.3  # linewidth for lines on overlay
    thumb_predicted_pos_marker = 'x'  # marker style for predicted position
    thumb_predicted_pos_size = 20  # marker size for predicted position
    thumb_predicted_pos_color = 'cornflowerblue'  # marker color


class ConfCombine(Conf):
    """configuration setup for pp_combine"""
    pass


class ConfMPCReport(Conf):
    """configuration setup for pptool_mpcreport"""
    pass


confprepare = ConfPrepare()
confcalibrate = ConfCalibrate()
confdistill = ConfDistill()
confdiagnostics = ConfDiagnostics()

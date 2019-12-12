""" DIAGNOSTICS - diagnostic routines for photometry pipeline
    v1.0: 2016-02-25, mommermiscience@gmail.com
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

import os
import sys
import numpy as np
import logging
import subprocess

from astropy.io import fits
from astropy import wcs
from astropy.visualization import (ZScaleInterval, ImageNormalize,
                                   LogStretch, LinearStretch)
from astropy.time import Time

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pylab as plt
    matplotlib.rcdefaults()  # restore default parameters
except ImportError:
    print('Module matplotlib not found. Please install with: pip install '
          'matplotlib')
    sys.exit()

try:
    from skimage.transform import resize
except ImportError:
    print('Module skimage not found. Please install with: pip install '
          'scikit-kimage')
    sys.exit()

# pipeline-specific modules
import _pp_conf
import toolbox
from catalog import *

# setup logging
logging.basicConfig(filename=_pp_conf.log_filename,
                    level=_pp_conf.log_level,
                    format=_pp_conf.log_formatline,
                    datefmt=_pp_conf.log_datefmt)


class Diagnostics_Html():
    """basis class for building pp html diagnostic output"""

    from pp_setup import confdiagnostics as conf

    def create_website(self, filename, content=''):
        """
        create empty website for diagnostics output
        """

        html = ("<!DOCTYPE html PUBLIC '-//W3C//DTD HTML 4.01//EN'>\n"
                "<HTML>\n"
                "<HEAD>\n"
                "<TITLE>Photometry Pipeline - Diagnostics</TITLE>\n"
                "<LINK rel=\"stylesheet\" href=\"{:s}"
                "diagnostics_stylesheet.css\">\n"
                "</HEAD>\n"
                "<BODY>\n"
                "<script>\n"
                "function toggledisplay(elementID)\n"
                "{{\n"
                "(function(style) {{\n"
                "style.display = style.display === 'none' ? '' :"
                "'none';\n"
                "}})(document.getElementById(elementID).style);\n"
                "}}\n"
                "</script>\n\n"
                "{:s}\n"
                "</BODY>\n"
                "</HTML>\n").format(os.getenv('PHOTPIPEDIR'), content)

        outf = open(filename, 'w')
        outf.writelines(html)
        outf.close()

    def append_website(self, filename, content,
                       replace_from='X?!do not replace anything!?X',
                       keep_at='</BODY>',):
        """append content to an existing website: replace content starting
        at line containing `replace_from` until line containin `keep_at`;
        by default, all content following `replace_from` is
        replaced
        """
        # read existing code
        existing_html = open(filename, 'r').readlines()

        # insert content into existing html
        outf = open(filename, 'w')
        delete = False
        for line in existing_html:
            if replace_from in line:
                delete = True
                continue
            if keep_at in line:
                outf.writelines(content)
                delete = False
            if delete:
                continue
            outf.writelines(line)
        outf.close()

    def abort(self, where):
        """
        use this function to add information to index.html that the
        pipeline crashed and where
        """
        logging.info('adding pipeline crash to diagnostics')

        html = ("<P><FONT COLOR=\"RED\">Pipeline crashed "
                "unexpectedly in module {:s}; refer to <A "
                "HREF=\"{:s}\">log</A> "
                "for additional information</FONT>\n").format(
                    _pp_conf.log_filename, where)

        self.append_website(os.path.join(self.conf.diagnostics_path,
                                         self.conf.main_html), html)
        logging.info('pipeline crash added')


class Prepare_Diagnostics(Diagnostics_Html):
    """diagnostics run as part of pp_prepare"""

    function_tag = "<!-- pp_prepare -->"

    def frame_table(self, filenames, obsparam):

        logging.info('create data summary table')

        if self.conf.individual_frame_pages:
            self.frame_pages(filenames, obsparam)

        # create frame information table
        html = "<P><TABLE CLASS=\"gridtable\">\n"
        html += ("<TR><TH>Idx</TH>"
                 "<TH>Filename</TH>"
                 "<TH>Observation Midtime (UT)</TH>"
                 "<TH>Object Name</TH>"
                 "<TH>Airmass</TH>"
                 "<TH>Exptime (s)</TH>"
                 "<TH>Pixel Size (\")"
                 "<TH>Binning</TH>"
                 "<TH>FoV (')</TH></TR>\n")

        for idx, filename in enumerate(filenames):
            hdulist = fits.open(filename, ignore_missing_end=True)
            header = hdulist[0].header
            binning = toolbox.get_binning(header, obsparam)
            try:
                objectname = header[obsparam['object']]
            except KeyError:
                objectname = 'Unknown Target'

            if self.conf.individual_frame_pages:
                framename = "<A HREF=\"{:s}\">{:s}</A>".format(
                    os.path.join(self.conf.diagnostics_path,
                                 '.diagnostics', filename+'.html'),
                    filename)
                if self.conf.show_quickview_image:
                    self.quickview_image(filename)

                    # update frame page
                    framehtml = ("<!-- Quickview -->\n"
                                 "<A HREF=\"#quickview\" "
                                 "ONCLICK=\"toggledisplay"
                                 "('quickview');\"><H2>Quickview Image</H2>"
                                 "</A>\n"
                                 "<IMG ID=\"quickview\" SRC=\"{:s}\" "
                                 "STYLE=\"display: none\"\>\n\n").format(
                                     filename+'.' +
                                     self.conf.image_file_format)
                    self.append_website(
                        os.path.join(self.conf.diagnostics_path,
                                     '.diagnostics',
                                     '{:s}.html'.format(filename)),
                        framehtml, replace_from='<!-- Quickview -->')
            else:
                framename = filename

            html += ("<TR><TD>{:d}</TD>"
                     "<TD>{:s}</TD>"
                     "<TD>{:s}</TD>"
                     "<TD>{:s}</TD>"
                     "<TD>{:4.2f}</TD>"
                     "<TD>{:.1f}</TD>"
                     "<TD>{:.2f} x {:.2f}</TD>"
                     "<TD>{:d} x {:d}</TD>"
                     "<TD>{:.1f} x {:.1f}</TD>\n"
                     "</TR>\n").format(
                         idx+1, framename,
                         Time(header["MIDTIMJD"], format='jd').iso,
                         str(objectname),
                         float(header[obsparam['airmass']]),
                         float(header[obsparam['exptime']]),
                         obsparam['secpix'][0],
                         obsparam['secpix'][1],
                         int(binning[0]), int(binning[1]),
                         float(header[obsparam['extent'][0]]) *
                         obsparam['secpix'][0]*binning[0]/60.,
                         float(header[obsparam['extent'][1]]) *
                         obsparam['secpix'][1]*binning[1]/60.)

        html += '</TABLE>\n'

        logging.info('data summary table created')

        return html

    def quickview_image(self, filename):
        """create quickview image for one frame"""

        logging.info('create image preview for file {:s}'.format(
            filename))

        hdulist = fits.open(filename, ignore_missing_end=True)

        # create frame image
        imgdat = hdulist[0].data.astype(np.float64)

        # normalize imgdat to pixel values 0 < px < 1
        imgdat[np.where(np.isnan(imgdat))[0]] = np.nanmedian(imgdat)
        imgdat = np.clip(imgdat, np.percentile(imgdat, 1),
                         np.percentile(imgdat, 99))
        imgdat = (imgdat-np.min(imgdat)) / np.max(imgdat-np.min(imgdat)+0.1)

        # resize image larger than lg_image_size_px on one side
        imgdat = resize(imgdat,
                        (min(imgdat.shape[0], self.conf.image_size_lg_px),
                         min(imgdat.shape[1], self.conf.image_size_lg_px)))

        plt.figure(figsize=(self.conf.image_size_lg_in,
                            self.conf.image_size_lg_in))

        norm = ImageNormalize(
            imgdat, interval=ZScaleInterval(),
            stretch={'linear': LinearStretch(),
                     'log': LogStretch()}[self.conf.image_stretch])

        img = plt.imshow(imgdat, cmap='gray', norm=norm,
                         origin='lower')
        # remove axes
        plt.axis('off')
        img.axes.get_xaxis().set_visible(False)
        img.axes.get_yaxis().set_visible(False)

        framefilename = os.path.join(self.conf.diagnostics_path,
                                     '.diagnostics', filename + '.' +
                                     self.conf.image_file_format)
        plt.savefig(framefilename, format=self.conf.image_file_format,
                    bbox_inches='tight',
                    pad_inches=0, dpi=self.conf.image_dpi)
        logging.info('image preview for file {:s} written to {:s}'.format(
            filename, os.path.join(self.conf.diagnostics_path,
                                   '.diagnostics', filename + '.' +
                                   self.conf.image_file_format)))

        plt.close()
        hdulist.close()

    def frame_pages(self, filenames, obsparam):
        """build information page for each individual frame"""

        logging.info('setting up individual frame diagnostics report pages')

        for idx, filename in enumerate(filenames):
            header = fits.open(filename)[0].header
            html = ("<H1>{:s} Diagnostics</H1>"
                    "<P><TABLE CLASS=\"gridtable\">\n"
                    "<TR><TH>Telescope/Instrument</TH><TD>{:s} ({:s})</TD>"
                    "</TR>\n"
                    "<TR><TH>Target/Field Identifier</TH><TD>{:s}</TD>"
                    "</TR>\n"
                    "<TR><TH>RA</TH><TD>{:s}</TD></TR>\n"
                    "<TR><TH>Dec</TH><TD>{:s}</TD></TR>\n"
                    "<TR><TH>Exposure Time (s)</TH><TD>{:s}</TD></TR>\n"
                    "<TR><TH>Observation Midtime</TH><TD>{:s}</TD></TR>\n"
                    "</TABLE><P>\n"
                    "<A HREF=\"{:s}\">"
                    "&laquo; previous frame &laquo;</A> | "
                    "<A HREF=\"../diagnostics.html\">run overview</A> | "
                    "<A HREF=\"{:s}\">"
                    "&raquo; next frame &raquo;</A></P>\n\n").format(
                        filename,
                        obsparam['telescope_instrument'],
                        obsparam['telescope_keyword'],
                        header[obsparam['object']],
                        str(header[obsparam['ra']]),
                        str(header[obsparam['dec']]),
                        str(header[obsparam['exptime']]),
                        str(Time(header['MIDTIMJD'], format='jd').iso),
                        filenames[(idx-1) % len(filenames)]+'.html',
                        filenames[(idx+1) % len(filenames)]+'.html')

            self.create_website(
                os.path.join(self.conf.diagnostics_path,
                             '.diagnostics', '{:s}.html'.format(filename)),
                html)
            logging.info(('diagnostics report page for file {:s} '
                          'written to {:s}').format(
                              filename,
                              os.path.join(self.conf.diagnostics_path,
                                           '.diagnostics',
                                           '{:s}.html'.format(filename))))

    def add_index(self, filenames, datadirectory, obsparam):
        """
        create index.html
        diagnostic root website
        """
        logging.info('create frame table')

        os.mkdir(self.conf.diagnostics_path) if not os.path.exists(
            self.conf.diagnostics_path) else None
        os.mkdir(os.path.join(self.conf.diagnostics_path,
                              '.diagnostics')) if not os.path.exists(
                                  os.path.join(self.conf.diagnostics_path,
                                               '.diagnostics')) else None

        # create header information
        refheader = fits.open(filenames[0],
                              ignore_missing_end=True)[0].header
        raw_filtername = refheader[obsparam['filter']]
        translated_filtername = obsparam['filter_translations'][
            refheader[obsparam['filter']]]

        html = ("{:s}\n<H1>Photometry Pipeline Diagnostic Output</H1>\n"
                "<TABLE CLASS=\"gridtable\">\n"
                "  <TR><TH>Data Directory</TH><TD>{:s}</TD></TR>\n"
                "  <TR><TH>Telescope/Instrument</TH><TD>{:s}</TD></TR>\n"
                "  <TR><TH>Number of Frames</TH><TD>{:d}</TD></TR>\n"
                "  <TR><TH>Raw Filter Identifier</TH><TD>{:s}</TD></TR>\n"
                "  <TR><TH>Translated Filter Identifier</TH>"
                "<TD>{:s}</TD></TR>\n"
                "  <TR><TH>Log File</TH>"
                "      <TD><A HREF=\"{:s}\">available here</A></TD></TR>"
                "</TABLE>\n").format(
                    self.function_tag,
                    datadirectory,
                    obsparam['telescope_instrument'],
                    len(filenames),
                    str(raw_filtername),
                    str(translated_filtername),
                    os.path.join(datadirectory, 'LOG'))

        html += "<H3>Data Summary</H3>\n"
        html += self.frame_table(filenames, obsparam)

        self.create_website(os.path.join(self.conf.diagnostics_path,
                                         self.conf.main_html), html)
        logging.info('frame table created')


# registration results website

class Registration_Diagnostics(Diagnostics_Html):

    function_tag = "<!-- pp_register -->"

    def registration_table(self, data, extraction_data, obsparam):
        """build overview table with astrometric registration results"""

        logging.info('creating image registration overview table')

        html = ("<TABLE CLASS=\"gridtable\">\n<TR>\n"
                "<TH>Filename</TH><TH>C<SUB>AS</SUB></TH>"
                "<TH>C<SUB>XY</SUB></TH>"
                "<TH>&sigma;<SUB>RA</SUB> (arcsec)</TH>"
                "<TH>&sigma;<SUB>DEC</SUB> (arcsec)</TH>"
                "<TH>&chi;<SUP>2</SUP><SUB>Reference</SUB></TH>"
                "<TH>&chi;<SUP>2</SUP><SUB>Internal</SUB></TH>\n</TR>\n")

        for dat in data['fitresults']:
            framefilename = os.path.join(self.conf.diagnostics_path,
                                         '.diagnostics',
                                         '{:s}.html'.format(dat[0]))
            if self.conf.individual_frame_pages:
                filename = '<A HREF=\"{:s}\">{:s}</A>'.format(
                    framefilename, dat[0])
            else:
                filename = dat[0]

            html += ("<TR><TD>{:s}</TD>"
                     + "<TD>{:4.1f}</TD><TD>{:4.1f}</TD>"
                     + "<TD>{:5.3f}</TD><TD>{:5.3f}</TD>"
                     + "<TD>{:e}</TD><TD>{:e}</TD>\n</TR>\n").format(
                         filename, dat[1], dat[2], dat[3],
                         dat[4], dat[5], dat[6])
        html += "</TABLE>\n"
        html += ("<P CLASS=\"caption\"><STRONG>Legend</STRONG>: "
                 "C<SUB>AS</SUB>: position "
                 "angle/scale contrast (values >{:.1f} are ok); ").format(
                     _pp_conf.scamp_as_contrast_limit)
        html += ("C<SUB>XY</SUB>: xy-shift contrast "
                 "(values >{:.1f} are ok); ").format(
            _pp_conf.scamp_xy_contrast_limit)
        html += ("&sigma;<SUB>RA</SUB> and &sigma;<SUB>DEC</SUB> "
                 "refer to the internal astrometric uncertainties as "
                 "provided by SCAMP; &chi;<SUP>2</SUP><SUB>Reference</SUB> "
                 "and &chi;<SUP>2</SUP><SUB>Internal</SUB> refer to the "
                 "&chi;<SUP>2</SUP> statistics based on the reference "
                 "catalog and the respective frame as provided by SCAMP."
                 "</P>\n")

        logging.info('image registration overview table created')
        return html

    def registration_maps(self, data, extraction_data, obsparam):
        """build overlays for image maps indicating astrometric reference
        stars"""

        logging.info('create registration overlays with reference stars')

        # load reference catalog
        refcat = catalog(data['catalog'])
        for filename in os.listdir('.'):
            if data['catalog'] in filename and '.cat' in filename:
                refcat.read_ldac(filename)
                break

        # create overlays
        for dat in extraction_data:
            framefilename = os.path.join(self.conf.diagnostics_path,
                                         '.diagnostics',
                                         '{:s}_astrometry.{:s}'.format(
                                             dat['fits_filename'],
                                             self.conf.image_file_format))
            imgdat = fits.open(dat['fits_filename'],
                               ignore_missing_end=True)[0].data
            resize_factor = min(
                1.,
                self.conf.image_size_lg_px/np.max(imgdat.shape))

            header = fits.open(dat['fits_filename'],
                               ignore_missing_end=True)[0].header

            # turn relevant header keys into floats
            # astropy.io.fits bug
            for key, val in list(header.items()):
                if 'CD1_' in key or 'CD2_' in key or \
                   'CRVAL' in key or 'CRPIX' in key or \
                   'EQUINOX' in key:
                    header[key] = float(val)

            plt.figure(figsize=(self.conf.image_size_lg_in,
                                self.conf.image_size_lg_in))
            # create fake image to ensure image dimensions and margins
            img = plt.imshow(np.ones((self.conf.image_size_lg_px,
                                      self.conf.image_size_lg_px))*np.nan,
                             origin='lower')

            # remove axes
            plt.axis('off')
            img.axes.get_xaxis().set_visible(False)
            img.axes.get_yaxis().set_visible(False)

            # plot reference sources
            if refcat.shape[0] > 0:
                try:
                    w = wcs.WCS(header)
                    world_coo = np.array(list(zip(refcat['ra_deg'],
                                                  refcat['dec_deg'])))
                    img_coo = w.wcs_world2pix(world_coo, True)
                    img_coo = [c for c in img_coo
                               if (c[0] > 0 and c[1] > 0 and
                                   c[0] < header[obsparam['extent'][0]] and
                                   c[1] < header[obsparam['extent'][1]])]
                    plt.scatter([c[0]*resize_factor for c in img_coo],
                                [c[1]*resize_factor for c in img_coo],
                                s=5, marker='o', edgecolors='red',
                                linewidth=self.conf.overlay_lg_linewidth,
                                facecolor='none')
                except astropy.wcs._wcs.InvalidTransformError:
                    logging.error('could not plot reference sources due to '
                                  'astropy.wcs._wcs.InvalidTransformError; '
                                  'most likely unknown distortion '
                                  'parameters.')

            plt.savefig(framefilename, bbox_inches='tight',
                        pad_inches=0, dpi=self.conf.image_dpi,
                        transparent=True)
            logging.info(('registration map image file for image {:s} '
                          'written to {:s}').format(
                              filename, os.path.abspath(framefilename)))

            plt.close()

        logging.info('create registration overlays with reference stars')

    def add_registration(self, data, extraction_data):
        """
        add registration results to website
        """
        logging.info('adding registration information')

        obsparam = extraction_data[0]['parameters']['obsparam']

        html = self.function_tag+'\n'
        html += ('<H2>Registration</H2>\n'
                 '<P>Registration based on {:s} catalog: \n').format(
            data['catalog'])
        if len(data['badfits']) == 0:
            html += ('<STRONG><FONT COLOR="GREEN">All frames registered '
                     'successfully</FONT></STRONG></P>\n')
        else:
            html += ('<STRONG><FONT COLOR="RED">{:d} files could not be '
                     'registered</FONT></STRONG></P>\n').format(
                         len(data['badfits']))

        if self.conf.show_registration_table:
            html += self.registration_table(data, extraction_data, obsparam)

        if (self.conf.individual_frame_pages and
            self.conf.show_quickview_image and
                self.conf.show_registration_star_map):
            self.registration_maps(data, extraction_data, obsparam)

            for framedata in data['fitresults']:
                # update frame page
                filename = framedata[0]
                if filename in data['goodfits']:
                    resultstring = ('<P><FONT COLOR="GREEN">Registration '
                                    'successful</FONT></P>')
                else:
                    resultstring = ('<P><FONT COLOR="RED">Registration '
                                    'failed</FONT></P>')

                framehtml = (
                    "<!-- Registration -->\n"
                    "<A HREF=\"#registration\" "
                    "ONCLICK=\"toggledisplay('registration');\">"
                    "<H2>Astrometric Registration</H2></A>\n"
                    "<DIV ID=\"registration\" STYLE=\"display: none\">\n"
                    "<TABLE CLASS=\"gridtable\">\n<TR>\n"
                    "<TH>Filename</TH><TH>C<SUB>AS</SUB></TH>"
                    "<TH>C<SUB>XY</SUB></TH>"
                    "<TH>&sigma;<SUB>RA</SUB> (arcsec)</TH>"
                    "<TH>&sigma;<SUB>DEC</SUB> (arcsec)</TH>"
                    "<TH>&chi;<SUP>2</SUP><SUB>Reference</SUB></TH>"
                    "<TH>&chi;<SUP>2</SUP><SUB>Internal</SUB></TH>\n</TR>\n"
                    "<TR><TD>{:s}</TD>"
                    "<TD>{:4.1f}</TD><TD>{:4.1f}</TD>"
                    "<TD>{:5.3f}</TD><TD>{:5.3f}</TD>"
                    "<TD>{:e}</TD><TD>{:e}</TD>\n</TR>\n"
                    "</TABLE>\n"
                    "<STRONG>{:s}</STRONG>"
                    "<DIV CLASS=\"parent_image\">\n"
                    "  <IMG CLASS=\"back_image\" SRC=\"{:s}\" />\n"
                    "  <IMG CLASS=\"front_image\" SRC=\"{:s}\" />\n"
                    "</DIV>\n</DIV>\n\n").format(
                        filename, framedata[1], framedata[2], framedata[3],
                        framedata[4], framedata[5], framedata[6],
                        resultstring,
                        filename+'.'+self.conf.image_file_format,
                        filename+"_astrometry."+self.conf.image_file_format)
                self.append_website(
                    os.path.join(self.conf.diagnostics_path,
                                 '.diagnostics',
                                 '{:s}.html'.format(filename)),
                    framehtml, replace_from='<!-- Registration -->')

        self.append_website(os.path.join(self.conf.diagnostics_path,
                                         self.conf.main_html), html,
                            replace_from=self.function_tag)

        logging.info('registration information added')


class Photometry_Diagnostics(Diagnostics_Html):

    function_tag = "<!-- pp_photometry -->"

    def curve_of_growth_plot(self, data):
        """ create curve of growth plot"""

        logging.info('create curve-of-growth plot')

        parameters = data['parameters']
        growth_filename = os.path.join(self.conf.diagnostics_path,
                                       '.diagnostics', 'curve_of_growth.' +
                                       self.conf.image_file_format)

        f, (ax1, ax2) = plt.subplots(2, sharex=True)

        ax1.set_xlim([min(parameters['aprad']), max(parameters['aprad'])])
        ax1.set_ylabel('Fractional Combined Flux')
        if not parameters['target_only']:
            ax1.plot(parameters['aprad'], data['background_flux'][0],
                     color='black', linewidth=1,
                     label='background sources')
            ax1.fill_between(parameters['aprad'],
                             (data['background_flux'][0] -
                              data['background_flux'][1]),
                             (data['background_flux'][0] +
                              data['background_flux'][1]),
                             color='black', alpha=0.2)
        if not parameters['background_only']:
            ax1.plot(parameters['aprad'], data['target_flux'][0],
                     color='red', linewidth=1,
                     label='target')
            ax1.fill_between(parameters['aprad'],
                             (data['target_flux'][0] -
                              data['target_flux'][1]),
                             (data['target_flux'][0] +
                              data['target_flux'][1]),
                             color='red', alpha=0.2)
        ax1.set_ylim([0, ax1.get_ylim()[1]])
        ax1.plot([data['optimum_aprad'], data['optimum_aprad']],
                 [ax1.get_ylim()[0], ax1.get_ylim()[1]],
                 linewidth=2, color='blue')
        ax1.plot([plt.xlim()[0], plt.xlim()[1]],
                 [data['fluxlimit_aprad'], data['fluxlimit_aprad']],
                 color='black', linestyle='--')
        ax1.grid()
        ax1.legend(loc=4)

        ax2.set_ylim([-0.1, 1.1])
        ax2.set_ylabel('SNR')
        if not parameters['target_only']:
            ax2.errorbar(parameters['aprad'], data['background_snr'],
                         color='black', linewidth=1)
        if not parameters['background_only']:
            ax2.errorbar(parameters['aprad'], data['target_snr'],
                         color='red', linewidth=1)
        ax2.plot([data['optimum_aprad'], data['optimum_aprad']],
                 [plt.ylim()[0], plt.ylim()[1]],
                 linewidth=2, color='blue')
        ax2.grid()
        ax2.set_xlabel('Aperture Radius (px)')
        plt.savefig(growth_filename, format=self.conf.image_file_format,
                    dpi=self.conf.plot_dpi)
        plt.close()
        data['growth_filename'] = growth_filename

        logging.info('curve-of-growth plot created')

    def fwhm_vs_time_plot(self, extraction, data):
        """create fwhm plot"""

        logging.info('create FWHM plot')

        fwhm_filename = os.path.join(self.conf.diagnostics_path,
                                     '.diagnostics',
                                     'fwhm.'+self.conf.image_file_format)

        frame_midtimes = np.array([frame['time'] for frame in extraction])
        fwhm = [np.median(frame['catalog_data']['FWHM_IMAGE'])
                for frame in extraction]
        fwhm_sig = [np.std(frame['catalog_data']['FWHM_IMAGE'])
                    for frame in extraction]

        fig, ax = plt.subplots()

        ax.set_title('Median PSF FWHM per Frame')
        ax.set_xlabel('Minutes after {:s} UT'.format(
            Time(frame_midtimes.min(), format='jd',
                 out_subfmt='date_hm').iso))
        ax.set_ylabel('Point Source FWHM (px)')
        ax.scatter((frame_midtimes-frame_midtimes.min())*1440,
                   fwhm, marker='o',
                   color='black')
        xrange = [plt.xlim()[0], plt.xlim()[1]]
        ax.plot(xrange, [data['optimum_aprad']*2, data['optimum_aprad']*2],
                color='blue')
        ax.set_xlim(xrange)
        ax.set_ylim([0, max([data['optimum_aprad']*2+1, max(fwhm)])])

        ax.grid()
        fig.savefig(fwhm_filename, dpi=self.conf.plot_dpi,
                    format=self.conf.image_file_format)
        data['fwhm_filename'] = fwhm_filename

        # create html map
        if self.conf.individual_frame_pages:
            data['fwhm_map'] = ""
            for i in range(len(extraction)):
                x, y = ax.transData.transform_point(
                    [((frame_midtimes-frame_midtimes.min())*1440)[i],
                     fwhm[i]])
                filename = extraction[i]['fits_filename']
                data['fwhm_map'] += (
                    '<area shape="circle" coords="{:.1f},{:.1f},{:.1f}" '
                    'href="{:s}#{:s}" alt="{:s}" title="{:s}">\n').format(
                        x, fig.bbox.height - y, 5,
                        os.path.join(self.conf.diagnostics_path,
                                     '.diagnostics', filename+'.html'),
                        '',
                        filename, filename)

        logging.info('FWHM plot created')

    def add_photometry(self, data, extraction):
        """
        add photometry results to website
        """
        logging.info('adding photometry information')

        # create curve-of-growth plot
        self.curve_of_growth_plot(data)

        # create fwhm vs time plot
        self.fwhm_vs_time_plot(extraction, data)

        # update index.html
        html = self.function_tag+'\n'
        html += ("<H2>Instrumental Photometry</H2>\n"
                 "<TABLE CLASS=\"gridtable\">\n"
                 "<TR><TH>Photometry Method</TH><TD>{:s}</TD></TR>\n"
                 "<TR><TH>Source Extractor MINAREA (px)</TH>"
                 "<TD>{:.1f}</TD></TR>\n"
                 "<TR><TH>Source Extractor Detection Threshold (&sigma;)"
                 "</TH><TD>{:.1f}</TD></TR>\n").format(
                     {'APER': 'Aperture Photometry'}[_pp_conf.photmode],
                     extraction[0]['parameters']['source_minarea'],
                     extraction[0]['parameters']['sex_snr'])

        if _pp_conf.photmode == 'APER':

            if data['n_target'] > 0 and data['n_bkg'] > 0:
                apsrc = ("{:d} target detections and {:d} "
                         "background detections").format(
                    data['n_target'], data['n_bkg'])
            elif data['n_target'] == 0 and data['n_bkg'] > 0:
                apsrc = "{:d} frames with background detections".format(
                    data['n_bkg'])
            elif data['n_bkg'] == 0 and data['n_target'] > 0:
                apsrc = "{:d} frames with target detections".format(
                    data['n_target'])
            else:
                apsrc = "manually defined"

            html += ("<TR><TH>Aperture Radius (px)</TH>"
                     "<TD>{:.2f}</TD></TR>\n"
                     "<TR><TH>Aperture Radius Basis</TH>"
                     "<TD>{:s}</TD></TR>\n"
                     "<TR><TH>Aperture Radius Strategy</TH>"
                     "<TD>{:s}</TD></TR>\n").format(
                         data['optimum_aprad'],
                         apsrc,
                         data['aprad_strategy']
            )

        html += "</TABLE>\n"
        html += "<P><IMG SRC=\"{:s}\">\n".format(data['growth_filename'])
        if self.conf.individual_frame_pages:
            html += "<IMG SRC=\"{:s}\" USEMAP=\"#FWHM\">\n".format(
                data['fwhm_filename'])
            html += "<MAP NAME=\"#FWHM\">\n{:s}</MAP>\n".format(
                data['fwhm_map'])
        else:
            html += "<IMG SRC=\"{:s}\">\n".format(
                data['fwhm_filename'])

        self.append_website(os.path.join(self.conf.diagnostics_path,
                                         self.conf.main_html), html,
                            replace_from=self.function_tag)
        logging.info('photometry information added')


class Calibration_Diagnostics(Diagnostics_Html):

    function_tag = "<!-- pp_calibrate -->"

    def zeropoint_overview_plot(self, data):
        """produce a plot of magnitude zeropoint as a function of time"""

        logging.info('create zeropoint overview plot')

        times = np.array([dat['obstime'][0] for dat in data['zeropoints']])
        zp = [dat['zp'] for dat in data['zeropoints']]
        zperr = [dat['zp_sig'] for dat in data['zeropoints']]

        fig, ax = plt.subplots()
        ax.errorbar((times-times.min())*1440, zp, yerr=zperr, linestyle='',
                    color='blue', marker='s', capsize=3)
        ax.set_xlabel('Minutes after {:s} UT'.format(
            Time(times.min(), format='jd',
                 out_subfmt='date_hm').iso))
        ax.set_ylabel(
            '{:s}-Band Magnitude Zeropoints (mag)'.format(
                data['filtername']))
        ax.set_ylim([ax.get_ylim()[1], ax.get_ylim()[0]])
        ax.grid()
        fig.savefig(os.path.join(self.conf.diagnostics_path, '.diagnostics',
                                 'zeropoints.'+self.conf.image_file_format),
                    format=self.conf.image_file_format,
                    dpi=self.conf.plot_dpi)
        logging.info('zeropoint overview plot written to {:s}'.format(
            os.path.abspath(os.path.join(self.conf.diagnostics_path,
                                         '.diagnostics', 'zeropoints.' +
                                         self.conf.image_file_format))))
        data['zpplot'] = 'zeropoints.' + self.conf.image_file_format

        # create html map
        if self.conf.individual_frame_pages:
            data['zpplotmap'] = ""
            for i in range(len(times)):
                x, y = ax.transData.transform_point(
                    [((times-times.min())*1440)[i],
                     [dat['zp'] for dat in data['zeropoints']][i]])
                filename = data['zeropoints'][i]['filename'][:-4]+'fits'
                data['zpplotmap'] += (
                    '<area shape="circle" coords="{:.1f},{:.1f},{:.1f}" '
                    'href="{:s}#{:s}" alt="{:s}" title="{:s}">\n').format(
                        x, fig.bbox.height - y, 5,
                        os.path.join(self.conf.diagnostics_path,
                                     '.diagnostics', filename+'.html'),
                        'calibration_overview',
                        filename, filename)
        logging.info('zeropoint overview plot created')

    def phot_calibration_plot(self, data, idx):
        """produce diagnostic plot for each frame"""

        logging.info('create photometric calibration overview plot')

        f, (ax1, ax3) = plt.subplots(2)
        plt.subplots_adjust(hspace=0.3)

        ax1.set_title('{:s}: {:s}-band from {:s}'.format(
            data['catalogs'][idx].catalogname,
            data['filtername'],
            data['ref_cat'].catalogname))
        ax1.set_xlabel('Number of Reference Stars')
        ax1.set_ylabel('Magnitude Zeropoint', fontdict={'color': 'red'})

        zp_idx = data['zeropoints'][idx]['zp_idx']
        clipping_steps = data['zeropoints'][idx]['clipping_steps']

        x = [len(clipping_steps[i][3]) for i
             in range(len(clipping_steps))]

        ax1.errorbar(x, [clipping_steps[i][0] for i
                         in range(len(clipping_steps))],
                     yerr=[clipping_steps[i][1] for i
                           in range(len(clipping_steps))], color='red')
        ax1.set_ylim(ax1.get_ylim()[::-1])  # reverse y axis
        ax1.plot([len(clipping_steps[zp_idx][3]),
                  len(clipping_steps[zp_idx][3])],
                 ax1.get_ylim(), color='black')

        ax1.grid(linestyle='--')

        ax2 = ax1.twinx()
        ax2.plot(x, [clipping_steps[i][2] for i
                     in range(len(clipping_steps))],
                 color='blue')
        ax2.set_ylabel(r'reduced $\chi^2$', fontdict={'color': 'blue'})
        ax2.set_yscale('log')
        ax2.plot([min(x), max(x)], [1, 1], linestyle='dotted', color='blue')

        # residual plot
        ax3.set_xlabel('Reference Star Magnitude')
        ax3.set_ylabel('Calibration-Reference (mag)')

        match = data['zeropoints'][idx]['match']
        x = match[0][0][clipping_steps[zp_idx][3]]
        residuals = (match[1][0][clipping_steps[zp_idx][3]]
                     + clipping_steps[zp_idx][0]
                     - match[0][0][clipping_steps[zp_idx][3]])
        residuals_sig = np.sqrt(match[1][1][clipping_steps[zp_idx][3]]**2
                                + clipping_steps[zp_idx][1]**2)

        ax3.errorbar(x, residuals, yerr=residuals_sig, color='black',
                     marker='o', linestyle='')
        x_range = ax3.get_xlim()
        ax3.plot(x_range, [0, 0], color='black', linestyle='--')
        ax3.set_xlim(x_range)
        ax3.set_ylim(ax3.get_ylim()[::-1])  # reverse y axis

        ax3.grid(linestyle='--')

        plotfilename = os.path.join(self.conf.diagnostics_path,
                                    '.diagnostics',
                                    '{:s}_photcal.{:s}'.format(
                                        data['catalogs'][idx].catalogname,
                                        self.conf.image_file_format))
        plt.savefig(plotfilename, dpi=self.conf.plot_dpi,
                    format=self.conf.image_file_format)
        data['zeropoints'][idx]['plotfilename'] = plotfilename
        logging.info('create photometric calibration overview plot')

    def calibration_raw_data_tables(self, dat):
        """build table with photometric calibration raw data"""

        logging.info('create calibration data table')

        html = "<TD><TABLE CLASS=\"gridtable\">\n<TR>\n"
        html += ("<TH>Idx</TH><TH>Source Name</TH><TH>RA</TH><TH>Dec</TH>"
                 "<TH>Catalog (mag)</TH>"
                 "<TH>Instrumental (mag)</TH><TH>Calibrated (mag)</TH>"
                 "<TH>Residual (mag</TH>\n</TR>\n")
        for i, idx in enumerate(dat['zp_usedstars']):
            name = str(dat['match'][0][2][idx])
            if isinstance(name, bytes):
                name = name.decode('utf8')
            html += ("<TR><TD>{:d}</TD><TD>{:s}</TD><TD>{:12.8f}</TD>"
                     "<TD>{:12.8f}</TD><TD>{:.3f}+-{:.3f}</TD>"
                     "<TD>{:.3f}+-{:.3f}</TD>"
                     "<TD>{:.3f}+-{:.3f}</TD><TD>{:.3f}</TD>"
                     "</TR>").format(
                         i+1, name,
                         dat['match'][0][3][idx],
                         dat['match'][0][4][idx],
                         dat['match'][0][0][idx],
                         dat['match'][0][1][idx],
                         dat['match'][1][0][idx],
                         dat['match'][1][1][idx],
                         dat['zp']+dat['match'][1][0][idx],
                         np.sqrt(dat['zp_sig']**2 +
                                 dat['match'][1][1][idx]**2),
                         (dat['zp']+dat['match'][1][0][idx]) -
                         dat['match'][0][0][idx])
        html += "</TABLE>\n"

        logging.info('calibration data table created')
        return html

    def calibration_star_maps(self, dat):
        """create thumbnail images with calibration stars marked"""

        logging.info('create calibration reference star overlay')

        fits_filename = (dat['filename'][:dat['filename'].find('.ldac')]
                         + '.fits')

        imgdat = fits.open(fits_filename,
                           ignore_missing_end=True)[0].data
        resize_factor = min(
            1.,
            self.conf.image_size_lg_px/np.max(imgdat.shape))

        header = fits.open(fits_filename,
                           ignore_missing_end=True)[0].header

        # turn relevant header keys into floats
        # astropy.io.fits bug
        for key, val in list(header.items()):
            if 'CD1_' in key or 'CD2_' in key or \
               'CRVAL' in key or 'CRPIX' in key or \
               'EQUINOX' in key:
                header[key] = float(val)

        plt.figure(figsize=(self.conf.image_size_lg_in,
                            self.conf.image_size_lg_in))
        # create fake image to ensure image dimensions and margins
        img = plt.imshow(np.ones((self.conf.image_size_lg_px,
                                  self.conf.image_size_lg_px))*np.nan,
                         origin='lower')

        # remove axes
        plt.axis('off')
        img.axes.get_xaxis().set_visible(False)
        img.axes.get_yaxis().set_visible(False)

        # plot reference sources
        if len(dat['match'][0][3]) > 0 and len(dat['match'][0][4]) > 0:
            try:
                w = wcs.WCS(header)
                world_coo = [[dat['match'][0][3][idx],
                              dat['match'][0][4][idx]]
                             for idx in dat['zp_usedstars']]
                img_coo = w.wcs_world2pix(world_coo, True)
                plt.scatter([c[0]*resize_factor for c in img_coo],
                            [c[1]*resize_factor for c in img_coo],
                            s=10, marker='o', edgecolors='red',
                            linewidth=0.3,
                            facecolor='none')
                for i in range(len(dat['zp_usedstars'])):
                    plt.annotate(str(i+1),
                                 xy=((img_coo[i][0]*resize_factor)+15,
                                     img_coo[i][1]*resize_factor),
                                 color='red',
                                 horizontalalignment='left',
                                 verticalalignment='center')
            except astropy.wcs._wcs.InvalidTransformError:
                logging.error('could not plot reference sources due to '
                              'astropy.wcs._wcs.InvalidTransformError; '
                              'most likely unknown distortion '
                              'parameters.')

        catframe = os.path.join(
            self.conf.diagnostics_path, '.diagnostics',
            '{:s}.fits_reference_stars.{:s}'.format(
                dat['filename'][:dat['filename'].find('.ldac')],
                self.conf.image_file_format))
        plt.savefig(catframe, format=self.conf.image_file_format,
                    bbox_inches='tight',
                    pad_inches=0, dpi=self.conf.image_dpi, transparent=True)
        plt.close()
        logging.info('calibration reference star overlay created')

    def add_calibration(self, data, instrumental=False):
        """
        wrapper to add calibration results to diagnostics website
        """
        logging.info('adding photometric calibration information')

        html = self.function_tag+'\n'
        html += "<H2>Photometric Calibration</H2>\n"

        if not instrumental:
            # create zeropoint overview plot
            self.zeropoint_overview_plot(data)

            # main diagnostics website content
            html += ("<TABLE CLASS=\"gridtable\">\n"
                     "<TR><TH>Reference Catalog</TH><TD>{:s}</TD></TR>\n"
                     "<TR><TH>Reference Catalog History</TH>"
                     "<TD>{:s}</TD></TR>\n"
                     "<TR><TH>Target Filter</TH><TD>{:s}</TD></TR>\n"
                     "</TABLE>\n").format(
                         data['ref_cat'].catalogname,
                         data['ref_cat'].history,
                         data['filtername'])

            # build overview table
            html += ("<P><TABLE CLASS=\"gridtable\">\n<TR>\n"
                     "<TH>Filename</TH><TH>Zeropoint (mag)</TH>"
                     "<TH>&sigma; (mag)</TH>"
                     "<TH>N<SUP>*</SUP><SUB>used</SUB></TH>"
                     "<TH>N<SUP>*</SUP><SUB>matched</SUB></TH>\n</TR>\n")

            for idx, dat in enumerate(data['zeropoints']):
                # update frame pages
                if self.conf.individual_frame_pages:
                    framename = "<A HREF=\"{:s}\">{:s}</A>".format(
                        os.path.join(self.conf.diagnostics_path,
                                     '.diagnostics',
                                     dat['filename'][:-4]+'fits.html'),
                        dat['filename'][:-4]+'fits')

                    framehtml = ("<!-- Calibration -->\n"
                                 "<A HREF=\"#calibration_overview\" "
                                 "ONCLICK=\"toggledisplay"
                                 "('calibration_overview');\">"
                                 "<H2>Photometric Calibration</H2></A>\n"
                                 "<DIV ID=\"calibration_overview\" "
                                 "STYLE=\"display: none\"\>\n")

                    if dat['success']:
                        framehtml += (
                            "<P><TABLE CLASS=\"gridtable\">\n"
                            "<TR><TH>Reference Catalog</TH>"
                            "<TD>{:s}</TD></TR>\n"
                            "<TR><TH>Reference Catalog History</TH>"
                            "<TD>{:s}</TD></TR>\n"
                            "<TR><TH>Target Filter</TH>"
                            "<TD>{:s}</TD></TR>\n"
                            "<TR><TH>Zeropoint (mag)</TH>"
                            "<TD>{:7.4f}+-{:.4f}</TD></TR>\n"
                            "<TR><TH>N<SUP>*</SUP><SUB>used</SUB>"
                            "</TH>"
                            "<TD>{:d}</TD></TR>\n"
                            "<TH>N<SUP>*</SUP><SUB>matched</SUB></TH>"
                            "<TD>{:d}</TD></TR>\n").format(
                                data['ref_cat'].catalogname,
                                data['ref_cat'].history,
                                data['filtername'],
                                dat['zp'], dat['zp_sig'],
                                dat['zp_nstars'],
                                len(dat['match'][0][0]))
                    else:
                        framehtml += (
                            "<P><TABLE CLASS=\"gridtable\">\n"
                            "<TR><TH>Reference Catalog</TH>"
                            "<TD>{:s}</TD></TR>\n"
                            "<TR><TH>Reference Catalog History</TH>"
                            "<TD>{:s}</TD></TR>\n"
                            "<TR><TH>Target Filter</TH>"
                            "<TD>{:s}</TD></TR>\n"
                            "<TR><TH>Zeropoint (mag)</TH>"
                            "<TD>{:7.4f}+-{:.4f}</TD></TR>\n"
                            "<TR><TH>N<SUP>*</SUP><SUB>used</SUB>"
                            "</TH>"
                            "<TD>{:d}</TD></TR>\n"
                            "<TH>N<SUP>*</SUP><SUB>matched</SUB></TH>"
                            "<TD>{:d}</TD></TR>\n").format(
                                data['ref_cat'].catalogname,
                                data['ref_cat'].history,
                                data['filtername'],
                                np.nan, np.nan, 0, len(dat['match'][0][0]))

                    framehtml += "</TABLE></P>\n"

                    # frame calibration data
                    catframe = os.path.join(
                        self.conf.diagnostics_path, '.diagnostics',
                        '{:s}.fits_reference_stars.{:s}'.format(
                            dat['filename'][:dat['filename'].find('.ldac')],
                            self.conf.image_file_format))

                    # build individual calibration plots
                    if (self.conf.show_phot_calibration_plots and
                            data['zeropoints'][idx]['success']):
                        self.phot_calibration_plot(data, idx)
                        framehtml += (
                            "<A HREF=\"#calibration_plot\" "
                            "ONCLICK=\"toggledisplay"
                            "('calibration_plot');\">"
                            "<H3>Calibration Analysis</H3></A>\n"
                            "<DIV ID=\"calibration_plot\">\n"
                            "<P><IMG SRC={:s} \></DIV>\n").format(
                            dat['plotfilename'].split(os.path.sep)[-1])

                    # build individual catalog maps
                    if (self.conf.individual_frame_pages and
                        self.conf.show_quickview_image and
                            self.conf.show_calibration_star_map and
                            data['zeropoints'][idx]['success']):
                        self.calibration_star_maps(dat)
                        framehtml += (
                            "<A HREF=\"#calibration_starmap\" "
                            "ONCLICK=\"toggledisplay"
                            "('calibration_starmap');\">"
                            "<H3>Calibration Map</H3></A>\n"
                            "<DIV ID=\"calibration_starmap\" "
                            "STYLE=\"display: none\"\>\n"
                            "<DIV CLASS=\"parent_image\">\n"
                            "<IMG CLASS=\"back_image\" SRC=\"{:s}\" />\n"
                            "<IMG CLASS=\"front_image\" SRC=\"{:s}\" />\n"
                            "</DIV></DIV>\n").format(
                                dat['filename'][:-4]+'fits.' +
                            self.conf.image_file_format,
                                catframe.split(os.path.sep)[-1])

                    # build individual catalog data table websites
                    if (self.conf.show_calibration_star_table and
                            data['zeropoints'][idx]['success']):
                        framehtml += (
                            "<A HREF=\"#calibration_table\" "
                            "ONCLICK=\"toggledisplay"
                            "('calibration_table');\">"
                            "<H3>Calibration Data Table</H3></A>\n"
                            "<DIV ID=\"calibration_table\" "
                            "STYLE=\"display: none\"\>\n"
                            "<P>{:s}</DIV>\n").format(
                            self.calibration_raw_data_tables(dat))

                    framehtml += "</DIV>\n\n"

                    self.append_website(
                        os.path.join(self.conf.diagnostics_path,
                                     '.diagnostics', '{:s}.html'.format(
                                         dat['filename'][:-4]+'fits')),
                        framehtml, replace_from='<!-- Calibration -->')

                else:
                    framename = dat['filename'][:-4]+'fits'

                if data['zeropoints'][idx]['success']:
                    html += ("<TR><TD>{:s}</TD>"
                             "<TD>{:7.4f}</TD><TD>{:7.4f}</TD><TD>{:d}</TD>"
                             + "<TD>{:d}</TD>\n</TR>").format(
                                 framename, dat['zp'],
                                 dat['zp_sig'], dat['zp_nstars'],
                                 len(dat['match'][0][0]))
                else:
                    html += ("<TR><TD>{:s}</TD>"
                             "<TD>{:7.4f}</TD><TD>{:7.4f}</TD><TD>{:d}</TD>"
                             + "<TD>{:d}</TD>\n</TR>").format(
                                 framename, np.nan, np.nan, 0,
                                 len(dat['match'][0][0]))

            html += "</TABLE></P>\n"

            if self.conf.individual_frame_pages:
                html += ("<P><IMG SRC=\"{:s}\" "
                         "USEMAP=\"#Zeropoints\">\n").format(
                             os.path.join(self.conf.diagnostics_path,
                                          '.diagnostics', data['zpplot']))
                html += "<MAP NAME=\"#Zeropoints\">\n{:s}</MAP>\n".format(
                    data['zpplotmap'])
            else:
                html += "<P><IMG SRC=\"{:s}\">\n".format(
                    os.path.join(self.conf.diagnostics_path,
                                 '.diagnostics', data['zpplot']))
        else:
            html += ("Instrumental magnitudes are reported "
                     "(filter used: {:s})\n").format(
                         str(data['filtername']))

        self.append_website(os.path.join(self.conf.diagnostics_path,
                                         self.conf.main_html), html,
                            replace_from=self.function_tag)
        logging.info('photometric calibration information added')


class Distill_Diagnostics(Diagnostics_Html):

    function_tag = "<!-- pp_distill -->"

    def lightcurve_plots(self, data):
        """build lightcurve plots for targets"""
        logging.info('create lightcurve plots for all targets')

        data['lightcurveplots'] = {}
        data['lightcurveplots']['maps'] = {}
        for target in data['targetnames']:

            if len(data[target]) == 0:
                continue

            logging.info('create lightcurve plot for {:s}'.format(target))

            midtimes = np.array([dat[9][0] for dat in data[target]])

            fig, ax = plt.subplots()
            ax.set_title(target.replace('_', ' '))
            ax.set_xlabel('Minutes after {:s} UT'.format(
                Time(midtimes.min(), format='jd',
                     out_subfmt='date_hm').iso))
            ax.set_ylabel('Magnitude')
            ax.errorbar((midtimes-midtimes.min())*1440,
                        [dat[7] for dat in data[target]],
                        yerr=[dat[8] for dat in data[target]],
                        linestyle='', color='red',
                        marker='o', capsize=3)
            ax.set_ylim([ax.get_ylim()[1], ax.get_ylim()[0]])
            ax.set_xticklabels = [Time(t, format='jd').iso
                                  for t in plt.xticks()[0]]
            ax.grid()

            fig.savefig(os.path.join(self.conf.diagnostics_path,
                                     '.diagnostics', '{:s}.{:s}'.format(
                                         target.translate(
                                             _pp_conf.target2filename),
                                         self.conf.image_file_format)),
                        format=self.conf.image_file_format,
                        dpi=self.conf.plot_dpi)
            logging.info('lightcurve plot for {:s} written to {:s}'.format(
                target, os.path.abspath(
                    os.path.join(self.conf.diagnostics_path,
                                 '.diagnostics', '{:s}.{:s}'.format(
                                     target.translate(
                                         _pp_conf.target2filename),
                                     self.conf.image_file_format)))))
            data['lightcurveplots'][target] = os.path.join(
                self.conf.diagnostics_path, '.diagnostics',
                '{:s}.{:s}'.format(
                    target.translate(_pp_conf.target2filename),
                    self.conf.image_file_format))

            # create html map
            if self.conf.individual_frame_pages:
                data['lightcurveplots']['maps'][target] = ""
                for i in range(len(midtimes)):
                    x, y = ax.transData.transform_point(
                        [((midtimes-midtimes.min())*1440)[i],
                         [dat[7] for dat in data[target]][i]])
                    data['lightcurveplots']['maps'][target] += (
                        '<area shape="circle" '
                        'coords="{:.1f},{:.1f},{:.1f}" '
                        'href="{:s}#{:s}" alt="{:s}" '
                        'title="{:s}">\n').format(
                            x, fig.bbox.height - y, 5,
                            os.path.join(
                                self.conf.diagnostics_path,
                                '.diagnostics', data[
                                    target][i][10][:-4]+'fits.html'),
                            target,
                            data[target][i][10][:-4] + 'fits',
                            data[target][i][10][:-4] + 'fits')

        logging.info('lightcurve plots for all targets created')

    def thumbnail_images(self, data):
        """build thumbnail images for each frame/target"""
        logging.info('create thumbnail images and overlays for all targets')

        data['thumbnailplots'] = {}
        data['thumbnailoverlays'] = {}

        for target in data['targetnames']:

            data['thumbnailplots'][target] = []
            data['thumbnailoverlays'][target] = []

            if sys.version_info < (3, 0):
                target = str(target)

            data['thumbnailplots'][target] = []
            for dat in data[target]:
                for fitsfilename in ['.fits', '.fit']:
                    fitsfilename = (dat[10][:dat[10].find('.ldac')] +
                                    fitsfilename)
                    if os.path.isfile(fitsfilename):
                        break
                hdulist = fits.open(fitsfilename, ignore_missing_end=True)

                logging.info('create thumbnail image for {:s}/{:s}'.format(
                    target, fitsfilename))

                # turn relevant header keywords into floats
                # should be fixed in astropy.wcs
                for key, val in list(hdulist[0].header.items()):
                    if 'CD1' in key or 'CD2' in key or \
                       'CRVAL' in key or 'CRPIX' in key or \
                       'EQUINOX' in key:
                        hdulist[0].header[key] = float(val)

                w = wcs.WCS(hdulist[0].header)
                obj_x, obj_y = dat[11], dat[12]
                image_coords = w.wcs_world2pix(np.array([[dat[1], dat[2]]]),
                                               True)
                exp_x, exp_y = image_coords[0][0], image_coords[0][1]

                # check if thumbnail area is close to image edge
                if (exp_x < self.conf.image_size_thumb_px/2 or
                    exp_x > (hdulist[0].data.shape[0] -
                             self.conf.image_size_thumb_px/2) or
                    exp_y < self.conf.image_size_thumb_px/2 or
                        exp_y > (hdulist[0].data.shape[1] -
                                 self.conf.image_size_thumb_px/2)):
                    # create margin around image allowing for any cropping
                    composite = np.ones(
                        (hdulist[0].data.shape[0] +
                         2*self.conf.image_size_thumb_px,
                         hdulist[0].data.shape[1] +
                         2*self.conf.image_size_thumb_px))*np.nan
                    # insert image
                    composite[
                        self.conf.image_size_thumb_px:
                        self.conf.image_size_thumb_px +
                        hdulist[0].data.shape[0],
                        self.conf.image_size_thumb_px:
                        self.conf.image_size_thumb_px +
                        hdulist[0].data.shape[1]] = (
                            hdulist[0].data)

                    # extract thumbnail data accordingly
                    thumbdata = composite[
                        int(self.conf.image_size_thumb_px+obj_y -
                            self.conf.image_size_thumb_px/2):
                        int(self.conf.image_size_thumb_px+obj_y +
                            self.conf.image_size_thumb_px/2),
                        int(self.conf.image_size_thumb_px+obj_x -
                            self.conf.image_size_thumb_px/2):
                        int(self.conf.image_size_thumb_px+obj_x +
                            self.conf.image_size_thumb_px/2)]
                    del composite
                else:
                    thumbdata = hdulist[0].data[
                        int(obj_y-self.conf.image_size_thumb_px/2):
                        int(obj_y+self.conf.image_size_thumb_px/2),
                        int(obj_x-self.conf.image_size_thumb_px/2):
                        int(obj_x+self.conf.image_size_thumb_px/2)]

                norm = ImageNormalize(
                    thumbdata, interval=ZScaleInterval(),
                    stretch={'linear': LinearStretch(),
                             'log': LogStretch()}[
                                 self.conf.image_stretch])

                # create plot
                fig = plt.figure(figsize=(self.conf.image_size_thumb_in,
                                          self.conf.image_size_thumb_in))
                img = plt.imshow(thumbdata, cmap='gray',
                                 origin='lower', norm=norm)
                # remove axes
                plt.axis('off')
                img.axes.get_xaxis().set_visible(False)
                img.axes.get_yaxis().set_visible(False)

                thumbfilename = (
                    os.path.join(self.conf.diagnostics_path,
                                 '.diagnostics', target.translate(
                                     _pp_conf.target2filename) + '_' +
                                 fitsfilename[:fitsfilename.
                                              find('.fit')] +
                                 '_thumb.'+self.conf.image_file_format))

                plt.savefig(thumbfilename,
                            format=self.conf.image_file_format,
                            bbox_inches='tight',
                            pad_inches=0, dpi=self.conf.image_dpi)
                plt.close()

                # create overlay

                fig = plt.figure(figsize=(self.conf.image_size_thumb_in,
                                          self.conf.image_size_thumb_in))
                img = plt.imshow(
                    np.ones((self.conf.image_size_thumb_px,
                             self.conf.image_size_thumb_px))*np.nan,
                    origin='lower')

                # remove axes
                plt.axis('off')
                img.axes.get_xaxis().set_visible(False)
                img.axes.get_yaxis().set_visible(False)

                # add image filename
                plt.annotate('{:s}'.format(fitsfilename), (3, 5),
                             color='white',
                             fontsize=self.conf.thumb_fontsize)
                # add target name
                plt.annotate('{:s}'.format(target.replace('_', ' ')),
                             (3, self.conf.image_size_thumb_px-10),
                             color='white',
                             fontsize=self.conf.thumb_fontsize)
                # add scales
                pixelscales = (np.fabs(w.pixel_scale_matrix[0][0])*3600,
                               np.fabs(w.pixel_scale_matrix[1][1])*3600)
                plt.plot([self.conf.image_size_thumb_px-3,
                          self.conf.image_size_thumb_px-3],
                         [3, 3+self.conf.thumb_scalelength/pixelscales[1]],
                         color='white')
                plt.plot([self.conf.image_size_thumb_px-3,
                          self.conf.image_size_thumb_px-3 -
                          self.conf.thumb_scalelength/pixelscales[0]],
                         [3, 3],
                         color='white')
                plt.annotate('{:d}\"'.format(self.conf.thumb_scalelength),
                             xy=(self.conf.image_size_thumb_px-3 -
                                 0.5*self.conf.thumb_scalelength /
                                 pixelscales[0],
                                 3+0.5*self.conf.thumb_scalelength /
                                 pixelscales[1]),
                             color='white',
                             fontsize=self.conf.thumb_fontsize,
                             horizontalalignment='center',
                             verticalalignment='center')
                # add compass?

                # place aperture
                if _pp_conf.photmode == 'APER':
                    aprad = float(hdulist[0].header['APRAD'])
                    targetpos = plt.Circle(
                        (self.conf.image_size_thumb_px/2,
                         self.conf.image_size_thumb_px/2),
                        aprad, ec='red', fc='none',
                        linewidth=self.conf.thumb_linewidth)
                else:
                    targetpos = plt.Rectangle(
                        (self.conf.image_size_thumb_px/2-7,
                         self.conf.image_size_thumb_px/2-7),
                        15, 15, ec='red', fc='none',
                        linewidth=self.conf.thumb_linewidth)
                plt.gca().add_patch(targetpos)

                # place predicted position (if within thumbnail)
                if ((abs(exp_x-obj_x) <=
                     self.conf.image_size_thumb_px/2) and
                        (abs(exp_y-obj_y) <=
                         self.conf.image_size_thumb_px/2)):
                    plt.scatter(exp_x-obj_x+self.conf.image_size_thumb_px/2,
                                exp_y-obj_y+self.conf.image_size_thumb_px/2,
                                marker=self.conf.thumb_predicted_pos_marker,
                                s=self.conf.thumb_predicted_pos_size,
                                color=self.conf.thumb_predicted_pos_color)

                thumbfileoverlay = os.path.join(
                    self.conf.diagnostics_path, '.diagnostics',
                    target.translate(_pp_conf.target2filename) + '_' +
                    fitsfilename[:fitsfilename.find('.fit')] +
                    '_thumb_overlay.' + self.conf.image_file_format)

                plt.savefig(thumbfileoverlay,
                            format=self.conf.image_file_format,
                            bbox_inches='tight', transparent=True,
                            pad_inches=0, dpi=self.conf.image_dpi)
                plt.close()
                hdulist.close()
                data['thumbnailplots'][target].append((fitsfilename,
                                                       thumbfilename))
                data['thumbnailoverlays'][target].append((fitsfilename,
                                                          thumbfileoverlay))
        logging.info(('thumbnail images and overlays for all targets '
                      'created'))

    def target_animations(self, data):
        """assemble thumbnails to gif animations"""
        logging.info('create target thumbnail animations')

        data['gifs'] = {}

        for target in data['targetnames']:
            gif_filename = '{:s}.gif'.format(
                target.translate(_pp_conf.target2filename))
            logging.info('converting images to gif: {:s}'.format(
                gif_filename))
            root = os.getcwd()
            os.chdir(os.path.join(self.conf.diagnostics_path,
                                  '.diagnostics'))
            try:
                convert = subprocess.Popen(
                    ['convert', '-delay', '50',
                     ('{:s}*thumb.{:s}'.format(target.translate(
                         _pp_conf.target2filename),
                         self.conf.image_file_format)),
                     '-loop', '0',
                     ('{:s}'.format(gif_filename))])

                convert.wait()
            except:
                logging.warning('could not produce gif animation for '
                                + 'target {:s}'.format(target))
            data['gifs'][target] = os.path.join(
                self.conf.diagnostics_path, '.diagnostics', gif_filename)
            os.chdir(root)

        logging.info('target thumbnail animations created')

    def add_results(self, data, imagestretch='linear'):
        """
        add results to website
        """
        logging.info('adding distill results')

        self.lightcurve_plots(data)

        self.thumbnail_images(data)

        if self.conf.show_target_animations:
            self.target_animations(data)

        # self.add_frame_report(data)

        html = self.function_tag+'\n'
        html += "<H2>Photometry Results</H2>\n"

        for idx, target in enumerate(data['targetnames']):
            html += ("<A HREF=\"#{:s}\" "
                     "ONCLICK=\"toggledisplay"
                     "('{:s}');\">"
                     "<H3>{:s}</H3></A>\n"
                     "<DIV ID=\"{:s}\" "
                     "STYLE=\"display: none\"\>\n"
                     "<TABLE BORDER=\"0\"><TR><TD>\n").format(
                         target, target,
                         target.replace('_', ' '), target)
            if self.conf.individual_frame_pages:
                try:
                    html += ("<IMG SRC=\"{:s}\" USEMAP=\"#{:s}\">\n"
                             "<MAP NAME=\"{:s}\">\n{:s}</MAP>\n\n").format(
                                 data['lightcurveplots'][target], target,
                                 target,
                                 data['lightcurveplots']['maps'][target])
                except KeyError:
                    pass
            else:
                try:
                    html += ("<IMG SRC=\"{:s}\">\n"
                             "<MAP NAME=\"{:s}\">\n{:s}</MAP>\n\n").format(
                                 data['lightcurveplots'][target],
                                 target, target)
                except KeyError:
                    pass
            html += "</TD>\n<TD>"
            if self.conf.show_target_animations:
                html += "<IMG SRC=\"{:s}\" \>".format(
                    data['gifs'][target])
            html += "</TD></TR></TABLE>\n</DIV>\n"

            # update framepages
            if self.conf.individual_frame_pages:
                for fidx, framedat in enumerate(data[target]):
                    fitsfilename = framedat[10][:-4]+'fits'

                    # identify next/previous frame in cyclical
                    # data['targetframes'] list
                    previousframe = data['targetframes'][target][
                        (np.where(np.array(data['targetframes'][target]) ==
                                  fitsfilename)[0][0]-1) %
                        len(data['targetframes'][target])]
                    nextframe = data['targetframes'][target][
                        (np.where(np.array(data['targetframes'][target]) ==
                                  fitsfilename)[0][0]+1) %
                        len(data['targetframes'][target])]

                    assert (data['thumbnailplots']
                            [target][fidx][0].strip() ==
                            fitsfilename.strip())

                    if self.conf.show_target_animations:
                        animation_button = (
                            "<P ALIGN=\"center\">"
                            "<BUTTON ONCLICK=\""
                            "window.open('{:s}', 'targetWindow', "
                            "'height=376,width=376,status=no,location=no,"
                            "scrollbars=no,toolbar=no,menubar=no')\">"
                            "open animation in new window</BUTTON>"
                            "</P>\n").format(
                                data['gifs'][target].split(os.path.sep)[-1])
                    else:
                        animation_button = ""

                    framehtml = (
                        "<!-- Results {:s} -->\n"
                        "<A HREF=\"#{:s}\" "
                        "ONCLICK=\"toggledisplay"
                        "('{:s}');\">"
                        "<H1>{:s}  Photometry</H1></A>\n"
                        "<DIV ID=\"{:s}\"\>\n"
                        "<TABLE BORDER=\"0\"><TR><TD>"
                        "<TABLE CLASS=\"gridtable\">\n"
                        "<TR><TH>Apparent Magnitude</TH>"
                        "<TD>{:.2f} +- {:.2f}</TD></TR>\n"
                        "<TR><TH>Target RA (deg)</TH>"
                        "<TD>{:7.5f}</TD></TR>\n"
                        "<TR><TH>Target Dec (deg)</TH>"
                        "<TD>{:6.4f}</TD></TR>\n"
                        "<TR><TH>RA Offset from Prediction (\")</TH>"
                        "<TD>{:.2f}</TD></TR>\n"
                        "<TR><TH>Dec Offset from Prediction (\")</TH>"
                        "<TD>{:.2f}</TD></TR>\n"
                        "<TH>Target Source Flag</TH>"
                        "<TD>{:d}</TD></TR>\n"
                        "</TABLE>\n"
                        "{:s}"
                        "<P ALIGN=\"center\">"
                        "<BUTTON ONCLICK=\""
                        "window.open('{:s}', 'targetWindow', "
                        "'height=480,width=640,status=no,location=no,"
                        "scrollbars=no,toolbar=no,menubar=no')\">"
                        "open lightcurve in new window</BUTTON></P>"
                        "<P ALIGN=\"center\">"
                        "<BUTTON ONCLICK=\"toggledisplay"
                        "('{:s}');\">toggle overlay</BUTTON></P>"
                        "</TD><TD>\n"
                        "<DIV CLASS=\"parent_image\">\n"
                        "<IMG CLASS=\"back_image\" SRC=\"{:s}\" />\n"
                        "<IMG ID=\"overlay_{:s}\" "
                        "CLASS=\"front_image\" SRC=\"{:s}\" />\n"
                        "</DIV>\n"
                        "<DIV ALIGN=\"center\">"
                        "<A HREF=\"{:s}#{:s}\">"
                        "&laquo; previous frame &laquo;</A> | "
                        "<A HREF=\"{:s}#{:s}\">"
                        "&raquo; next frame &raquo;</A></DIV>\n"
                        "</TD></TR></TABLE></DIV>\n").format(
                            target, target, target,
                            target.replace('_', ' '), target,
                            framedat[7], framedat[8],
                            framedat[3], framedat[4],
                            (framedat[1]-framedat[3])*3600,
                            (framedat[2]-framedat[4])*3600,
                            int(framedat[14]),
                            animation_button,
                            data['lightcurveplots'][
                                target].split('/')[-1],
                            'overlay_'+target,
                            data['thumbnailplots'][
                                target][fidx][1].split('/')[-1],
                            target,
                            data['thumbnailoverlays'][
                                target][fidx][1].split('/')[-1],
                            previousframe+'.html', target,
                            nextframe+'.html', target)

                    self.append_website(os.path.join(
                        self.conf.diagnostics_path, '.diagnostics',
                        '{:s}.html'.format(
                            fitsfilename)),
                        framehtml,
                        replace_from='<!-- Results {:s} -->'.format(target))

        self.append_website(os.path.join(self.conf.diagnostics_path,
                                         self.conf.main_html), html,
                            replace_from=self.function_tag)
        logging.info('distill results added')


preparation = Prepare_Diagnostics()
registration = Registration_Diagnostics()
photometry = Photometry_Diagnostics()
calibration = Calibration_Diagnostics()
distill = Distill_Diagnostics()

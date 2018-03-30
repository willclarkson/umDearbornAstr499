#
# testGaia.py
#

# 2018-03-30 WIC - utility to extract gaia data about selected
# location and convert to astropy table. Hacked out of code originally
# for project GO-12020.

# for low-tech timing of the query
import time, os

# useful imports
from astropy.table import Table, join
import numpy as np

# for querying an external catalog
from astroquery.vizier import Vizier
from astropy.coordinates import Angle

# for parsing box region for co-ordinate import
import astropy.units as u
import astropy.coordinates as coord

class GaiaCat(object):

    """Object to hold the gaia query"""

    def __init__(self, ra=250., dec=-18., \
                     wid=0.01, hgt=0.01, \
                     cooFrame='icrs', \
                     boxUnits='deg', \
                     maxRows=10000, \
                     Verbose=True, runOnInit=True, \
                     srcCat='gaia', \
                     forceQuery=True):

        # our field center and box dimensions
        self.ra  = ra
        self.dec = dec
        
        # ra, dec in degrees. Will be overwritten when the coords are
        # parsed.
        self.raDeg = 250.
        self.decDeg = -18.

        # coordinate units for field center
        self.cenUnit='deg'

        self.boxWidth = np.copy(wid)
        self.boxHeight = np.copy(hgt)
        self.boxUnits = boxUnits[:]

        # coord frame
        self.cooFrame = cooFrame[:]

        # astropy-coordinates object for the field center
        self.fieldCen = None

        # a few sensible defaults - don't let the user query the
        # entire sky!
        self.maxBoxWidth = 0.1
        self.maxBoxHeight = 0.1
        self.maxUnitWidth = 'deg'
        self.maxUnitHeight = 'deg'

        # convert the widths and bounds to astropy angles
        self.widthsToAngles()

        # a few housekeeping quantities for the Vizier query
        self.vizRowLimit = np.copy(maxRows)

        # which catalog?
        self.srcCat=srcCat[:]

        # redo the query even if the output file already exists?
        self.forceQuery = forceQuery

        # For this object - print to screen much?
        self.Verbose=Verbose

        # temporary file name
        self.filCat = '_tmpQuery.fits'

        # table holding the queried catalog
        self.tCat = Table()

        # do a few bits of checking on setup
        if runOnInit:
            self.checkBoxBounds()
            self.setupFieldCenter()
            self.convertCenToDegrees()

    def widthsToAngles(self):

        """Convert the given widths to astropy Angles, both for the
        values given and the bounds"""

        self.widAng = Angle(self.boxWidth, self.boxUnits)
        self.hgtAng = Angle(self.boxHeight, self.boxUnits)

        self.widAngMax = Angle(self.maxBoxWidth, self.maxUnitWidth)
        self.hgtAngMax = Angle(self.maxBoxHeight, self.maxUnitHeight)

    def checkBoxBounds(self):

        """Shrinks the box width to the max if greater."""

        # we could probably refactor this into a single method to call
        # twice, but at the moment I like having the screen output
        # different.

        if np.bool(self.widAng > self.widAngMax):
            if self.Verbose:
                print("GaiaCat INFO - replacing width %.2f %s w/ max %.2f %s" \
                    % (self.widAng.value, str(self.widAng.unit), \
                           self.widAngMax.value, str(self.widAngMax.unit)))
            self.widAng = self.widAngMax

        if np.bool(self.hgtAng > self.hgtAngMax):
            if self.Verbose:
                print("GaiaCat INFO - replacing height %.2f %s w/ max %.2f %s" \
                    % (self.hgtAng.value, str(self.hgtAng.unit), \
                           self.hgtAngMax.value, str(self.hgtAngMax.unit)))
            self.hgtAng = self.hgtAngMax

    def setupFieldCenter(self):

        """Sets up the field center, dependent on the datatype of the
        coordinates given."""

        # I trust this method more than simply putting this all into a
        # try/except clause.

        # If this was passed a string, astropy.SkyCoord can use
        # defaults. That's the easy case. If not, we have to tell
        # SkyCoord what units to use (I think!!). This is usually
        # going to be degrees, so let's stick with that for now.
        if not isinstance(self.ra, basestring):
            self.fieldCen = coord.SkyCoord(self.ra, self.dec, \
                                               unit=self.cenUnit, \
                                               frame=self.cooFrame)
        else:
            self.fieldCen = coord.SkyCoord(self.ra, self.dec, \
                                               frame=self.cooFrame)

    def convertCenToDegrees(self):

        """Sets the coordinate center in degrees from the field center
        object"""

        self.raDeg = np.copy(self.fieldCen.icrs.ra.degree)
        self.decDeg = np.copy(self.fieldCen.icrs.dec.degree)

    def doQueryWithCheck(self):

        """Does the query only if the catalog file doesn't already exist"""

        # This is really a wrapper for other routines.

        self.setOutFilename(informative=True)

        if os.access(self.filCat, os.R_OK) and not self.forceQuery:
            print "GaiaCat.doQueryWithCheck INFO - catalog %s already exists and forceQuery=False. Loading the old table instead of re-querying. (Set forceQuery=True to override.)"

            self.tCat = Table.read(self.filCat)
            return

        self.doQuery()
        self.writeCatToDisk()

        if self.Verbose:
            if len(self.tCat) >= Vizier.ROW_LIMIT:
                print "GaiaCat.doQuery WARN - returned table at maxrows. Consider increasing maxRows above %i" % (Vizier.ROW_LIMIT)

    def doQuery(self):

        """Does the query"""

        Vizier.ROW_LIMIT = int(self.vizRowLimit)

        t0 = time.time()
        if self.Verbose:
            print("testGaia.GaiaCat.doQuery INFO - starting query of %s..." \
                % (self.srcCat))

        # no, the flipping of "height" and "width" below is NOT a bug
        # - Vizier seems to flip these by default...
        result = Vizier.query_region(self.fieldCen, \
                                         width=self.hgtAng, \
                                         height=self.widAng, \
                                         catalog=self.srcCat)

        # this produces a one-element list... we write the result to
        # the tmp query
        if len(result) < 1:
            if self.Verbose:
                print "testGaia.GaiaCat.doQuery WARN - no results"
            return

        if self.Verbose:
            t1 = time.time()
            print "testGaia.GaiaCat.doQuery INFO - got %i rows in %.2e sec" \
                % (len(result[0]), t1 - t0)


        # pass to the table...
        self.tCat = result[0]

        # ... and assign metadata
        self.assignMetadata()

    def assignMetadata(self):

        """Sets up the table metadata with query parameters"""
        
        # types specified due to astropy.Table odd behavior for these
        # quantities:

        # query information
        self.tCat.meta['lonQ'] = float(self.ra)
        self.tCat.meta['latQ'] = float(self.dec)
        self.tCat.meta['frameQ'] = str(self.cooFrame) 
        self.tCat.meta['unitQ'] = str(self.cenUnit)

        # field center information
        self.tCat.meta['fieldRA'] = float(self.raDeg)
        self.tCat.meta['fieldDE'] = float(self.decDeg)

        # box information, in arcmin
        self.tCat.meta['boxwid'] = float(self.widAng.arcmin)
        self.tCat.meta['boxhgt'] = float(self.hgtAng.arcmin)
        self.tCat.meta['boxunit'] = str('arcmin')

        # maxrows
        self.tCat.meta['maxrows'] = int(self.vizRowLimit)

        # just for info, the table length
        self.tCat.meta['nresult'] = len(self.tCat)

    def setOutFilename(self, informative=False):

        """Generates output filename"""

        if not informative:
            self.filCat='_tmpQuery.fits'
            return

        # otherwise, set a temporary filename from the coordinates.
        sRa = '%.1f' % (self.raDeg)
        sDe = '%.1f' % (self.decDeg)
        sWi = '%.1f' % (self.widAng.arcmin)  # do in arcmin
        sHg = '%.1f' % (self.hgtAng.arcmin)

        self.filCat = '_tmp_%s_%s_%sx%s.fits' \
            % (sRa.zfill(5), sDe.zfill(5), sWi, sHg)

    def writeCatToDisk(self):

        """Writes the queried table to disk. Does not output a blank
        table."""

        if len(self.tCat) < 1:
            if self.Verbose:
                print "testGaia.GaiaCat.doQuery WARN - zero-length catalog returned."
            return

        # 2018-03-30 put the explicit existence-test here because the
        # overwrite argument changed between two astropy versions that
        # are different on machines at UM-D.
        try:
            self.wipeCatFil()
            self.tCat.write(self.filCat)
        except:
            print "problem writing out to file %s" % (self.filCat)

    def wipeCatFil(self):

        """Removes the output catalog file if it exists."""

        if os.access(self.filCat, os.W_OK):
            os.remove(self.filCat)

def queryGaia(ra=2.65, dec=-3.25, srcCat='gaia', \
                  maxRows=100000, \
                  wid=5., hgt=5., \
                  boxUnits='arcmin', \
                  frame='galactic', \
                  wipeAfter=False, Verbose=True):

    """Wrapper to query the Gaia source catalog. 

    RETURNS an astropy table with the Gaia query about the specified region. 

    Examples of usage:
    =================

    # galactic coordinates, default widths in arcmin
    tDum = testGaia.queryGaia(2.65, -3.25)

    # equatorial coordinates, width in degrees, short table for debug
    tDum = testGaia.queryGaia(271.1, -28.3, wid=0.1, hgt=0.1, boxUnits='deg', frame='fk5', maxRows=1000)

    # no screen output, output NOT written to disk:
    tDum = testGaia.queryGaia(2.65, -3.25, wipeAfter=False, Verbose=False)

    """

    # The object is written to do the coordinate conversion on
    # initialization. To test those individual pieces, set
    # runOnInit=False in the call below.
    GC = GaiaCat(ra, dec, wid, hgt, frame, maxRows=maxRows, \
                     boxUnits=boxUnits, \
                     srcCat=srcCat, Verbose=Verbose)
    
    # now do the query - UPDATE: replaced with a wrapper
    GC.doQueryWithCheck()

    # remove the temporary file?
    if wipeAfter:
        GC.wipeCatFil()

    # return the Gaia catalog as an astropy table
    return GC.tCat


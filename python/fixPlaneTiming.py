#
# fixPlaneTiming.py
#

# WIC 2018-02-17

# Hack to correct the time keywords in v404 Cyg single-plane fits
# files (there was an error converting to JD when we split datacubes
# into planes on the mountain).

# Note to self: do this procedurally (rather than object oriented) to
# be easier to follow.

# We could break this up into multiple methods for flexibility with
# other datasets. For the moment, however, let's confine this to one
# or two monolithic methods to make the thread easier to follow. 

# for the wrapper
import glob

from astropy.io import fits
import os

# for the barycentric correction
from astropy import time, coordinates as coord, units as u

def fixAll(srchString='tmp*_proc_0??.fits*', srchDir='./', \
               strAvoid='TFIX'):

    """Finds fits files in the current directory and corrects
    them. Example call:

    fixPlaneTiming.fixAll('tmp_*proc_???.fits*', './', 'TFIX')

    The default search string tmp*_proc_???.fits* should avoid all the
    actually processed files.

    srchDir = directory to search for input files

    strAvoid = ignore files containing this string (useful if you
    don't want to re-correct output found in the input directory)"""

    lToFix = glob.glob('%s/%s' % (srchDir, srchString))
    
    # correct all of the files
    iDone = 0 
    for pathIn in lToFix:

        if pathIn.find(strAvoid) > -1:
            continue

        sStatus = fixHeader(pathIn)

        if len(sStatus) < 1:
            print "fixPlaneTiming.fixAll WARN - problem with file %s" \
                % (pathIn)
        else:
            iDone = iDone + 1

    print "Corrected %i of %i frames in matched-list" \
        % (iDone, len(lToFix))

def fixHeader(pathIn='./tmp_V404Cyg_52_proc_001.fits', \
                  iPlane=-1, \
                  pathOut='', stemOut='TFIX', dirOut='./', \
                  moveBackup=True, dirBackup='./UNCORRECTED', \
                  keyStart='HJD', keyEnd='HJDEND', keyExptime='KCT', \
                  keyJD='JD', dtOverhead=0.0, \
                  keyOverhead='', \
                  keyRA='RA', keyDEC='DEC', \
                  keyPlaneOut = 'IPLANE', \
                  stemCube='c', \
                  extRead=0, \
                  whichObservatory='mdm', \
                  Verbose=True):

    """Corrects the header for an MDM tmp_plane.fits file. Written to
    be used conveniently AFTER the files have been moved into
    place. Returns a string giving the name of the output file if
    successful, a blank string if not.

    Example call:
    sOut = fixPlaneTiming.fixHeader('./tmp_V404Cyg_52_proc_001.fits.gz', moveBackup=True)

    ARGUMENTS:
    ==========

    pathIn = path to input fits file. Can be gzipped.

    iPlane = which plane this is (counting from 1). If <0, will work
    this out from the filename.

    pathOut = path for corrected FITS file. If blank, will be
    constructed from the input filename and the 'stemOut' argument

    stemOut = string to append to the output filename

    dirOut = output directory 

    moveBackup = True -- move the original file to a backup subdirectory?

    dirBackup = subdirectory to move the uncorrected files. If not
    present, will be created.

    keyStart, keyExptime, keyEnd are the header keywords for the time
    at start, the cycle time, and the end HDD (the latter is read in
    case we want to do some sanity checking later).

    keyJD = header keyword for the un-barycenter corrected julian date

    dtOverheads = overhead time per cycle to add to the keyExptime to
    get cycle time (default 0.0s)

    keyOverhead = header keyword for any overhead time to add to the
    cycle time. Defaults to blank (no overhead to add)

    keyRA, keyDEC = header keywords for object co-ordinates (needed
    for the barycenter correction)

    keyPlane = keyword in the output plane for which plane we decided
    this was

    stemCube = 'c' - if we want to preserve the HJD, HJDEND but for
    the datacube, create a new header keyword with this stem at the
    end (so the default would be 'HJDC' and 'HJDENDC'). Note that FITS
    header keywords have an 8-character limit. 

    extRead = FITS extension of the input file to read. Defaults to 0
    (the first extension).

    whichObservatory = 'mdm' : observatory location for HJD correction
    (to list options, coord.EarthLocation.get_site_names())

    Verbose = True - print screen output?

    """

    # return status string
    statusRet = ''

    # is the input file readable?
    if not os.access(pathIn, os.R_OK):
        if Verbose:
            print "fixPlaneTiming.fixHeader FATAL - cannot read path %s" \
                % (pathIn)

        return statusRet

    # get the file and the header
    try:
        aImg, header = fits.getdata(pathIn, extRead, header=True)
    except:
        if Verbose:
            print "fixPlaneTiming.fixHeader FATAL - problem reading %s extension %i" % (pathIn, extRead)
        return statusRet

    # read the HJD keywords. For the moment we hardcode this to assignments.

    # DAY NUMBER OF EXPOSURE START
    try:
        hjdStart = header[keyStart]
    except:
        print "fixPlaneTiming.fixHeader FATAL - keyword %s not present in %s:%i" % (keyStart, pathIn, extRead)
        return statusRet

    # EXPOSURE TIME
    try:
        expTime = header[keyExptime]
    except:
        print "fixPlaneTiming.fixHeader FATAL - keyword %s not present in %s:%i" % (keyExptime, pathIn, extRead)
        return statusRet

    # OVERHEAD FROM FILE
    dtOverheadHeader = 0.
    if len(keyOverhead) > 0:
        try:
            dtOverheadHeader = header[keyOverhead]
        except:
            if Verbose:
                print "fixPlaneTiming.fixHeader WARN - overhead keyword %s not present in %s:%i" % (keyOverhead, pathIn, extRead)

    # EXPOSURE END-TIME FROM HEADER
    hjdEnd = -99.9  # easy-to-find default
    if len(keyEnd) > 0:
        try:
            hjdEnd = header[keyEnd]
        except:
            if Verbose:
                print "fixPlaneTiming.fixHeader WARN - exposure-end keyword %s not present in %s:%i" % (keyEnd, pathIn, extRead)
        
    # Compute the time interval between the start of each cycle
    tCycle = expTime + dtOverheadHeader + dtOverhead

    # Which plane of the cube was this?
    if iPlane > -1:
        iFound = iPlane
    else:
        iFound = scrapePlaneFromFilename(pathIn)

    if iFound < 0:
        if Verbose:
            print "fixPlaneTiming.fixHeader WARN - bad plane number: %i" \
                % (iFound)
        return statusRet

    # compute the start and end time of the exposure
    dayStart = hjdStart + ( float(iFound) * tCycle ) / 86400.
    dayEnd = dayStart + tCycle / 86400.

    # create the cube header keywords.
    keyStartCube = appendKeyword(keyStart, stemCube)
    keyEndCube = appendKeyword(keyEnd, stemCube)

    # Now we remove the HJD correction. This isn't going to be quite
    # right since we're given the HJD and want to go to JD - so we'll
    # be wrong about where the Earth was at the time of observation by
    # a minute or two - but this shouldn't make much of a difference
    # compared to the Earth's orbit (500,000 minutes or so).

    # Note that we can only do this if we can read the
    # coordinates. Rather than send lots of arguments to another
    # method, for the moment we'll do everything here...
    sRA = ''
    sDE = ''
    doHJDtoJD = True

    jdStart = -99.9  # obvious default

    try:
        sRA = header[keyRA]
        sDE = header[keyDEC]
    except:
        doHJDtoJD = False

    didBary = False  # flag: did we do the barycentric correction?

    if doHJDtoJD:
        targLocation = coord.SkyCoord(sRA, sDE, \
                                          unit=(u.hourangle, u.deg), \
                                          frame='icrs')
        siteLocation = coord.EarthLocation.of_site('mdm')

        # populate the astropy time object to compute the heliocentric
        # correction. Note that we don't have the photon arrival time
        # at the observatory in this file! (We might later decide to
        # re-produce these from the datacube. Let's cross that bridge
        # if we come to it.)
        timeObj = time.Time(dayStart, format='jd', scale='utc', \
                                location=siteLocation)

        # estimate the barycentric correction at this location. Note
        # there's an inconsistency between the keyword name ("HJD")
        # and the comment ("barycentric"). Trust the comment for now.

        # 2018-02-22 - if the machine uses an older astropy version,
        # the light_travel_time() method might not be present. If so,
        # exit gracefully
        if hasattr(timeObj,'light_travel_time'):
            ltt_bary = timeObj.light_travel_time(targLocation, 'barycentric')

            # we SUBTRACT this - in days - off the HJDStart to get the JDstart
            timeObj_atObs_approx = timeObj - ltt_bary

            jdStart = timeObj_atObs_approx.jd
            didBary = True

        else:
            jdStart = np.copy(dayStart)

    # now update the FITS header. First ensure the datacube parameters
    # are written
    header[keyStartCube] = hjdStart
    if hjdEnd > 0:
        header[keyEndCube] = hjdEnd
    
    # now update the HJD and HJDEND keywords
    header[keyStart] = dayStart
    header[keyEnd] = dayEnd
    header[keyPlaneOut] = iFound

    header['didBARY'] = ( int(didBary), 'HJD barycentric un-done for JD' )
    

    # If we made it, update the JD keyword
    if jdStart > 0:
        header[keyJD] = jdStart


    # If we've got here, we can write output. So - ensure the output
    # path is set properly, and move the original file to backup
    # subdirectory if desired.
    if len(pathOut) < 2:

        # hack for gzipped input
        isGZ = os.path.splitext(pathIn)[-1].find('gz')
        splitIn = os.path.splitext(pathIn.split('.gz')[0])

        stemIn = os.path.split(splitIn[0])[-1]

        pathOut = '%s/%s_%s%s' % (dirOut,stemIn, stemOut, splitIn[-1])

        if isGZ > -1:
            pathOut = '%s.gz' % (pathOut)

    # by this point we can write output!
    fits.writeto(pathOut, aImg, header, overwrite=True, output_verify="fix")

    # ensure the output directory is present
    if moveBackup:
        if len(dirBackup) > 1:
            if not os.access(dirBackup, os.R_OK):
                os.makedirs(dirBackup)

            # output path
            stemBak = os.path.split(pathIn)[-1]
            pathBak = '%s/%s' % (dirBackup, stemBak)

            os.rename(pathIn, pathBak)

    return pathOut

def scrapePlaneFromFilename(pathIn='', separator='_', position=-1, \
                             inputCountsFromOne=True, \
                             outputCountsFromOne=False, \
                             Verbose=True):

    """Scrapes the filename to return the plane number. Uses a simple
    heuristic: the filename is split by the 'separator' character and
    the character at position 'position' in the resulting array is
    evaluated as an integer.

    inputCountsFromOne = does the INPUT name count from one?

    outputCountsFromOne = should the OUTPUT count from zero?

    Defaults to -1 (i.e. no plane)"""

    whichPlane = -1
    
    # split the input string into an array by separator
    lSplit = os.path.splitext(pathIn)[0].split(separator) 
    
    try:
        # hack to deal with the .fits.gz case
        iFound = int(lSplit[position].split('.f')[0])
    except:
        if Verbose:
            print "fixPlaneTiming.scrapePlaneFromFilename WARN - problem extracting position %i from list length %i" % (position, len(lSplit))
            print lSplit

    # correct for the convention used coming in vs going out
    if not outputCountsFromOne:
        if inputCountsFromOne:
            iFound -= 1

    if outputCountsFromOne:
        if not inputCountsFromOne:
            iFound += 1

    return iFound

def appendKeyword(keyOld='HJD', keyAdd='c', maxLen=8):

    """Utility to produce a new header keyword that respects the FITS
    convention's maximum length convention.

    keyOld = old keyword

    keyAdd = string to append

    maxLen = maximum string length (don't change from 8 unless you
    know what you are doing)"""

    keyRet = '%s%s' % (keyOld[0:maxLen-len(keyAdd)], keyAdd)
    return keyRet[0:maxLen]
    


#
# csv2fits.py
#

# WIC 2019-08-27: simple utility to convert .csv.gz files into fits
# format for ingestion into lsd

import glob, os, sys
from astropy.table import Table
from astropy.time import Time

def go(dirCSV='./csv', dirFITS='./fits', srchStr='GaiaSource*csv*', \
           filLog='csv2fits.log', gzOutput=True, showProgress=True):

    """Convert csv to fits via astropy. Assumes the .csv[.gz] files
    are in subdirectory [dirCSV] and the output fits files will go
    into a different subdirectory [dirFITS].

    srchStr -- search string to use by glob when finding files in a
    given directory.

    gzOutput -- write output as .fits.gz instead of .fits

    showProgess -- prints progress to a single line in the terminal."""

    # Designed to run from the following location on Rigel:
    # /raid1/data/gaiadr2/cdn.gea.esac.esa.int/Gaia/gdr2/gaia_source

    # Much of the code here is actually boilerplate bells and
    # whistles, since this is going to be run on half a terabyte of
    # .csv.gz files and I want to have some sort of record for
    # progress.

    # return rather than try to guess user intention
    if not os.access(dirCSV, os.R_OK):
        print("csv2fits.go WARN - cannot read input directory %s" \
                  % (dirCSV))
        return

    # look for source files. Trust the user to invent a sensible
    # search string.
    searchString = '%s/%s' % (dirCSV, srchStr)
    lFiles = glob.glob(searchString)

    if len(lFiles) < 1:
        print("csv2fits.go WARN - no files match %s" % (searchString))

    if len(dirFITS) < 3:
        print("csv2fits.go WARN - dirFits < 3 characters. Not proceeding.")
        return

    # guard against some stupid logfile name like './'
    if len(filLog) < 3:
        filLog = 'csv2fits.log'

    with open(filLog, 'w') as wObj:

        tStarted = Time.now()
        sDate = '%s (UT)' % (tStarted.datetime)
        wObj.write('# csv2fits started at %s\n' % (sDate))
        wObj.write('# From directory %s\n' % (os.getcwd()))

        # ok IF we got here, then we have our list of files to
        # process. Proceed:
        for iFil in range(len(lFiles)):
            thisFil = lFiles[iFil]
        
            lStem = os.path.splitext(os.path.split(thisFil)[-1])

            # if the resulting extension has bz, gz, bz2, zip,
            # ... anything with 'z' in it, then it's compressed and we
            # need to split again.
            if lStem[-1].find('z') > -1:
                fStem = os.path.splitext(lStem[0])[0]
            else:
                fStem = lStem[0]

            thisFits = '%s.fits' % (fStem)
            if gzOutput:
                thisFits = '%s.gz' % (thisFits)

            outPath = '%s/%s' % (dirFITS, thisFits)

            if showProgress:
                sys.stdout.write("\r csv2fits: %4i/%-4i: %s --> %s" \
                                     % (iFil+1, len(lFiles), thisFil, outPath))
                sys.stdout.flush()

            # we wrap the reading and writing into a try/except clause
            # so that we can handle bad cases (we have tens of
            # thousands of input files)
            try:

                # now we read the input frame
                tThis = Table.read(thisFil, format='ascii.csv')
            
                # only if this was successful do we proceed
                if not os.access(dirFITS, os.R_OK):
                    os.makedirs(dirFITS)
            
                tThis.write(outPath, overwrite=True)

                wObj.write('%s ok\n' % (thisFil))

            # break out of the loop if a keyboard interrupt is
            # detected. (Here, the counter will assume the current
            # file failed.)
            except KeyboardInterrupt:
                print("csv2fits.go INFO - keyboard interrupt detected.")
                wObj.write('# Elapsed time before keyboard interrupt: %i files in %.2e seconds\n' \
                               % (iFil, (Time.now()-tStarted).sec))
                return


            except:                
                wObj.write('%s PROBLEM\n' % (thisFil))


        # let's record the elapsed time in our logfile
        wObj.write('# Elapsed time: %i files in %.2e seconds\n' \
                       % (iFil + 1, (Time.now()-tStarted).sec))

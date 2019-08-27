#
# csv2fits.py
#

# WIC 2019-08-27: simple utility to convert .csv.gz files into fits
# format for ingestion into lsd

import glob, os, sys
from astropy.table import Table
from astropy.time import Time

def go(dirCSV='./csv', dirFITS='./fits', srchStr='GaiaSource', \
           filLog='csv2fits.log', gzOutput=True, showProgress=True):

    """Convert csv to fits via astropy"""

    # return rather than try to guess user intention
    if not os.access(dirCSV, os.R_OK):
        print("csv2fits.go WARN - cannot read input directory %s" \
                  % (dirCSV))
        return

    # look for source files
    searchString = '%s/%s*csv*' % (dirCSV, srchStr)
    lFiles = glob.glob(searchString)

    if len(lFiles) < 1:
        print("csv2fits.go WARN - no files match %s" % (searchString))

    if len(dirFITS) < 3:
        print("csv2fits.go WARN - dirFits < 3 characters. Not proceeding.")
        return

    # guard against some stupid logfile name like './'
    if len(filLog) < 3:
        filLog = 'csv2fits.log'

    # open the logfile for writing
    #if os.access(filLog, os.R_OK):
    #    os.remove(filLog)

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
            # so that we can handle bad cases (we have hundreds of
            # input files)
            try:

                # now we read the input frame
                tThis = Table.read(thisFil, format='ascii.csv')
            
                # only if this was successful do we proceed
                if not os.access(dirFITS, os.R_OK):
                    os.makedirs(dirFITS)
            
                tThis.write(outPath, overwrite=True)

                wObj.write('%s ok\n' % (thisFil))

            except:
                
                wObj.write('%s PROBLEM\n' % (thisFil))

        # let's record the elapsed time in our logfile
        wObj.write('# Elapsed time: %i files in %.2e seconds\n' \
                       % (iFil + 1, (Time.now()-tStarted).sec))

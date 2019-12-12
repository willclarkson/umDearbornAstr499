#
# doAnim.py
#

# 2019-12-12 WIC & ZCS - call the animation maker using several images

import programday_4
import glob
import os

def doMany(srchString='Comet_R',catTail='ASC', Verbose=True):

    """Loops through the .asc catalogs and runs the programday4"""

    # first, find the .asc catalogs
    lCats = glob.glob('%s*.%s' % (srchString, catTail))

    # now run programday4 on each file.
    for iCat in range(len(lCats)):

        thisCat = lCats[iCat]
        
        # look for the corresponding image:
        stemCat = thisCat.split('.%s' % (catTail))[0]
        thisImg = '%s.fits' % (stemCat)

        if not os.access(thisImg, os.R_OK):
            print("doAnim.doMany WARN - cannot read image file %s" \
                  % (thisImg))
            return

        # now we're here, we can run "go"
        programday_4.go(thisCat, thisImg)
        

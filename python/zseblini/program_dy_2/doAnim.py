#
# Doanim.Py
#

# 2019-12-12 Wic& ZCS - call the animation maker using several images

import programday_4
import glob
import os
import subprocess
from astropy.table import Table, vstack

def doMany(srchString='Comet_R',catTail='ASC', Verbose=True):

    """Loops through the .asc catalogs and runs the programday4"""

    # first, find the .asc catalogs
    pattn='%s*.%s' % (srchString, catTail)
    lCats = glob.glob(pattn)
    MasterTable=Table()

    if len(lCats) < 1:
        print("doAnim.doMany WARN - no files match search string %s" % (pattn))
    
    # now run programday4 on each file.
    for iCat in range(len(lCats)):

        thisCat = lCats[iCat]
        
        # look for the corresponding image:
        stemCat = thisCat.split('.%s' % (catTail))[0]
        thisImg = '%s.fits' % (stemCat)

        if not os.access(thisImg, os.R_OK):
            if Verbose:
                print("doAnim.doMany WARN - cannot read image file %s" \
                      % (thisImg))
            continue

        # now we're here, we can run "go"
        ThisTable=programday_4.go(thisCat, thisImg)
        MasterTable=vstack([MasterTable,ThisTable])
        #print(ThisTable)
        # stack the table

        
       # thisMov='%s.mp4' % (srchString)
        
    os.system("ffmpeg -r 1 -i *.jpg -vcodec mpeg4 -y 'thisMov'.mp4")
    ffmpeg.view()
    #"ffmpeg.view")
    print(MasterTable)

#
# Doanim.Py
#

# 2019-12-12 Wic& ZCS - call the animation maker using several images

import programday_4
import getPositions
import glob
import os
import subprocess
from astropy.table import Table, vstack
from astropy.samp import SAMPIntegratedClient
#from astropy.io.votable.tree import VOTableFile, Resource, Table, Field

def doMany(srchString='Comet_R',catTail='ASC', Verbose=True, clobberCats=False, \
           renameFrames=True, tmpDirFrames='tmpFrames', tmpFrameStem='tmpFrame_'):

    """Loops through the .asc catalogs and runs the programday4

    Argument 'renameFrames' causes the routine to ensure that the frame jpg's are numbered sequentially. This means ffmpeg will know where to find them [its 'glob' pattern matching does not appear to work]."""

    # first, find the .asc catalogs
    pattn='%s*.%s' % (srchString, catTail)
    lCats = glob.glob(pattn)

    # how about this: if the length of lCats is zero, we call getPositions.go to ensure it DOES exist??
    if len(lCats) < 1 or clobberCats:
        getPositions.go(srchString,catTail, Verbose)
    print(lCats)
    
    MasterTable=Table()

    if len(lCats) < 1:
        print("doAnim.doMany WARN - no files match search string %s" % (pattn))

    # output subdirectory for the frames
    if len(tmpDirFrames) < 3:
        tmpDirFrames = os.getcwd()
    else:
        if not os.access(tmpDirFrames, os.R_OK):
            os.makedirs(tmpDirFrames)
            
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

        # generate the filename for the image frame
        thisFrameImg = ''
        if renameFrames:
            thisFrameImg = '%s/%s%s.jpg' % (tmpDirFrames, tmpFrameStem,str(iCat).zfill(3))

        ThisTable=programday_4.go(thisCat, thisImg, frameName=thisFrameImg)
        MasterTable=vstack([MasterTable,ThisTable])
        print(ThisTable)
        # stack the table

        
       # thisMov='%s.mp4' % (srchString)

    #os.system("ffmpeg -r 1 -i *.jpg -vcodec mpeg4 -y 'thisMov'.mp4")
    #print(MasterTable)

    # 2019-12-18 - updated the call to ffmpeg to use the temporary
    # frame names that we forced, so that we can tell ffmpeg where to
    # find them below. Since python doesn't like strings with '%' in
    # them, we do this in two steps:
    sSrch = tmpDirFrames+'/'+tmpFrameStem+'%03d.jpg'
    os.system("ffmpeg -r 1 -i %s -vcodec mpeg4 -y 'tmpMov'.mp4" % (sSrch))
    
    # let's save to votable - see if this works...
    MasterTable.write('MasterTable.xml',format='votable',overwrite=True)
    
    # asking users if we should save the master table
    
    ##Confirmation= input( "Save MasterTable into a .xml file? y or n?")
    ###if Confirmation == 'y':
    ###    votable.to_xml('MasterTable.xml')
    ##else:
    ###     return
     
        

    #client = SAMPIntegratedClient()
    #client.connect()
    #params = {}
    #params["url"] = 'file:///Users/zseblini/projects/umDearbornAstr499/python/zseblini/program_dy_2/MasterTable.xml'
    #params["name"] = "Seblini et al. (2019), Table 1"
    #message = {}
    #message["samp.mtype"] = "table.load.votable"
    #message["samp.params"] = params
    
    

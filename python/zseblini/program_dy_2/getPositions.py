#getPositions
#
# 2019-12-12 WIC & ZIS - use source extractor to get photometry of a
# list of images, where we control the filename

import glob
import os, shutil
import subprocess
import inspect

def testInspect():

    """Try the inspect module to learn where this routine is stored"""

    myPath = inspect.getfile(inspect.currentframe())
    myDirec = os.path.split(myPath)[0]

    return myDirec
    
    
def go(srchString='Comet_R', catTail='ASC', Verbose=True, \
       filPars='se_umdobs.param', filConfig='se_umdobs.se'):

    """Use source extractor to extract positions for a list of images"""

    # First, let's find the images
    lImgs = glob.glob('*%s*.fit?' % (srchString))

    # If there are no images, there's nothing to do
    if len(lImgs) < 1:
        if Verbose:
            print("getPosistions.go WARN - no images found")
        return

    # let's try constructing dirTop using the source for this module:
    dirTop = testInspect()
    if len(dirTop) < 1:
        dirTop = os.getcwd()
#    dirTop = '/Users/zseblini/projects/umDearbornAstr499/python/zseblini/program_dy_2'
    dirPars = 'params'

    # source path for the se and param file:
    pathSrcPars = '%s/%s/%s' % (dirTop, dirPars, filPars)
    pathSrcConf = '%s/%s/%s' % (dirTop, dirPars, filConfig)

    print("getPositions.go INFO - path for source params:", pathSrcPars, os.access(pathSrcPars, os.R_OK))
    if os.access(pathSrcPars, os.R_OK):
        shutil.copy2(pathSrcPars, os.getcwd())

    if os.access(pathSrcConf, os.R_OK):
        shutil.copy2(pathSrcConf, os.getcwd())

    # we could put in some sort of condition-trap here to see if the
    # files actually came across. For the moment we'll press on.

    # let's go through the images, and lift various bits of the
    # filename that we want
    for iImg in range(len(lImgs)):

        # extract a substring to identify this particular image
        thisImg = lImgs[iImg].split('.fit')[0]

        # insert this into various files we're going to want
        thisCat = '%s.%s' % (thisImg, catTail)

        # now we construct the command to call the extractor. We build
        # this argument-by argument, like so:
        commandProg = '/usr/local/bin/sex'
        argsProg = ' %s -c %s' \
                   % (lImgs[iImg], filConfig)

        # Now we split this into an input subprocess can actually
        # understand...
        argsList = argsProg.split()
        commandRun = [commandProg]+argsList
        returnCode = subprocess.call(commandRun)
        
        # after running, we move the output catalog to the name we wanted.
        if os.access('SE_UMDOBS.ASC', os.R_OK):
            shutil.move('SE_UMDOBS.ASC', thisCat)
        else:
            if Verbose:
                print("getPositions.go WARN - no extractor output for %s" \
                      % (lImgs[iImg]))
            

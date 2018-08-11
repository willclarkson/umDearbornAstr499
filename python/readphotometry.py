import sys
import os
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
plt.ion() # ensure we get the command line back

# the glob library contains a function, also called glob, that finds files and directories whose names match a pattern. We provide those patterns as strings: the character *
# matches zero or more characters, while ? matches any one character. We can use this to get the names of all files in the directory.
import glob

#from astropy.table import Table
#from astropy.io import fits

# for serializing to disk (will still be quicker to read than to re-generate)
import cPickle as pickle

# a test method to read-in one photometry file

def testread():
    #try:
        #image = np.genfromtxt('tmp_V404Cyg_51_proc_001.fits.flux')
        #images = np.genfromtxt('tmp_*.flux')
    image = np.genfromtxt('tmp_V404Cyg_51_proc_002_TFIX.fits.flux')
    print(image)
    plt.imshow(image)
    np.shape(image)
    image[0]
            #except:
#print('PROBLEM IMPORTING IMAGE')
        

#print('IMAGE READ-IN SUCCESSFUL.')
    return


# ^^^^ This works. ^^^^


def teststack():
    a=np.arange(1,29)
    filenames=sorted(glob.glob('tmp_*.flux'))
    filenames = filenames[0:]
    for f in filenames:
        data=np.genfromtxt(fname=f)
        np.hstack((a,data))

# ignore ^^^^


# read-in all photometry files
def read(prefix='',Verbose=False):
    
   
    hjdAll = np.array([])
    fluxAll = np.array([])
    unctAll = np.array([])

    print("read INFO - searching for files %s\*.flux ..." % (prefix))

    Lfound = glob.glob('%s*.flux' % (prefix))

    if Verbose:
        print("read INFO - file list with prefix %s :" % (prefix),
              Lfound)
            
    # filenames = sorted(glob.glob('%s*.flux' % (prefix)))
    filenames = sorted(Lfound)
    filenames = filenames[0:]
    for f in filenames:

# get the filename for the image

        # If error message, "IOError: [Errno 2] No such file or directory: '../ref'" shows up, change pathImg to either go up a directory
        # (pathImg = '../%s' ... ) or look in the current directory (pathImg = '%s' ... )

        #pathImg = '../%s' % (f.split('.flux')[0]) # use when in a folder inside Output_dir ../ goes up a directory
        pathImg = '../%s' % (f.split('.flux')[0]) # use when in Output_dir
        sys.stdout.write('\r %s, %s, %s' % (f, pathImg, os.access(pathImg, os.R_OK)))
        sys.stdout.flush()

        data = np.genfromtxt(fname=f)
    #print(data)
    #print np.shape(data)
        
        flux = data[:,0]
        error = data[:,1]
        #print(np.shape(flux))
        #print(np.shape(error))
    
        # if we don't yet have the flux array, copy in
        if np.size(fluxAll) < 1:
            fluxAll = np.copy(flux)
            unctAll = np.copy(error)
        else:
            # stuff comes here
            if Verbose:
                print('all, single:', 'shape of fluxAll array:', np.shape(fluxAll),'shape of flux array:', np.shape(flux))

            # fluxAll = STACK((fluxAll, flux))
            fluxAll = np.vstack((fluxAll, flux))    # gives ('all single:', (1xx, 58), (58))
            unctAll = np.vstack((unctAll, error))
            #fluxAll = np.hstack((fluxAll, flux))    # gives ('all single:', (xxxx,), (58,))
            #unctAll = np.hstack((unctAll, error))

#        print(flux)
#        print(error)

        # OK that gets us the data. Now let's read the HJD
        #print(" ")
        #print("FITS IMAGE:",pathImg)
        #print(" ")
        thisHdr = fits.getheader(pathImg)
        thisHJD = thisHdr['HJD']
        hjdAll = np.hstack((hjdAll, thisHJD))

    return hjdAll, fluxAll, unctAll

def loadAndPlot(iShow=[0], inFile='lightcurves.cPickle', convertTimes=True):

    """Imports the set of lightcurves from disk and plots them """

    try:
        Dcurves = pickle.load(open(inFile,'r'))
    except:
        print("importAndPlot FATAL - problem reading input file %s" %(inFile))
        return

    # now we unpack the dictionary into the variables we expect (this is really so that we can paste our commands from below)
    times = Dcurves['times']
#fluxAll = Dcurves['fluxAll']
#    unctAll = Dcurves['unctAll']

    sLabelX = 'HJD (days)'

    if convertTimes:
        times = (times - np.min(times))*1440.
        sLabelX = 'Minutes elapsed'

    # create the figure when we need it
    fig = plt.figure(1)
    fig.clf()
    ax = fig.add_subplot(111)

    for i in range(len(iShow)):
        iThis = iShow[i]
        sLabel = 'Object %i' % (iThis+1)
        #ax.plot(x,fluxAll[:,iThis], label=sLabel)
        dum = ax.errorbar(times, Dcurves['fluxAll'][:,iThis], Dcurves['unctAll'][:,iThis],
                          label=sLabel, lw=1)

        ax.set_xlabel(sLabelX)
        #ax.set_ylim(-2000, 1000)
        leg = ax.legend()

    plt.show()


def read_and_plot(iShow = [0], convertTimes=True, prefix='', Verbose=False, \
                  outFile='lightcurves.cPickle'):
    
    # call our read routine to bring in the data
    hjdAll, fluxAll, unctAll = read(prefix=prefix, Verbose=Verbose)
    
    fig = plt.figure(1)
    fig.clf()
    ax = fig.add_subplot(111)

    times = np.copy(hjdAll)
    sLabelX = 'HJD (days)'
    
    if convertTimes:
        times = (hjdAll - np.min(hjdAll))*1440. # 7/17/18: running this method on aligned pyDIA data gives an error pointing to this line
        sLabelX = 'Minutes elapsed'

# print the shape of the big arrays we've built
    print("read_and_plot INFO: array shapes:", np.shape(times), np.shape(fluxAll), np.shape(unctAll))

# now we build the output
    DOut = {'times':times, 'fluxAll':fluxAll,'unctAll':unctAll}
    pickle.dump(DOut,open(outFile,'w'))

#x=np.arange(len(fluxAll))
    
    for i in range(len(iShow)):
        iThis = iShow[i]
        sLabel = 'Object %i' % (iThis+1)
        #ax.plot(x,fluxAll[:,iThis], label=sLabel)
        dum = ax.errorbar(times, fluxAll[:,iThis], unctAll[:,iThis], label=sLabel, lw=1)

    ax.set_xlabel(sLabelX)
    ax.set_ylim(-2000, 1000)
    leg = ax.legend()

    # shows the lightcurve plot
    plt.show()
    


    return fluxAll, unctAll
    
# read-in star positions
def starPos(degFitX=1):

    jdAll = np.array([])
    xShiftAll = np.array([])
    yShiftAll = np.array([])
    print('Reading in files...')
    
    files = glob.glob('Aligned_*.fits')
                 
#Data=np.genfromtxt

                    #fileNames = sorted(files)
                    #fileNames = fileNames[0:]

    for f in files:
                    #Data=np.genfromtxt(fname=f)
        thisHdr = fits.getheader(f)
        thisJD = thisHdr['JD']
        jdAll = np.hstack((jdAll, thisJD))
        thisX = thisHdr['X_SHIFT']
        thisY = thisHdr['Y_SHIFT']
        xShiftAll = np.hstack((thisX, xShiftAll))
        yShiftAll = np.hstack((thisY, yShiftAll))

        #xyShift = np.vstack((thisX, thisY))
    #print jdAll, xShiftAll, yShiftAll

#   update the time vector in-place to convert from days --> seconds
    jdAll = (jdAll - np.min(jdAll))*86400.

    # fit low-order polynomial to the shifts
    parsX = np.polyfit(jdAll, xShiftAll, deg=degFitX)
    parsY = np.polyfit(jdAll, yShiftAll, deg=degFitX)

    tFine = np.linspace(np.min(jdAll)-0.01, np.max(jdAll)+0.01, 1000, endpoint=True)


#Plot data
    figx = plt.figure(1)
    figx.clf()
    axx = figx.add_subplot(211)
    dum1 = axx.scatter(jdAll, xShiftAll, marker='^', s=9, c='b', \
                       label=r"xShift")
    
    if np.abs(degFitX - 1) < 1e-3:
        axx.set_title(r'$\Delta x = (%.2e)t + (%.2e$)' % (parsX[0], parsX[1]))
    #axx.set_xlabel(r"$\Delta t$, seconds")
    
    # overplot the fit
    dum2 = axx.plot(tFine, np.polyval(parsX, tFine), 'k-')
    
    #plt.show(block=False)
    
#figy = plt.figure(2)
#figy.clf()
    axy = figx.add_subplot(212, sharex=axx)
    dum2=axy.scatter(jdAll, yShiftAll, marker='s', s=9, c='r', \
                label=r"yShift")
    #plt.show(block=False)
    dum3 = axy.plot(tFine, np.polyval(parsY, tFine), 'k-')

    #axy.set_xlabel('JD - min(JD), seconds')
    sLabelT =r"$\Delta t$, seconds"
    axy.set_xlabel(sLabelT)
    #axy.set_xlabel(r"$\Delta t$, seconds")

#    for ax in [axx, axy]:
#        ax.grid(which='both')
#        leg=ax.legend()

    axx.set_ylabel('xShift (pix)')
    axy.set_ylabel('yShift (pix)')


    figxy = plt.figure(3)
    figxy.clf()
    axxy = figxy.add_subplot(111)
    dum = axxy.scatter(xShiftAll, yShiftAll, c=jdAll, marker='o', label='path', edgecolor='0.5')
    axxy.set_xlabel('xShift')
    axxy.set_ylabel('yShift')

    cbar = figxy.colorbar(dum)
    
    # let's get the residuals from the straight line fit in X and the fit in Y
    residX = xShiftAll - np.polyval(parsX, jdAll)
    residY = yShiftAll - np.polyval(parsY, jdAll)

    fig4 = plt.figure(4)
    fig4.clf()
    ax4 = fig4.add_subplot(211)
    dum4 = ax4.plot(jdAll, residX, color='b', ms=3, label='X Residual', marker='^')
#ax4.set_xlabel(r"JD (seconds)")
    ax4.set_ylabel(r"$\Delta x$, pix")
    ax5 = fig4.add_subplot(212, sharex=ax4)
    dum5 = ax5.plot(jdAll, residY, color='r', ms=3, label='Y Residual', marker='s')
    #ax5.set_xlabel(r"JD(seconds)")
    ax5.set_xlabel(sLabelT)
    ax5.set_ylabel(r"$\Delta y$, pix")

    for ax in [axx, axy, axxy, ax4, ax5]:
        ax.grid(which='both')
        leg=ax.legend()

    # save the figures
    figx.savefig('fig_shiftsVsTime.pdf')
    figxy.savefig('fig_shiftsVsShift.pdf')

# plt.show(block=False)

#print xyShift

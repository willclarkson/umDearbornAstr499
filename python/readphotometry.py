import sys
import os
import numpy as np
from astropy.io import fits
from astropy.stats import LombScargle
import matplotlib.pyplot as plt
plt.ion() # ensure we get the command line back

# the glob library contains a function, also called glob, that finds files and directories whose names match a pattern. We provide those patterns as strings: the character *
# matches zero or more characters, while ? matches any one character. We can use this to get the names of all files in the directory.
import glob

#from astropy.table import Table
#from astropy.io import fits

# for serializing to disk (will still be quicker to read than to re-generate)
import cPickle as pickle


# a little method that when called will print the names of all the important methods within this routine 
# and instructions for each, if needed.
def readphotometryHelp():
    print ''
    print 'read()'
    print 'The method read() reads in all the photometry files.'
    print 'This method will be called within read_and_plot()'
    print 'INSTRUCTIONS: Use read_and_plot()'
    print ''
    print 'read_and_plot()'
    print 'The method read_and_plot() will call the read() method, which searches for .fits.flux files for lightcurve generation.'
    print 'INSTRUCTIONS: Enter in the argument for read_and_plot() the prefix of the .fits.flux file names. (e.g. prefix="alignedCropped")'
    print ''
    print 'starPos()'
    print 'The method starPos() will read the FITS header information of a FITS image that has been aligned and cropped in AstroImageJ.'
    print 'INSTRUCTIONS: There are no specific arguments needed for this method to work.'

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
    #Lfound = glob.glob('%s*.fits' % (prefix))
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
        #pathImg = '../%s' % (f.split('.fits')[0]) # for difference imagges ? 
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
    
    plt.style.use('classic') # changing the plot style
    # CALLS READ() ROUTINE TO BRING IN THE DATA
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
        plt.scatter(times, fluxAll[:,iThis], alpha=0.5)
    ax.set_xlabel(sLabelX)
    ax.set_ylim(-5000, 5000)
    leg = ax.legend()

    # shows the lightcurve plot
    plt.show()
    


    return fluxAll, unctAll
 
# -------------------------------------------------------------------------- O ------------------------------------------------------------------------

# read-in star positions
def starPos(degFitX=2, sTitl=''):

    jdAll = np.array([])
    xShiftAll = np.array([])
    yShiftAll = np.array([])
    print('Reading in files...')
    
    files = glob.glob('aligned*.fits')
    #files = glob.glob('*.fits')            
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
    #print jdAll
#   update the time vector in-place to convert from days --> seconds
    #print jdAll
    jdAll = (jdAll - np.min(jdAll))*86400.

    # fit low-order polynomial to the shifts
    parsX = np.polyfit(jdAll, xShiftAll, deg=degFitX)
    parsY = np.polyfit(jdAll, yShiftAll, deg=degFitX)

    tFine = np.linspace(np.min(jdAll)-0.01, np.max(jdAll)+0.01, 1000, endpoint=True)

#----------------- FIGURE 1 AND 2 -------------------#
#Plot data
    figx = plt.figure(1)
    figx.clf()
    axx = figx.add_subplot(211)
    dum1 = axx.scatter(jdAll, xShiftAll, marker='^', s=9, c='b', \
                       label=r"xShift %s" % (sTitl))
    
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
                label=r"yShift %s" % (sTitl))
    #plt.show(block=False)
    dum3 = axy.plot(tFine, np.polyval(parsY, tFine), 'k-')

    #axy.set_xlabel('JD - min(JD), seconds')
    sLabelT =r"$\Delta t$, seconds"
    axy.set_xlabel(sLabelT)
    #axy.set_xlabel(r"$\Delta t$, seconds")

#    for ax in [axx, axy]:
#        ax.grid(which='both')
#        leg=ax.legend()

    axx.set_ylabel('xShift (pixels)')
    axy.set_ylabel('yShift (pixels)')
    plt.suptitle('X-Shift and Y-Shift vs. Time')

#---------------- FIGURE 3 -----------------------#

    figxy = plt.figure(3)
    figxy.clf()
    axxy = figxy.add_subplot(111)
    dum = axxy.scatter(xShiftAll, yShiftAll, c=jdAll, marker='o', label='path %s' % (sTitl), \
        edgecolor='0.5')
    axxy.set_xlabel('xShift (pixels)')
    axxy.set_ylabel('yShift (pixels)')
    plt.title('X-Shift vs. Y-Shift')
    cbar = figxy.colorbar(dum)
    
    # let's get the residuals from the straight line fit in X and the fit in Y
    residX = xShiftAll - np.polyval(parsX, jdAll)
    residY = yShiftAll - np.polyval(parsY, jdAll)

#--------------- FIGURE 4 -----------------------#

    fig4 = plt.figure(4)
    fig4.clf()
    ax4 = fig4.add_subplot(211)
    dum4 = ax4.plot(jdAll, residX, color='b', ms=3, label='X Residual %s' %(sTitl), marker='^')
#ax4.set_xlabel(r"JD (seconds)")
    ax4.set_ylabel(r"$\Delta x$, pix")
    ax5 = fig4.add_subplot(212, sharex=ax4)
    dum5 = ax5.plot(jdAll, residY, color='r', ms=3, label='Y Residual', marker='s')
    #ax5.set_xlabel(r"JD(seconds)")
    ax5.set_xlabel(sLabelT)
    ax5.set_ylabel(r"$\Delta y$, pix")
    plt.suptitle('X and Y Residuals vs. Time')

#-------------- FIGURE 5 ------------------------#

    # Creating Lomb-Scargle Periodograms
    # imported from astropy.stats import LombScargle
    fig5 = plt.figure(5)
    fig5.clf()
    ax6 = fig5.add_subplot(212)
    # fig.suptitle() (could also use plt.suptitle()) displays a tilte above both subplots
    plt.suptitle('Lomb-Scargle Power Spectrum for X and Y Residuals')


# let's generate the frequencies we want
    pDesired = np.linspace(500., 1500., 1000, endpoint=True)
    freq = 1.0/pDesired
    power = LombScargle(jdAll, residY).power(freq)

# let's get the maximium value
    iMax = np.argmax(power)
    periodMax = pDesired[iMax]

    sPeakY = 'Peak period = %.2fs' % (periodMax)

#freq, power = LombScargle(jdAll,residY).autopower()
    ax6.plot(1.0/freq, power, 'ro', ls='-', ms=2, label=sPeakY)
#ax6.set_xlabel('Period (s)')
    ax6.set_xlim(500., 1500.)
    # leg6 = ax6.legend()

# Plot for X residuals
# sPeakX = 'Peak period = %.2f s' % (periodMax)
    powerX = LombScargle(jdAll, residX).power(freq)
    iMaxX = np.argmax(powerX)
    periodMaxX = pDesired[iMaxX]
    
    sPeakX = 'Peak period = %.2fs %s' % (periodMaxX, sTitl)

    ax7=fig5.add_subplot(211)
    ax7.plot(1.0/freq, powerX, 'bo', ls='-', ms=2, label=sPeakX)
    ax6.set_xlabel('Period (s)')
#leg7 = ax7.legend()
    avgPeriod = np.average([periodMax,periodMaxX])
    phase = jdAll/np.float(avgPeriod)
    print np.shape(phase)
    phase = phase - np.floor(phase)

# let's fit the residuals against each other
    parsXY = np.polyfit(residX, residY, 1)
# ------------------ FIGURE 6 --------------- #
    fig6 = plt.figure(6)
    fig6.clf()
    axdxdy = fig6.add_subplot(111)

# let's try ordering by phase
    iSort = np.argsort(phase)

# increase the plot symbol size with phase
    sPhs = 25.0 + 50.*phase**2
#sPhs = np.repeat(25., np.size(phase))

    ax8 = axdxdy.scatter(residX[iSort], residY[iSort], c=phase[iSort], \
                        marker='o', label='path %s' % (sTitl), edgecolor='0.5', s=sPhs, \
                            vmin=0., vmax=1., zorder=2)


#dumLine = axdxdy.plot(residX[iSort], residY[iSort], ls='-.', lw=1, color='0.8')
    axdxdy.set_xlabel('residX (pixels)')
    axdxdy.set_ylabel('residY (pixels)')
    cbar1 = fig6.colorbar(ax8)
#axdxdy.set_xlim(500., 1500.)

#a = np.array([residX])
#    b = np.array([residY])

    c = np.hstack((residX, residY))
    cMax = np.max(np.abs(c))

    # let's overplot the trend HERE
    xFine = np.linspace(-cMax, cMax, 1000)
    dumTrend = axdxdy.plot(xFine, np.polyval(parsXY, xFine), color='k', ls='--', zorder=1, \
                           label=r'Trend angle: %.2f$^{\circ}$' % (np.degrees(np.arctan(parsXY[0]))) )


    axdxdy.set_xlim(-cMax,cMax)
    axdxdy.set_ylim(-cMax,cMax)
    plt.title('Residuals in X vs. Residuals in Y')
#    axdxdy.set_xlim(-np.max(np.abs(c)), np.max(np.abs(c)))
#    axdxdy.set_ylim(-np.max(np.abs(c)), np.max(np.abs(c)))
#axdxdy.set_xlim(


    for ax in [ax6, ax7]:
        ax.set_ylabel('Lomb-Scargle power')
    



    for ax in [axx, axy, axxy, ax4, ax5, ax6, ax7, axdxdy]:
        ax.grid(which='both')
        leg=ax.legend()

    figList = [figx, figxy, fig4, fig6, fig5 ]
    figTails = ['shiftVsTime', 'shiftVsShift', 'residVsTime', 'residVsResid','powspec']

    # clean up the string to make it filename-appropriate
    sFil=''
    if len(sTitl) > 0:
        sFil = sTitl[:].replace("(","").replace(")","")
        sFil = sFil.replace(",","_").replace(" ","")
        sFil = sFil.replace("/","").replace("\\","")
        sFil = "%s_" % (sFil)
        #print "SFIL INFO: %s" % (sFil)
        print '---------Date---file name--'
        print ''

    for iFig in range(len(figList)):
        fileName = 'fig_%s%s.pdf' % (sFil, figTails[iFig])
        figList[iFig].savefig(fileName)

        
        print fileName, 'SAVED' 
            # save the figures
#figx.savefig('fig_shiftsVsTime.pdf')
#    figxy.savefig('fig_shiftsVsShift.pdf')

# plt.show(block=False)

#print xyShift

# finds the largest reisdual values and scales the x and y limits of the plot accordingly
#def resid():




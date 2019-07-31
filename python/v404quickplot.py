import sys, glob
import numpy as np
import matplotlib.pylab as plt
plt.ion()
from astropy.table import Table
from astropy.io import fits

def go(sSrch='FluxData_2*.txt', Debug=False, showMag=True, binMinutes=5.3):

    # let's look for our photom files
    lPhot = glob.glob(sSrch)

    iMin = 0
    iMax = len(lPhot)
    if Debug:
        iMax = 1

    # initialize master arrays
    tAll = np.array([])
    rAll = np.array([])
    eAll = np.array([])
    nAll = np.array([])

    # let's stack the point-by-point photometry too
    tFine = np.array([])
    rFine = np.array([])
    eFine = np.array([])

    for iFile in range(iMin, iMax):
        # thisTime, thisRate, thisError, thisNbin = readAndBin(lPhot[iFile])

        tRaw, rRaw, eRaw = np.genfromtxt(lPhot[iFile], unpack=True)
        thisTime, thisRate, thisError, thisNbin = \
            BinData(tRaw, rRaw, vError=eRaw, Verbose=False, BinTime=binMinutes/1440.)

        # now we np.hstack each 1D array onto its corresponding 1D master array
        tAll = np.hstack((tAll, thisTime))
        rAll = np.hstack((rAll, thisRate))
        eAll = np.hstack((eAll, thisError))
        nAll = np.hstack((nAll, thisNbin))

        # let's stack our raw data too
        tFine = np.hstack(( tFine, tRaw ))
        rFine = np.hstack(( rFine, rRaw ))
        eFine = np.hstack(( eFine, eRaw ))

        print lPhot[iFile], np.shape(tAll), np.shape(rAll), np.shape(eAll), np.shape(nAll)

    # all being well, after exitting the loop we should have the
    # master tAll, rAll etc. for the entire set of files. THESE we can
    # plot.

    #  tBin, rBin, eBin = readAndBin(lcFile)

    # we'll create separate variables to plot so that we can control what they do
    rSho = np.copy(rAll)
    eSho = np.copy(eAll)

    rShoFine = np.copy(rFine)
    eShoFine = np.copy(eFine)

    sYaxis = 'Flux relative to reference star'

    if showMag:
        rSho = -2.5*np.log10(rAll)
        eSho = 1.086 * eAll

        rShoFine = -2.5*np.log10(rShoFine)
        eShoFine = 1.086 * eShoFine

        sYaxis = r'$\Delta mag$ relative to reference star'

    # syntax to plot would come here...
    fig = plt.figure(1)
    fig.clf()
    plt.scatter(tAll, rSho, color='b', zorder=10)
    plt.errorbar(tAll, rSho, yerr=eSho, color='b', alpha=0.5, ls='none', zorder=10)
    plt.gca().set_ylabel(sYaxis)
    plt.gca().set_xlabel('MJD')

    # we get the axis limits here so that we can apply them below
    yAx = plt.gca().get_ylim()

    # let's underplot the raw data
    faintColor = '0.9'
    plt.scatter(tFine, rShoFine, color=faintColor, zorder=5, s=4)
    plt.errorbar(tFine, rShoFine, yerr=eShoFine, color=faintColor, zorder=5, ls='none', ms=2)
    plt.gca().set_ylim(yAx)

    plt.show()
    

def readAndBin(lcFile='blah.txt', Verbose=False):

    # tThis = [genfromtxt table reads here from lcFile]
    time, rate, error = np.genfromtxt(lcFile, unpack=True)

    if Verbose:
        print "INFO:", np.shape(time)

    # binData [on the time, rate, error lifted from the resulting table]

    # (You'll need to know the column names to assign to time, rate,
    # error, but that might just be 'col1', 'col2', 'col3' if there
    # are none in the AIJ output. print tThis.colnames might help!

    # [bin the data]
    vOutTime, vOutRate, vOutError, vNPerBin = BinData(time, rate, vError=error, Verbose=False)

    # return the binned data to the calling routine
    return vOutTime, vOutRate, vOutError, vNPerBin

def BinData(vTime=np.array([]), vRate=np.array([]), vError=np.array([]), \
                nMin=1, tStart=-1e9, tEnd=-1.0e9, BinTime=0.003703703703703704, \
                Verbose=True, binIsSeconds=False):

    """Bin data given a time-series.

    RETURNS binned time, rate, uncertainty, and the number per bin"""

    if binIsSeconds:
    	binTime = binTime / 86400.

	# Initialize output arrays:
    vOutTime = np.array([])
    vOutRate = np.array([])
    vOutError = np.array([])
    vNPerBin = np.array([])

    if np.size(vTime) < 1:
    	return vOutTime, vOutRate, vOutError, vNPerBin

    # tStart is needed; tEnd is not needed yet!
    if tStart < -1e8:
        tStart = np.min(vTime)

    if tEnd < -1e8:
        tEnd = np.max(vTime)

    # NOTE - user might have put tend < tstart...
    if tEnd < tStart:
        print "BINDATA INFO: tEnd < tStart: I assume you want tStart + tEnd"
        tEnd = tStart + tEnd

    # report the settings if Verbosity is "on":
    if Verbose:
        print "BinData settings:"
        print "tStart: %.3f , tEnd: %.3f" % (tStart, tEnd)
        print "nMin = %i" % (nMin)
        print "BinTime = %.2e" % (BinTime)
        print "=================="

    # Now we do the binning:
    vOfBins = (vTime - tStart) / BinTime  # Which time bin does each
                                          # point fall into?
    lOfBins = np.asarray(vOfBins,'int') # Bin ID of each data point

    # Now we calculate the average rate and error in each time
    # bin. For each time bin, we determine (i) if any data points fall
    # in that bin, and (ii) what the average t, rate, error values
    # actually are. Because we calculated the bin id's per datapoint,
    # there should not be any empty bins at all, but we'll keep the
    # conditional in there for robustness. 

    # We can also impose a minimum number of points for a
    # "trustworthy" bin.
    
    # Let's only loop through the unique bins (don't want to calculate
    # once for every point within every bin!)
    vBinIDs = np.unique(lOfBins)
    
    # vector of bin start times
    vBinTimes = vBinIDs * BinTime + tStart

    # Initialize output arrays:
    vOutTime = np.array([])
    vOutRate = np.array([])
    vOutError = np.array([])
    vNPerBin = np.array([])

    if Verbose:
        print np.shape(vBinIDs), np.min(vTime), np.max(vTime), np.min(vBinTimes), np.max(vBinTimes), nMin

    for iBin in range(0, np.size(vBinIDs)):
        
        # When does this bin start and end?
        ThisTStart = vBinTimes[iBin]
        ThisTEnd   = ThisTStart + BinTime

        # ignore data > tend
        if ThisTEnd >= tEnd:
            if Verbose:
        	print "HALT CONDITION: %i" % (iBin)
            break

        # all datapoints inside this bin, we average together
        gInThisBin = np.where( (vTime >= ThisTStart) & (vTime < ThisTEnd) )[0]

        if Verbose:
            sys.stdout.write("\r DBG: this bin: %.2f, %.2f, %i" \
                                 % (ThisTStart, ThisTEnd, np.size(gInThisBin) ) )
            sys.stdout.flush()
        
        # If there are *no* datapoints in this bin, junk it and move on
        if np.size(gInThisBin) < nMin:
            continue  # ignore all following instructions within this
                      # loop and go to the next bin.

        # if there *are* more than nMin, we can proceed! Find the
        # average t, r, e values.
        ThisTimeAverage = np.mean(vTime[gInThisBin])
        ThisRateAverage = np.mean(vRate[gInThisBin])
        
        # Error combination is slightly more involved - add the quad
        # sum of the errors (remember propagation of errors from
        # classes)
        ThisErrorAverage = np.sqrt(np.sum(vError[gInThisBin]**2) / np.size(gInThisBin)**2 )

        # Having found the average in all the inputs, stick them onto
        # the end of the output
        vOutTime  = np.hstack(( vOutTime, ThisTimeAverage ))
        vOutRate  = np.hstack(( vOutRate, ThisRateAverage ))
        vOutError = np.hstack(( vOutError, ThisErrorAverage ))
        vNPerBin = np.hstack(( vNPerBin, np.size(gInThisBin) ))
            
    if Verbose:
        print "BinData INFO - number of non-empty bins: %i" % (np.size(vOutTime) )
        print "BinData INFO - number of input datapoints %i" % (np.size(vTime) )

    return vOutTime, vOutRate, vOutError, vNPerBin




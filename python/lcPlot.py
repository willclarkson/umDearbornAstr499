#
# lcPlot.py
# test for github 2018-05-09

# 2018-04-05 WIC - plot night-by-night lightcurve a la Z04

from astropy.table import Table
import numpy as np
import v404_test
import os

import matplotlib.pylab as plt
plt.ion()

#plt.style.use('seaborn-white')
plt.style.use('ggplot')

# a few utilities for doing the lomb-scargle
from astroML.time_series import lomb_scargle, lomb_scargle_bootstrap

def go(times=np.array([]), mags=np.array([]), unctys=np.array([]), \
           dayNo=np.array([]), \
           nCols=4, \
           test=True, nNights=8, logPer=True,
           nNoise=1000, pctile=5., \
           lsOnAll=True, degPoly=2, errScale=1.0, showFit=True, \
           filName='UnknownYear.txt', write=False, binLS=False, select=True, writeLS=False):
    
    """Plots (time, mag, error) data partitioned by day number. 

    if test=True, data are generated up to nNights.

    If degPoly < 0, no detrending is performed."""


    # if test is set, generate the data
    if test:
        times, mags, unctys = genData(nNights=nNights)
        dayNo = np.asarray(times, 'int')

    # now break into an n x 4 grid
    nPanels = int(np.max(dayNo)+1)  # can tweak this if we're counting from
                               # 1
    nRows = int(np.ceil(nPanels/np.float(nCols)))

    # The lightcurve figure
    fig1=plt.figure(1)
    fig1.clf()

    # The lomb-scargle figure
    figLS = plt.figure(3)
    figLS.clf()

    # now for the white-noise test
    figLSnoise = plt.figure(4)
    figLSnoise.clf()

    # remove spacing between plots - for all the figures
    for fig in [fig1, figLS, figLSnoise]:
        fig.subplots_adjust(hspace=0.02, wspace=0.05)
    
    # list of axes (so that we can get them back later)
    lAxes = []
    hrMin = 99.
    hrMax = -99.

    # we can generalize this later. For the moment, let's do separate
    # variables per figure
    perMin = 99.
    perMax = -99.
    lsMin =  999.
    lsMax = -999.

    if logPer:
        perMin = 1e-1
        perMax = 10.

    # axes list for the l-s and the gaussian noise l-s
    lAxesLS = []
    lAxesLSnoise = []

    coloSim='g'
    print "nPanels: %i, nRows:%i, nCols:%i" % (nPanels, nRows, nCols)

    sigTable = Table([[0.], [0.]], names=('JD', 'std'))

    ## 2019-04-11 - let's set up a master 'detrended' array for the LS
    detMJD= np.array([])
    detMag = np.array([])
    detUnc = np.array([])

    # 2019-04-18 WIC - hacked to handle less orthodox panel ordering
    rowNotEmpty = False

    # loop through the panels
    # hack for 17:
    #iPanels = [3,5]
    #for iPlot in iPanels:
    for iPlot in range(nPanels):
        bThis = dayNo == iPlot


        # 2019-04-18 WIC - hack for labeline
        if iPlot % nCols < 1:
            rowNotEmpty = False

        if np.sum(bThis) < 10:
            print "Bad interval: %i, %i points" % (iPlot, np.sum(bThis))
            continue
        
        # add the panel
        ax=fig1.add_subplot(nRows, nCols, iPlot+1)
        
        # times from the beginning of the dataset. NOTE - when doing
        # this with real data, we can use the actual decimal day
        # number. With simulated data, though, we use the first time
        # in each interval.
        hrs = 24.* (times[bThis] - np.min(times[bThis]))

        # do the plot
        dum = ax.errorbar(hrs, mags[bThis], unctys[bThis], \
                              ls='None', ecolor='b', \
                              marker='o', color='b', \
                              ms=2, alpha=0.5, elinewidth=1)



        # Set the label for the day number
        sLabel = 'Night %i' % (iPlot + 1)

        # Label the day number
        dumAnno = ax.annotate(sLabel, (0.95, 0.95), \
                                  xycoords='axes fraction', \
                                  ha='right', va='top', \
                                  fontsize=12, color='k')

        # if we're more than nCols from the end, hide the horizontal
        # axis
        if nPanels - iPlot > nCols and iPlot > 0:
            ax.tick_params(labelbottom='off')
        else:
            ax.set_xlabel('Time (hours)')

            # 2019-04-18 WIC - another hack to suppress the last time label to avoid overlap
            
            #xticks = ax.xaxis.get_major_ticks()
            #xticks[-1].set_visible(False) 

        #2019-04-19 WIC - that hack again
        if not rowNotEmpty:
            ax.set_ylabel(r'$\Delta R$ (mag)')
            rowNotEmpty = True

            # 2019-04-18 WIC - 
            # if this is anything other than the upper-left corner, hide the uppermost tick since it'll overlap
            if iPlot > 0:
                yticks = ax.yaxis.get_major_ticks()
                yticks[0].set_visible(False) 

        else:

#        if iPlot % nCols > 0:
            ax.tick_params(labelleft='off')
#        else:
#            ax.set_ylabel(r'$\Delta R$')
        
        

        # record the maximum axis length so that we can access it
        # later
        hrMax = np.max([hrMax, np.max(hrs) ])
        hrMin = np.min([hrMin, np.min(hrs) ])

        lAxes.append(ax)

        # 2019-04-02: Find the standard deviation and output to the terminal

        sigz, pars = findResidualSigma(hrs, mags[bThis], unctys[bThis], errScale=errScale, degPoly=degPoly)

        # let's re-produce the detrended points and stack them back on to the master array
        thisDetMag = mags[bThis] - np.polyval(pars, hrs)
        detMJD = np.hstack(( detMJD, times[bThis]))
        detMag = np.hstack(( detMag, thisDetMag))
        detUnc = np.hstack(( detUnc, unctys[bThis]))

        # Show the poly-fit in Figure 1

        if showFit:
            fit = ax.plot(hrs, np.polyval(pars, hrs), 'k--', lw=1)

        print "Night %i: stddev %.5f" % (iPlot+1, sigz)

        # Create a table with the format (JD, sigma)

        meanTime = np.mean(times[bThis])

        sigTable.add_row([meanTime, sigz])

        # now - IF there are enough datapoints - we do the
        # Lomb-Scargle
        if np.size(hrs) < 40:
            continue
        
        dtMin = np.min(hrs[1::] - hrs[0:-1])
        dtRange = np.max(hrs) - np.min(hrs)
        
        # Generate periods array out of this and do the L-S
        #pers = np.logspace(np.log10(dtMin*2.), np.log10(dtRange*1.0),1000) 
        pers = np.logspace(np.log10(dtMin*2.), np.log10(7.), 200)
        omegas = 2.0 * np.pi / pers
        
        lsBin = lomb_scargle(hrs, mags[bThis], unctys[bThis], \
                                 omegas, \
                                 generalized=False)

        # Do the same exercise for gaussian white-noise
        # print "Starting %i trials..." % (nNoise)
        lsNoiseUpper, lsNoise = whiteNoiseLS(hrs, unctys[bThis], \
                                                 omegas, nNoise, \
                                                 pctile)
        
        #aNoiseLS = np.zeros((nNoise, np.size(omegas)))
        #for iNoise in range(nNoise):
        #    magNoise = np.random.normal(size=np.size(hrs))*unctys[bThis]
        #    lsNoise = lomb_scargle(hrs, magNoise, unctys[bThis], \
        #                               omegas, \
        #                               generalized=False)

        #    # slot this in to the simulation array
        #    aNoiseLS[iNoise] = lsNoise
            
        ## do upper and lower bounds from the noise
        #lsNoiseUpper = np.percentile(aNoiseLS, 100.-pctile, axis=0)
        #lsNoiseLower = np.percentile(aNoiseLS, pctile, axis=0)

        # OK now we plot the LS for this particular axis. Just like
        # before, this time for a different figure.
        axLS = figLS.add_subplot(nRows, nCols, iPlot+1)
        axLSnoise = figLSnoise.add_subplot(nRows, nCols, iPlot+1)
        

        ### 2018-05-07 log-log from semilogx
        dumLS = axLS.loglog(pers, lsBin, 'bo', ms=1, ls='-', color='b')
        dumLSnoise = axLSnoise.loglog(pers, lsNoise, \
                                            ms=1, ls='-', color='0.2', \
                                            alpha=0.5)
        
        coloSim = 'darkmagenta'
        dum1 = axLSnoise.plot(pers, lsNoiseUpper, alpha=0.7, \
                                  color=coloSim)
        #dum2 = axLSnoise.plot(pers, lsNoiseLower)

        dumAnno = axLSnoise.annotate('%.0f percent, %i trials' \
                                         % (100.-pctile, nNoise), \
                                         (0.05, 0.85), \
                                         xycoords='axes fraction', \
                                         ha='left', va='top', \
                                         color=coloSim, \
                                         alpha=0.7)

        # LS powers for standardization
        perMin = np.min([perMin, np.min(pers)])
        perMax = np.max([perMax, np.max(pers)])
        lsMin = np.min([lsMin, np.min(lsBin)])
        lsMax = np.max([lsMax, np.max(lsBin)])

        # do the axes as before
        if nPanels - iPlot > nCols:
            axLS.tick_params(labelbottom='off')
            axLSnoise.tick_params(labelbottom='off')
        else:
            axLS.set_xlabel('Period (hours)')
            axLSnoise.set_xlabel('Period (hours)')

        if iPlot % nCols > 0:
            axLS.tick_params(labelleft='off')
            axLSnoise.tick_params(labelleft='off')
        else:
            axLS.set_ylabel('LS Power')
            axLSnoise.set_ylabel('LS Power')

        # annotate the day number on both LS axes
        for axAnno in [axLS, axLSnoise]:
            dumdum = axAnno.annotate(sLabel, (0.05, 0.95), \
                                         xycoords='axes fraction', \
                                         ha='left', va='top', \
                                         fontsize=12, color='k')

        # append the LS axis onto the list
        lAxesLS.append(axLS)
        lAxesLSnoise.append(axLSnoise)

        ### 2018-05-07 hardcode the vertical axis
        ### ax.set_ylim(np.max(mags[bThis]), np.max(mags[bThis])-0.4)

    # now set the same x-range for all the axes
    for thisAx in lAxes:
        thisAx.set_xlim(hrMin, hrMax)

        # set the y limit
        ### 2018-05-07 
        thisAx.set_ylim(np.max(mags + unctys), np.min(mags - unctys))

        # 2019-04-09: Set the y-range for all the axes to be the same as Z04 Figure 3

        thisAx.set_ylim([0.05, -0.3])
        
    # do the same for the lomb-scargle figures
    for iAx in range(len(lAxesLS)):
        thisLS = lAxesLS[iAx]
        try:
            thisLS.set_xlim(perMin, perMax)
            thisLS.set_ylim(lsMin, lsMax)
        except:
            print("WARN - LS axes badval??")
        #thisLS.set_ylim(lsMin, lsMax)
        thisLS.grid(which='both', visible=True)


        thisLSnoise = lAxesLSnoise[iAx]
        try:
            thisLSnoise.set_xlim(perMin, perMax)
            thisLSnoise.set_ylim(lsMin, lsMax)
        except:
            print("WARN - LS Noise axes badval??")
        thisLSnoise.grid(which='both', visible=True)
        
    # output the plot to a figure
    fig1.savefig('test_grid_lc.png', transparent=True)
    figLS.savefig('test_grid_LS.png', transparent=True)
    figLSnoise.savefig('test_grid_LSnoise.png', transparent=True)

    fig2=plt.figure(2)
    fig2.clf()
    ax2 = fig2.add_subplot(111)
    dum = ax2.scatter(times, mags, c=dayNo)

    cbar = fig2.colorbar(dum, ax=ax2)

    if not lsOnAll:
        return

    print "lcPlot.go INFO - starting on the full dataset..."

    # do the lomb-scargle on the entire dataset
    dtAll = np.max(times) - np.min(times)
    dtAll = 100./24.  # (100 hours)
    dtAll = 12.0/24.0 # (12 hours max)
    dtMin = np.min(times[1::] - times[0:-1])
    perAll = np.logspace(np.log10(dtMin*2.), np.log10(dtAll), 250)
    omeAll = 2.0 * np.pi / perAll

    # let's do the LS on the detrended data

    # lsAll = lomb_scargle(times, mags, unctys, omeAll, generalized=False)

    lsAll = lomb_scargle(detMJD, detMag, detUnc, omeAll, generalized=False)
    lsUp, lsNoise = whiteNoiseLS(times, unctys, omeAll, nNoise, pctile)

    fig5 = plt.figure(5, figsize=(8,4))
    fig5.clf()

    # try replacing perAll * 24. with omeAll / (2.0*np.pi)
    fAll = omeAll / (2.0*np.pi * 86400.)

    # take log10 of frequency and power

    logFreq = np.log10(fAll)
    logPower = np.log10(lsAll)

    # Save this to a table

    tbl = Table([fAll, lsAll])
    if writeLS:
        tbl.write(filName, format='ascii')
    # Try to select only the data that have a frequency between -3.9 and -2.8

    pwrFit = np.polyfit(logFreq, logPower, 1)

    ax51 = fig5.add_subplot(121)
    ax52 = fig5.add_subplot(122, sharex=ax51, sharey=ax51)

#    LS = ax51.loglog(perAll * 24., lsAll, 'bo', ms=1, ls='-', color='b')
#    LSc = ax52.loglog(perAll * 24., lsNoise, \
#                            ms=1, ls='-', color='0.2', \
#                                            alpha=0.5)
#    LSn = ax52.loglog(perAll * 24., lsUp, alpha=0.7, color=coloSim)
    LS = ax51.loglog(fAll, lsAll, 'bo', ms=1, ls='-', color='b')
    LSc = ax52.loglog(fAll, lsNoise, \
                            ms=1, ls='-', color='0.2', \
                                            alpha=0.5)
    LSn = ax52.loglog(fAll, lsUp, alpha=0.7, color=coloSim)


    fig5.subplots_adjust(hspace=0.02, wspace=0.05, bottom=0.15)

    # standardize
    for thisAx in [ax51, ax52]:
        thisAx.grid(which='both', visible=True)
#        thisAx.set_xlabel('Period (hours)')
        thisAx.set_xlabel(r'Frequency (s$^{-1}$)')
        
    ax51.set_ylabel('LS Power')
    ax52.tick_params(labelleft='off')

    # annotate the noise axis
    dumAnno = ax52.annotate('%.0f percent, %i trials' \
                                         % (100.-pctile, nNoise), \
                                         (0.05, 0.85), \
                                         xycoords='axes fraction', \
                                         ha='left', va='top', \
                                         color=coloSim, \
                                         alpha=0.7)

    # save this figure to disk
    fig5.savefig('test_alldata_LS.png', transparent=True)


    # Remove the first row from the table

    sigTable.remove_row(0)

    # See what the table looks like

    print sigTable

    # Write the table to disk

    if write:
        sigTable.write(filName, format='ascii')

    # Create a "Figure 7" to plot the logs of the LS in linear space (really logspace, since we're plotting logs)

    f7 = plt.figure(7)
    f7.clf()
    ax7 = f7.add_subplot(111)
    ax7.scatter(logFreq, logPower, color='b')
    ax7.plot(logFreq, logPower, color='b')
    ax7.plot(logFreq, np.polyval(pwrFit, logFreq), c='r', label=r"$\beta$" +": %5.2f" % pwrFit[0])
    ax7.legend()
    ax7.set_xlabel("log(Frequency)")
    ax7.set_ylabel("log(LS Power)")

    if binLS:


        nBins = 20

        freqBin, pwrBin, errBin, nPerBin = v404_test.BinData(fAll, lsAll, np.ones(len(fAll)), \
            tStart=10**-3.9, tEnd=10**-2.8, BinTime=(0.00146 / nBins))

        print np.shape(pwrBin)

        logFreqBin = np.log10(freqBin)
        logPwrBin = np.log10(pwrBin)
        pwrFitBin, covarBin = np.polyfit(logFreqBin, logPwrBin, 1, cov=True)

        print "Line, covariance matrix:", pwrFitBin, covarBin
        stdResid = np.std(logPwrBin - np.polyval(pwrFitBin, logFreqBin))
        print "Stddev of residuals: %.2e" % (stdResid)

        f8 = plt.figure(8)
        f8.clf()
        ax8 = f8.add_subplot(111)
        ax8.scatter(logFreqBin, logPwrBin, color='b')
        ax8.plot(logFreqBin, np.polyval(pwrFitBin, logFreqBin), c='r', label=r"$\beta$" +": %5.2f" % pwrFitBin[0])
        ax8.legend()
        ax8.set_xlabel("log(Frequency)")
        ax8.set_ylabel("log(LS Power)")

        # Monte-Carlo for beta values varying bin numbers

        #if MC:
        #mcTbl = Table()

def genData(nPoints=1000, nNights=6):

    """Utility to generate fake datapoints"""

    # let's assume a 7-hour night
    phs = np.random.uniform(size=nPoints) % 0.3
    lNight = np.random.random_integers(0,nNights-1, size=nPoints)

    times = np.asarray(lNight, 'float') + phs 

    # just generate gaussian random noise
    unctys = np.random.uniform(size=nPoints)*0.4
    mags = np.random.normal(size=np.size(unctys))*unctys

    # sort by time
    lSor = np.argsort(times)

    return times[lSor], unctys[lSor], mags[lSor]

def whiteNoiseLS(t=np.array([]), u=np.array([]), \
                     omegas=np.array([]), nTrials=100, \
                     pctile=1., Verbose=True):

    """Wrapper to do the lomb-scargle on white noise datasets"""

    if Verbose:
        print "whiteNoiseLS INFO - Starting %i trials..." % (nTrials)
    
    nData = np.size(t)
    aNoiseLS = np.zeros((nTrials, np.size(omegas)))
    for iNoise in range(nTrials):
        magNoise = np.random.normal(size=nData)*u
        lsNoise = lomb_scargle(t,magNoise, u, \
                                   omegas, \
                                   generalized=False)

        aNoiseLS[iNoise] = lsNoise

    lsNoiseUpper = np.percentile(aNoiseLS, 100.-pctile, axis=0)

    # return the lomb-scargle of the upper percentile and also an
    # example LS periodogram.
    return lsNoiseUpper, lsNoise

def showBinnedLC(filTable='v404_binSub.fits', nCols=3, \
                     nTrials=6, pctile=5., stopAfterChunk=False, \
                     degPoly=0, errScale=1.0, filName='UnknownYear.txt', write=False, \
                     binLS=False, select=True, writeLS=False, forProposal=False):

    """Loads photometry file and plots in our nice grid"""

    if not os.access(filTable, os.R_OK):
        print "showBinned WARN - cannot read input path %s" \
            % (filTable)
        return

        # 2019-04-use style sheet depending on what we're doing
    if not forProposal:
        plt.style.use('ggplot')
    else:
        plt.style.use('classic')
        #plt.style.use('seaborn-poster')
        #plt.style.use('seaborn-white')

    tPho = Table.read(filTable)

    # 2019-03-25: let's try cleaning the entire table by non-isolated times:
    bNonIsol = findNonIsolatedMJD(tPho['tBin'])

    # now we apply the isolation index to the table as we read it in:
    tPho = tPho[bNonIsol]

    times = tPho['tBin']
    mags = tPho['fBinSub']
    unctys = tPho['uBin'] # / 3.0
 
    print "showBinned INFO: minmax times: %.2f, %.2f" % (np.min(times), np.max(times))

    # sanity check - let's view the times then quit

    # 2019-03-25 - re-chunk into days (if this works, might put back into v404Test)
    # dayno = assignDays(times)

    # apply our brand new chunker here:
    dayno2 = chunkByMJD(times)

    if not stopAfterChunk:
        go(times, mags, unctys, dayno2, test=False, nCols=nCols, \
           nNoise=nTrials, pctile=pctile, errScale=errScale, write=write, filName=filName, \
           degPoly=degPoly, binLS=binLS, select=select, writeLS=writeLS)
        return

    tDiff = times-np.min(times)

    #if stopAfterChunk:
    #    bHi = (tDiff > 0.20)
    #    tDiff=tDiff[bHi]
    #    mags = mags[bHi]
    #    dayno = dayno[bHi]

    #    tDiff = times-np.min(times)
    iDiff = np.asarray(tDiff, 'int')

    dayno = np.floor(times - np.min(times))+0.
    #print np.min(dayno), np.max(dayno)

    # CHUNK METHOD 1: loop through the nights
    for iChunk in range(0, np.max(iDiff)+1):
        bThisChunk = (tDiff > iChunk - 0.4) & (tDiff < iChunk + 0.4)
        dayno[bThisChunk] = iChunk

    # CHUNK METHOD 2: go by the largest gap
    dayno2 = np.zeros(np.size(tDiff))
    deltat = tDiff - np.roll(tDiff,1)
    deltat[0] = np.max(deltat)

    # CHUNK 2B: let's remove points that are isolated in BOTH directions
    deltaLo = tDiff - np.roll(tDiff, -1)
    deltaLo[-1] = 0.
    isol = 10. * np.abs(deltat) * np.abs(deltaLo)
    bNonIsolated = isol < 0.5

    #times = times[bNonIsolated]
    #tDiff = tDiff[bNonIsolated]
    #mags = mags[bNonIsolated]
    #deltat = deltat[bNonIsolated]
    #deltaLo = deltaLo[bNonIsolated]
    #isol = isol[bNonIsolated]
    #dayno = dayno[bNonIsolated]
    #dayno2 = dayno2[bNonIsolated]
    #unctys = unctys[bNonIsolated]

    # now we have to re-compute deltat!!
    #dayno2 = np.zeros(np.size(tDiff))

    # retained for backwards compatibility with the plot
    deltat = tDiff - np.roll(tDiff,1)
    deltat[0] = np.max(deltat)

    ## now we pick out the times for which deltat > 0.5 days (or whatever our threshold is)
    #gNewDayLo = np.where(deltat > 0.5)[0]
    #for iChunk in range(np.size(gNewDayLo)-1):
    #    rowLo = gNewDayLo[iChunk]
    #    rowHi = gNewDayLo[iChunk+1]
    #    dayno2[rowLo:rowHi] = int(np.round(tDiff[rowLo]))

    ## now we fill in that last chunk
    #dayno2[gNewDayLo[-1]::] = int(np.round(tDiff[gNewDayLo[-1]]))

    # dayno2 = chunkByMJD(times)

    figX = plt.figure(6)
    figX.clf()
    axX = figX.add_subplot(111)
    dum = axX.scatter(tDiff, mags, alpha=0.5, s=1, c=dayno, cmap=plt.cm.jet)
    dum2 = axX.plot(tDiff, tDiff)
    dum2 = axX.plot(tDiff, tDiff, 'bo', ms=4, alpha=0.75)

    dum3 = axX.plot(tDiff, deltat, 'r-')
    dum3b = axX.plot(tDiff, deltat, 'ro')
    dum3c = axX.plot(tDiff, dayno2, 'r-.',lw=2)

    dum4 = axX.plot(tDiff, deltaLo, 'g-')
    dum4b = axX.plot(tDiff, deltaLo, 'go')
    dum4c = axX.plot(tDiff, isol, 'b-')
    dum4d = axX.plot(tDiff, isol, 'bo')

def chunkByMJD(times=np.array([]), dtMin=0.5):

    """Method to break a set of times into chunks"""

    # initialize the return
    daysRet = np.asarray(times * 0., 'int')

    tDiff = times - np.min(times)
    deltat = tDiff - np.roll(tDiff,1)
    deltat[0] = np.max(deltat)

    gNewDayLo = np.where(deltat > 0.5)[0]
    for iChunk in range(np.size(gNewDayLo)-1):
        rowLo = gNewDayLo[iChunk]
        rowHi = gNewDayLo[iChunk+1]
        daysRet[rowLo:rowHi] = int(np.round(tDiff[rowLo]))

    # now we fill in that last chunk
    daysRet[gNewDayLo[-1]::] = int(np.round(tDiff[gNewDayLo[-1]]))


    return daysRet

def findNonIsolatedMJD(times=np.array([]), isolMax=0.5, isolScale=10.):

    """Returns a boolean index for all non-isolated times"""

    tDiff = times - np.min(times)
    deltat = tDiff - np.roll(tDiff,1)
    deltat[0] = np.max(deltat)    

    deltaLo = tDiff - np.roll(tDiff, -1)
    deltaLo[-1] = 0.
    isol = isolScale * np.abs(deltat) * np.abs(deltaLo)
    bNonIsolated = isol < isolMax

    return bNonIsolated

def findResidualSigma(times=np.array([]), mags=np.array([]), unctys=np.array([]), errScale=1.0, degPoly=2):

    """Takes the residual of ellipsoidally-subtracted flare data and finds the standard deviation in quadrature.

    Set degPoly < 0 to switch off detrending."""

    # Consistency check

    print "errScale: ", errScale

    unctys = unctys * errScale

    # Fit a polynomial to the nightly light curve

    # 2019-04-08 - if degPoly < 0, don't do the fitting (but retain a zero-array for parameters 
    # so that nothing downstream breaks).
    if degPoly > -1:
        pars = np.polyfit(times, mags, deg=degPoly)
    else:
        pars = np.array([0.])


    # Subtract the polynomial from the data to find the residuals

    sub = mags - np.polyval(pars, times)

    # Find the standard deviation sigma_z

    unctMean = np.mean(unctys)
    if np.std(sub) > unctMean:
        stddev = np.sqrt(np.std(sub)**2 - unctMean**2)
    else:
        print "WARN- unctMean > std(sub). Sigma has been given the dummy value of -1."
        stddev = -1

    #print "std(sub): ", np.std(sub)
    #print "std(mags): ", np.std(mags)

    return stddev, pars

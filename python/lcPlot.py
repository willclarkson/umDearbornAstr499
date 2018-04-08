#
# lcPlot.py
#

# 2018-04-05 WIC - plot night-by-night lightcurve a la Z04

from astropy.table import Table
import numpy as np

import os

import matplotlib.pylab as plt
plt.ion()

plt.style.use('seaborn-white')

# a few utilities for doing the lomb-scargle
from astroML.time_series import lomb_scargle, lomb_scargle_bootstrap

def go(times=np.array([]), mags=np.array([]), unctys=np.array([]), \
           dayNo=np.array([]), \
           nCols=4, \
           test=True, nNights=8, logPer=True, \
           nNoise=1000, pctile=5., \
           lsOnAll=True):
    
    """Plots (time, mag, error) data partitioned by day number. 

    if test=True, data are generated up to nNights"""


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

    # loop through the panels
    for iPlot in range(nPanels):
        bThis = dayNo == iPlot
        if np.sum(bThis) < 1:
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
                              ls='None', ecolor='0.6', \
                              marker='o', color='b', \
                              ms=4, alpha=0.5)


        # Set the label for the day number
        sLabel = 'Night %i' % (iPlot + 1)

        # Label the day number
        dumAnno = ax.annotate(sLabel, (0.95, 0.95), \
                                  xycoords='axes fraction', \
                                  ha='right', va='top', \
                                  fontsize=12, color='k')

        # if we're more than nCols from the end, hide the horizontal
        # axis
        if nPanels - iPlot > nCols:
            ax.tick_params(labelbottom='off')
        else:
            ax.set_xlabel('Time (hours)')

        if iPlot % nCols > 0:
            ax.tick_params(labelleft='off')
        else:
            ax.set_ylabel(r'$\Delta R$')
        
        # record the maximum axis length so that we can access it
        # later
        hrMax = np.max([hrMax, np.max(hrs) ])
        hrMin = np.min([hrMin, np.min(hrs) ])

        lAxes.append(ax)

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
        

        dumLS = axLS.semilogx(pers, lsBin, 'bo', ms=1, ls='-', color='b')
        dumLSnoise = axLSnoise.semilogx(pers, lsNoise, \
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

    # now set the same x-range for all the axes
    for thisAx in lAxes:
        thisAx.set_xlim(hrMin, hrMax)

        # set the y limit
        thisAx.set_ylim(np.max(mags + unctys), np.min(mags - unctys))

    # do the same for the lomb-scargle figures
    for iAx in range(len(lAxesLS)):
        thisLS = lAxesLS[iAx]
        thisLS.set_xlim(perMin, perMax)
        thisLS.set_ylim(lsMin, lsMax)
        thisLS.grid(which='both', visible=True)

        thisLSnoise = lAxesLSnoise[iAx]
        thisLSnoise.set_xlim(perMin, perMax)
        thisLSnoise.set_ylim(lsMin, lsMax)
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
    dtMin = np.min(times[1::] - times[0:-1])
    perAll = np.logspace(np.log10(dtMin*2.), np.log10(dtAll), 500)
    omeAll = 2.0 * np.pi / perAll

    lsAll = lomb_scargle(times, mags, unctys, omeAll, generalized=False)
    lsUp, lsNoise = whiteNoiseLS(times, unctys, omeAll, nNoise, pctile)

    fig5 = plt.figure(5, figsize=(8,4))
    fig5.clf()
    ax51 = fig5.add_subplot(121)
    ax52 = fig5.add_subplot(122, sharex=ax51, sharey=ax51)

    LS = ax51.semilogx(perAll * 24., lsAll, 'bo', ms=1, ls='-', color='b')
    LSc = ax52.semilogx(perAll * 24., lsNoise, \
                            ms=1, ls='-', color='0.2', \
                                            alpha=0.5)
    LSn = ax52.semilogx(perAll * 24., lsUp, alpha=0.7, color=coloSim)

    fig5.subplots_adjust(hspace=0.02, wspace=0.05, bottom=0.15)

    # standardize
    for thisAx in [ax51, ax52]:
        thisAx.grid(which='both', visible=True)
        thisAx.set_xlabel('Period (hours)')
        
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
                     nTrials=1000, pctile=5.):

    """Loads photometry file and plots in our nice grid"""

    if not os.access(filTable, os.R_OK):
        print "showBinned WARN - cannot read input path %s" \
            % (filTable)
        return

    tPho = Table.read(filTable)
    times = tPho['tBin']
    mags = tPho['fBinSub']
    unctys = tPho['uBin']
    dayno = np.floor(times - np.min(times)-0.2)+1
    print np.min(dayno), np.max(dayno)

    go(times, mags, unctys, dayno, test=False, nCols=nCols, \
           nNoise=nTrials, pctile=pctile)

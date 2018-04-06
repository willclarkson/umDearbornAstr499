#
# lcPlot.py
#

# 2018-04-05 WIC - plot night-by-night lightcurve a la Z04

from astropy.table import Table
import numpy as np

import matplotlib.pylab as plt
plt.ion()

def go(times=np.array([]), mags=np.array([]), unctys=np.array([]), \
           dayNo=np.array([]), \
           nCols=4, \
           test=True, nNights=8):
    
    """Plots (time, mag, error) data partitioned by day number. 

    if test=True, data are generated up to nNights"""


    # if test is set, generate the data
    if test:
        times, mags, unctys = genData(nNights=nNights)
        dayNo = np.asarray(times, 'int')

    # now break into an n x 4 grid
    nPanels = np.max(dayNo)+1  # can tweak this if we're counting from
                               # 1
    nRows = int(np.ceil(nPanels/np.float(nCols)))

    fig1=plt.figure(1)
    fig1.clf()

    # remove spacing between plots
    fig1.subplots_adjust(hspace=0.02, wspace=0.05)

    # list of axes (so that we can get them back later
    lAxes = []
    hrMin = 99.
    hrMax = -99.

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
                              ls='None', ecolor='0.8', \
                              marker='o', color='0.1', \
                              ms=2, alpha=0.5)

        # label the day number (we can refine this to dates later)
        sLabel = 'Day %i' % (iPlot)
        dumAnno = ax.annotate(sLabel, (0.95, 0.92), \
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

    # now set the same x-range for all the axes
    for thisAx in lAxes:
        thisAx.set_xlim(hrMin, hrMax)

        # set the y limit
        thisAx.set_ylim(np.max(mags + unctys), np.min(mags - unctys))

    # output the plot to a figure
    fig1.savefig('test_lcGrid.png')

    fig2=plt.figure(2)
    fig2.clf()
    ax2 = fig2.add_subplot(111)
    dum = ax2.scatter(times, mags, c=dayNo)

    cbar = fig2.colorbar(dum, ax=ax2)

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

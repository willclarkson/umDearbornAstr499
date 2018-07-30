#
# loadOld.py
#

# 2018-05-07

from astropy.table import Table, Column
import numpy as np

import matplotlib.pylab as plt

def loadPhot(fileIn='199208elip_R.dat', fix98=True, fudge92=True):

    """Quick utility to read in the zurita photometry, using the ad hoc
rules for reading in the casares collection of zurita04 data

    fix98 = use fudge to convert 1998 data (for which we only have
    phase information) into MJD data

    fudge92 = assign uncertainties read off Z04 figure 3a?

    """

    # what might the comments strings be...
    lCommens = ['%', '!', '#']

    aData = np.array([])
    for sCommen in lCommens:
        try:
            aDum = np.genfromtxt(fileIn, comments=sCommen, unpack=False)
        except:
            badRead = True

    # Now we have to decide what this all means...
    tPho = Table()

    # are we dealing with times or phases here?
    nCols = np.shape(aDum)[-1]
    timeOrPhase = aDum[:,0]
    magn = aDum[:,1]

    sCol = 'time'
    if np.max(timeOrPhase) < 1.0:
        sCol = 'phase'

    tPho[sCol] = Column(timeOrPhase)

    # Flux or magnitude?
    # 2018-07-30 UPDATE - 1992 already IS in magnitude, just relative!!
    #if np.median(magn) > 15.0:
    #    tPho['mag'] = Column(magn)
    #else:
    #    tPho['flux'] = Column(magn)

    tPho['mag'] = Column(magn)

    # now for the other columns
    if nCols > 2:
        # 2018-07-30 OBSOLETE - 1992 is also in mag!!!
        #if np.median(magn) > 15.0:
        tPho['magErr'] = Column(aDum[:,2])
        #else:
        #    tPho['fluxErr'] = Column(aDum[:,2])

    if nCols > 3:
        tPho['magC'] = Column(aDum[:,3])
        tPho['magCerr'] = Column(aDum[:,4])

    # the 1992 data do not have uncertainties, but Z04 do provide an
    # estimate in their figure [we might consider going back to the
    # Pavlenko source paper to estimate the uncertainty]
    if len(tPho.colnames) < 3 and fudge92:
        tPho = estUncty92(tPho)
        
    # if we're fixing 1998 data, change things here, AFTER the table
    # has been completed (since we may be changing in-place).
    if np.min(timeOrPhase) < 1. and \
            np.max(timeOrPhase) - np.min(timeOrPhase) < 1.1 and \
            fix98:

        tPho = phaseToMJD(tPho)

    # attach metadata
    tPho.meta['file'] = fileIn[:]
        
    return tPho

def phaseToMJD(tPho=Table(), MJD0=51000.750, \
                   per=6.4714, tZer=48813.873, \
                   trimWeird=True, \
                   plotDBG=False):

    """Unpacks phase to MJD, with various assumptions about the
    observation date and row order. MJD0 is the MJD at the start of
    the run (will be used to find the nearest whole-number number of
    cycles for the ephemeris). Defaults to 1998 July 6th at 6pm UT.

    trimWeird: Two of the Zurita points appear to be
    mis-ordered. Remove them if trimWeird is set to "True."

    """
    if not 'phase' in tPho.colnames:
        return tPho # return unchanged

    phase = tPho['phase']

    # WATCHOUT - The Casares et al. ephemeris expresses tZer as JD - 2
    # 400 000.0 which means MJD + 0.5. For ease of reading, we add
    # that 0.5 days back on here.
    tZ = tZer + 0.5
    
    N_start = np.floor((MJD0 - tZ)/per)

    # which means the next lowest phase zero must occur at MJD...
    mjdPrevZero = tZ + N_start*per

    # it will be handy to have an array giving "weird" points for
    # which the time appears to have moved backwards...
    bWeird = np.repeat(False, np.size(phase))

    # we loop through these since it's the row number that preserves
    # the phase ordering
    mjdCalc = np.zeros(np.size(phase))
    iOrbit = 0.0
    for iRow in range(np.size(phase)):

        # which orbit are we on now?
        if iRow > 0:
            if phase[iRow-1] - phase[iRow] > 0.5:
                iOrbit = iOrbit + 1.0

        # now compute the mjd calc
        nOrbs = phase[iRow]+iOrbit
        mjdCalc[iRow] = mjdPrevZero + nOrbs * per
    
        if iRow > 0:
            if mjdCalc[iRow] < mjdCalc[iRow-1]:
                bWeird[iRow] = True

    # update in-place
    tPho['time'] = Column(mjdCalc)


    if plotDBG:
        plt.figure(1)
        plt.clf()
        lCount = np.arange(np.size(phase))
        plt.plot(lCount, phase, 'bo', ls='-')
        plt.plot(lCount[bWeird], phase[bWeird], 'rx', zorder=25)
        plt.xlabel('Row number')
        plt.ylabel('Phase')

        plt.figure(2)
        plt.clf()
        plt.scatter(tPho['time'], tPho['mag'], c='g', \
                        edgecolor='0.5')
        plt.plot(tPho['time'],tPho['mag'],c='g', lw='1')

        plt.plot(tPho['time'][bWeird], tPho['mag'][bWeird], 'rx', zorder=25)
        plt.xlabel('MJD')
        plt.ylabel('Magn')

    if trimWeird:
        tPho = tPho[~bWeird]
        tPho.meta['trimOoO'] = int(np.sum(bWeird))
        

    # return the table with MJDCalc
    return tPho

def estUncty92(tIn=Table(), unctAvg=0.04, magOffset = 16.08):

    """Adds an uncertainty column for forwards compatibility"""

    # 2018-07-20 - for the moment, 1992 is also the only data-table
    # for which the magnitudes are relative. So we add the offset here
    # too.

    if not 'mag' in tIn.colnames:
        return tIn

    if np.median(tIn['mag']) < 10.0: # if we haven't already corrected!
        tIn['mag'] = tIn['mag'] + magOffset

    # now we add the uncty
    tIn['magErr'] = Column(np.repeat(unctAvg, len(tIn)) )

    return tIn

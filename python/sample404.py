#
# sample404.py
#

# Started 2019-06-19 by WIC

# Purpose: sample fake lightcurve for V404 Cyg to test the hypothesis
# that our 201x datasets are drawn from the same underlying
# distribution as the Zurita et al. 2004 datasets to which we have
# access.

# 2019-06-19 11:00 pm - to add: 
#
# (i)  read in the Z04 lightcurve with ellipsoidals subtracted
#
# (ii) implement the figure of merit for lightcurve comparison

import os, sys, time
import numpy as np
import matplotlib.pylab as plt

from astropy.table import Table
from astropy.stats import LombScargle # for characterization
from sympy import factorint
import DELCgen

class FakeLC(object):

    """Object to hold fake lightcurve parameters and samples"""

    def __init__(self, filSamples='DIA2017.csv'):

        """INIT"""

        # power spectrum model
        self.modelChoice = 'BendingPL'
        self.methPSD = DELCgen.BendingPL
        self.parseModelChoice()

        # external file containing at least the times of observation
        self.filSamples = filSamples[:]
        self.keyObsTime = 'time' # from CJF's dia
        self.keyObsFlux = 'flux' # placeholder
        self.keyObsUnct = 'error' # placeholder
        self.useObsFlux = False
        self.useObsUncty = False
        self.tSamplesFactor = 1.0/1440. # scale factor for time file

        # time-chunks object for the observations
        self.ChunksObs = None
        self.minDtChunk = 0.01 # in days

        # the sample lightcurve
        self.tSample = np.array([])
        self.ySample = np.array([])
        self.eSample = np.array([])

        # the "blank" lightcurve to set up the sampling. This might be
        # time, rate, uncertainty from the 2017A dataset (we will want
        # to preserve the uncertainty if not the actual rates).
        self.yBlank = np.array([])
        self.eBlank = np.array([])

        # Model PSD parameters: defaults for self-testing
        #self.PSDpars = [10., 10., 0.00, 1., 50.]
        self.PSDpars = [100., 1000, 1.0, 1.0, 0.]

        # the following two are from tests of DELCgen's own
        # PSD-fitting on the raw 1992 data. I suspect the fitter isn't
        # working...
        #self.PSDpars = [0.10, 2.80, 19.1, 2.04, 14.02] # from a pass of fitter
        #self.PSDpars = [2.2e-6, 1.57, 1.57, 1.57, 3.3]

        #self.PSDpars = [5.01345204e-03,\
        #                    1.93937705e-02,\
        #                    -2.39479873e-02, \
        #                    2.21425182e+00, \
        #                    1.04358075e-03]


        # rates and uncertainties for the template lightcurve
        self.tTemplate = np.array([])
        self.yTemplate = np.array([])
        self.eTemplate = np.array([])

        # some default parameters
        self.defaultUncty = 0.01
        self.defaultMean = 16.

        # control parameters for sampled lightcurve
        self.sampleMean = 16.0
        self.sampleStd = 0.15
        self.sampleLen = 0
        self.sampleTbin = 100000 # make this large for fine sampling

        # DELCgen objects
        self.LCtemplate = None  # template from previous data
        self.LCblank = None  # blank lc for passing to Simulate_TK
        self.LCsample = None # the sampled lightcurve
        self.LCobs = None # separate LC object for observations if
                          # we're stacking this onto a longer blank

        # boolean for objects in the output sample that actually
        # correspond to observations
        self.bObs = np.array([])

        # a few parameters for picking the red noise simulation length
        self.rnlMax=105
        self.rnlMin=89
        self.rnlFactor = 100 # default rednoise factor. Updated by
                             # pickRNLfactor()

        # lomb-scargle variables for assessing the result
        self.lsPmin = 1. # day
        self.lsPmax = 80. # day
        self.lsLogbin = True
        self.lsNper = 1e3
        self.lsPer = np.array([])
        self.lsPow = np.array([])

        # control variable
        self.Verbose=True

    def genBlankLC(self, genUncty=False, yScatt=0.01, yMed=16.):

        """Initializes an LC object to be populated with the current
        simulation values"""

        nSample = np.size(self.tSample)

        # only override the uncertainties if not already present or
        # incompatible with the sample lightcurve size.
        if genUncty or np.size(self.eBlank) <> nSample:
            self.eBlank = np.random.normal(size=nSample)*yScatt

        # generate white-noise gaussian
        self.yBlank = np.random.normal(size=nSample)*self.eBlank + yMed
        
    def blankLCfromArrays(self):

        """Passes the blank output arrays (however we built them) up
        to a TKLC object"""

        self.LCblank = DELCgen.Lightcurve(time=np.copy(self.tSample), \
                                              flux=self.yBlank, \
                                              tbin=self.sampleTbin,\
                                              errors=self.eBlank)

    def genSampleTimes(self, nData=2000, mjdMin=0., mjdMax=16.):

        """Creates fake times for sampling"""

        self.tSample = self.genFakeTimes(nData, mjdMin, mjdMax)

    def lcObsFromFile(self):

        """Imports observation lightcurve from file"""

        self.LCobs = None
        if not os.access(self.filSamples, os.R_OK):
            if self.Verbose:
                print("FakeLC.lcObsFromFile WARN - cannot read path %s" \
                          % (self.filSamples))
            return

        tablObs = Table.read(self.filSamples)

        if not self.keyObsTime in tablObs.colnames:
            if self.Verbose:
                print("FakeLC.lcObsFromFile WARN - time key not in table: %s" \
                          % (self.keyObsTime))
            return

        # parse the observations table
        tObs = tablObs[self.keyObsTime]*self.tSamplesFactor
        
        # default unctys to be replaced if appropriate
        eObs = np.repeat(self.defaultUncty, np.size(tObs))
        if self.useObsUncty and self.keyObsUnct in tObs.colnames:
            eObs = tablObs[self.keyObsUnct]

        yObs = np.random.normal(size=np.size(tObs))*eObs + self.defaultMean
        if self.useObsFlux and self.keyObsFlux in tObs.colnames:
            yObs = tObs[self.keyObsFlux]

        # now we have the time, flux, error, create the "Obs"
        # lightcurve object
        self.LCobs = DELCgen.Lightcurve(tObs, yObs, self.sampleTbin, eObs)    

    def findObsChunks(self):

        """Finds the chunks in the LCobs object"""

        try:
            timesObs = self.LCobs.time
        except:
            if self.Verbose:
                print("FakeLC.findObsChunks WARN - problem with LCobs.time")
            return

        self.ChunksObs = TimeChunks(timesObs, self.minDtChunk, runOnInit=True)

    def buildOvertimes(self, tBuffer=0.5, nFac=5):

        """Creates a hybrid lightcurve object, with the observed
        sampling times meshed with a larger set of sampiing times that
        fill the gaps.

        tBuffer = number of days early and late to build the lightcurve

        nFac = if > 50, the number of datapoints in the larger
        list. Otherwise, the number of datapoints to simulate as a
        multiple of the observation data-length

        """

        # find the minmax times for the larger dataset
        try:
            timesObs = self.LCobs.time
            fluxObs =  self.LCobs.flux
            unctyObs = self.LCobs.errors
        except:
            if self.Verbose:
                print("FakeLC.buildOvertimes WARN - problem with LCobs.time")
            return

        tMin = np.min(timesObs) - tBuffer
        tMax = np.max(timesObs) + tBuffer
        
        nSim = np.copy(nFac)
        if nFac < 50:
            nSim = np.size(timesObs)*nFac

        # median uncty and y values, taken from the observation object
        eMed = np.median(unctyObs)
        yMed = np.median(fluxObs)

        tLarge = np.linspace(tMin, tMax, nSim, endpoint=True)
        eLarge = np.repeat(eMed, np.size(tLarge))
        yLarge = np.random.normal(size=np.size(tLarge))*eLarge + yMed
                           
        # now fuse the large with the observations. We do this for the
        # time, flux, error for both objects. 
        bLargeInChunk = self.ChunksObs.timesInChunks(tLarge)

        # now build the combined arrays and argsort by times
        tCombo = np.hstack(( tLarge[~bLargeInChunk], timesObs ))
        yCombo = np.hstack(( yLarge[~bLargeInChunk], fluxObs ))
        eCombo = np.hstack(( eLarge[~bLargeInChunk], unctyObs ))
        
        lSor = np.argsort(tCombo)
        
        # pass this up to the blank LC object and set the chunk
        # indices
        self.LCblank = DELCgen.Lightcurve(tCombo[lSor], \
                                              yCombo[lSor], \
                                              tbin=self.sampleTbin, \
                                              errors=eCombo[lSor])

        # now we have this, set a boolean with the observations in the
        # chunk
        self.bObs = self.ChunksObs.timesInChunks(self.LCblank.time)

    def templateLCfromArrays(self):

        """Creates the DELCgen template lightcurve object from input
        arrays"""

        self.LCtemplate = DELCgen.Lightcurve(\
            time=self.tTemplate, \
                flux=self.yTemplate, \
                tbin=self.sampleTbin, \
                errors=self.eTemplate)

    def genTemplateLightcurve(self, nData=1500, \
                                  mjdMin=0., mjdMax=16., \
                                  yMed=16., yScatt=0.05):

        """Generates template lightcurve as arrays from arguments"""

        self.genTemplateTimes(nData, mjdMin, mjdMax)
        self.genTemplateFlux(yMed, yScatt)
        self.templateLCfromArrays()

    def genTemplateTimes(self, nData=1500, mjdMin=0., mjdMax=16.):

        """Creates template lightcurve times"""

        self.tTemplate = self.genFakeTimes(nData, mjdMin, mjdMax)

    def genGaussRates(self, nData=1000, yMed=16., yStd=0.05):

        """Generates gaussian random numbers w/ uncertainties. To be
        used initializing templates."""

        yScatt = np.repeat(yStd, nData)
        yValus = np.random.normal(size=nData) * yScatt + yMed

        return yScatt, yValus

    def genTemplateFlux(self, yMed=16., yStd=0.05):

        """Creates template lightcurve magnitudes. Array size
        inherited from self.tTemplate, so make sure that's set
        first."""

        self.yTemplate, self.eTemplate = \
            self.genGaussRates(np.size(self.tTemplate), yMed, yStd)

    def genFakeTimes(self, nData=1500, \
                         mjdMin=0.0, mjdMax=16.0):

        """Generates fake MJD for self-testing purposes"""
        
        return np.sort(np.random.uniform(size=nData, low=mjdMin, high=mjdMax))

    def getStatsFromTemplate(self):

        """Get the statistics for the sample lightcurve from the
        template object"""

        self.sampleMean = self.LCtemplate.mean
        self.sampleStd = self.LCtemplate.std
        self.sampleLen = self.LCtemplate.length

    def parseModelChoice(self):

        """Selects the chosen PSD model from DELCgen"""

        modelDefault='BendingPL'
        if not hasattr(DELCgen, self.modelChoice):
            print("FakeLC.parseModelChoice WARN - model %s not found." \
                      % (self.modelChoice))
            print("FakeLC.parseModelChoice WARN - defaulting to %s" \
                      % (modelDefault))
            self.methPSD = getattr(DELCgen, modelDefault)
            return

        self.methPSD = getattr(DELCgen, self.modelChoice)

    def sampleNoiseModel(self):

        """Samples the BPL into a new lightcurve object"""

        self.LCsample = None
        self.LCsample = \
            DELCgen.Simulate_TK_Lightcurve(\
            self.methPSD, \
                self.PSDpars, \
                lightcurve=self.LCblank, \
                RedNoiseL=self.rnlFactor, \
                mean=self.sampleMean, \
                std=self.sampleStd,\
                length=self.LCblank.length, \
                tbin=self.LCblank.tbin)

        # interesting bug: the errors don't seem to be sent through to
        # the simulated object. We'll graft them on here.
        self.LCsample.errors = np.copy(self.LCblank.errors)

    def lsNoiseModel(self, clobberPer=False):

        """Convenience-method to draw the lomb-scargle on the current
        sample"""

        try:
            tSampl = self.LCsample.time
            ySampl = self.LCsample.flux
            eSampl = self.LCsample.errors
        except:
            if self.Verbose:
                print("FakeLC.lsNoiseModel WARN - problem with noise model")
            return

        if np.size(self.lsPer) < 1 or clobberPer:
            self.setupLombScargle()

        self.lsPow = LombScargle(tSampl, ySampl, eSampl).power(1.0/self.lsPer)
        

    def pickRNLfactor(self):

        """Picks the appropriate rednoise factor to make simulations
        most speedy. Refactored a little from ASMLS.py developed for
        Dage et al. (2019 MNRAS)"""

        nData = 0 
        if len(self.tSample) > 1:
            nData = np.size(self.tSample)
        else:
            if len(self.LCblank.time) > 1:
                nData = np.size(self.LCblank.time)

        if nData < 1:
            if self.Verbose:
                print("FakeLC.pickRNLfactor WARN - sample times not yet set")
            return

        # trial RN values, set up the big factors and the lens
        rnlVals = np.arange(self.rnlMax, self.rnlMin, -1)

        # initialize the big-factor and len arrays by copy from the
        # rnlVals
        bigFacs = rnlVals * 0. + 1.0e7
        tklens = rnlVals * 0.

        for iRNL in range(np.size(rnlVals)):
            thisRNL = rnlVals[iRNL]
            tkLength = (thisRNL * nData) +1
            bigFacs[iRNL] = max(factorint(tkLength))
            tklens[iRNL] = tkLength

            if self.Verbose:
                print("FakeLC.pickRNLfactor INFO - RNL, bigFac %i, %i at length = %i" % (thisRNL, bigFacs[iRNL], tkLength))

        # Having populated the factors from the lengths, find the
        # shortest and set this as an instance-level property.
        iShortest = np.argmin(bigFacs)
        self.rnlFactor = rnlVals[iShortest]

        if self.Verbose:
            print("FakeLC.pickRNLfactor INFO - picked factor %i" \
                      % (self.rnlFactor))

    def setupLombScargle(self, scaleLimits=True):

        """Sets up lomb-scargle variables. Optionally scales the LS
        limits appropriately for the time sample simnulated"""

        if scaleLimits:
            self.lsPmax = 2.0*np.max(self.LCblank.time)
            self.lsPmin = self.lsPmax / 1000.

        self.lsPer = np.logspace(np.log10(self.lsPmin), \
                                     np.log10(self.lsPmax), \
                                     self.lsNper, endpoint=True)
        
    def showLC(self, figname='testLC.png'):

        """Debug routine - shows the simulated lightcurve"""

        try:
            tSampl = self.LCsample.time
            ySampl = self.LCsample.flux
            eSampl = self.LCsample.errors
        except:
            if self.Verbose:
                print("FakeLC.showLC WARN - problem getting the sampled arrays")
            return

        # ensure a boolean is set corresponding to the actual
        # observations. If we're not being clever and breaking this
        # up, use all the objects.
        bSho = self.bObs[:]
        if len(bSho) < 1:
            bSho = np.repeat(True, np.size(tSampl))

        # calculate the lomb scargle here
        if np.size(self.lsPer) < 1:
            self.setupLombScargle()
            
        # compute the lomb-scargle both in and out of observation
        lsFreq = 1.0/self.lsPer
        lsSho = LombScargle(tSampl[bSho], ySampl[bSho], \
                                eSampl[bSho]).power(lsFreq)

        # do the Lomb-Scargle for all the points
        lsAll = LombScargle(tSampl, ySampl, \
                                eSampl).power(lsFreq)

        # set colors upfront
        colorSho = 'k'
        colorOut = '0.5'

        fig1 = plt.figure(1)
        fig1.clf()

        fig1.subplots_adjust(hspace=0.3)

        ax1 = fig1.add_subplot(211)
        
        dumScatt = ax1.errorbar(tSampl[bSho], ySampl[bSho], \
                                    yerr=eSampl[bSho], ls='none', \
                                    ms=0.5, alpha=0.5, marker='o', zorder=2, \
                                    color=colorSho, ecolor=colorSho)

        dumOut = ax1.errorbar(tSampl[~bSho], ySampl[~bSho], \
                                  yerr=eSampl[~bSho], ls='none', \
                                  ms=0.5, alpha=0.5, marker='o', zorder=2, \
                                  color=colorOut, ecolor=colorOut)


        #dumPlot = ax1.plot(tSampl, ySampl, alpha=0.15, zorder=1, c='b')
        ax1.set_xlabel('MJD')
        ax1.set_ylabel('Flux')

        # now show the LS of this sample. We might or might not use
        # two panels for this...
        if np.sum(~bSho) < 1:
            ax2 = fig1.add_subplot(212)
        else:
            ax2 = fig1.add_subplot(223)

        # do the LS for the "observations"
        dumLS = ax2.loglog(self.lsPer, lsSho, color=colorSho, \
                               lw=1)
        ax2.set_xlabel('Period (d)')
        ax2.set_ylabel('LS power')
        ax2.set_xlim(self.lsPmin, self.lsPmax)

        # if the entire observation set is different, plot that too.
        if np.size(lsAll) > 0:
            ax3 = fig1.add_subplot(224, sharey=ax2)
            dumLS3 = ax3.loglog(self.lsPer, lsAll, color=colorOut, \
                                    lw=1)
            ax3.set_xlabel('Period (d)')
            #ax3.set_ylabel('LS power')
            ax3.set_xlim(self.lsPmin, self.lsPmax)

        # save the figure to disk
        fig1.savefig(figname, rasterized=False)

class TimeChunks(object):

    """Convenience object to identify chunks in a time-series"""

    def __init__(self, times=np.array([]), chunkMin=0.01, \
                     runOnInit=True):

        self.chunkMin=chunkMin
        self.times = np.copy(times)

        self.chunkIDs = np.array([])

        if runOnInit:
            self.findChunks()
            self.assignChunkIDs()

    def findChunks(self):

        """Finds the chunks"""
        
        bStarts = self.times - np.roll(self.times, 1) > self.chunkMin
        bEnds = np.roll(self.times, -1) - self.times > self.chunkMin

        # don't forget the endpoints
        bStarts[0] = True
        bEnds[-1] = True

        self.chunkStarts = self.times[bStarts]
        self.chunkEnds = self.times[bEnds]

    def assignChunkIDs(self):

        """Assign every point to a chunk"""
    
        self.chunkIDs = np.repeat(0, np.size(self.times))
        for iChunk in range(np.size(self.chunkStarts)):

            # conditional must be .le. and .ge. not .lt. or .gt.
            bChunkAbove = self.times >= self.chunkStarts[iChunk] 
            bChunkBelow = self.times <= self.chunkEnds[iChunk]
            bChunk = bChunkAbove & bChunkBelow
                                  
            self.chunkIDs[bChunk] = iChunk

    def timesInChunks(self, tTest=np.array([])):

        """For input times, return True for items that are within a
        chunk"""

        bWithin = np.repeat(False, np.size(tTest))

        # I don't know of a quick array-based way to do this, so let's
        # loop through the chunks instead.

        for iChunk in range(np.size(self.chunkStarts)):
            bHi = tTest >= self.chunkStarts[iChunk]
            bLo = tTest <= self.chunkEnds[iChunk]
            bThis = bHi & bLo

            bWithin[bThis] = True

        return bWithin

    def mergeWithChunks(self, tInpu=np.array([])):

        """Given an input time array, replaces any times within a
        chunk with the object times that are actually in the chunk"""

        bInChunk = self.timesInChunks(tInpu)
        
        # construct the output by stacking input times not in a chunk
        # with the times in the chunks, then sorting by time
        tCombo = np.hstack((tInpu[~bInChunk], self.times))
        return np.sort(tCombo)
            

def testFakeLC(sampleFil=''):

    """Tests the functionality of the fake lc generator"""

    FLC = FakeLC()

    # load sample times to send to the sampler
    if os.access(sampleFil, os.R_OK):
        tTimes = Table.read(sampleFil)
        FLC.tSample = tTimes['time']/1440.

        # test our chunkfinder
        TC = TimeChunks(FLC.tSample)

        # generate a uniformly sampled array over the whole time
        # baseline so that we can test our distinguisher
        tUnif = np.linspace(0., 5.5, 1000)
        tJoin = TC.mergeWithChunks(tUnif)
        yJoin = np.repeat(1.2, np.size(tJoin))
        bJoin = TC.timesInChunks(tJoin)
                          
        
        # generate a time array a little longer
        tRand = np.random.uniform(0., 5.5, 1000)
        yRand = np.repeat(0.90, np.size(tRand))
        bRand = TC.timesInChunks(tRand)
        

        plt.figure(2)
        plt.clf()
        yTimes = np.repeat(1.0, np.size(FLC.tSample))
        yStarts = np.repeat(1.05, np.size(TC.chunkStarts))
        yEnds = np.repeat(1.05, np.size(TC.chunkEnds))
        
        plt.scatter(FLC.tSample, yTimes, c=TC.chunkIDs, s=36, \
                        cmap='flag')
        
        for iChunk in range(np.size(TC.chunkStarts)):
            plt.plot(TC.chunkStarts[iChunk]*np.array([1.0, 1.0]), \
                         [0.5, 1.5], color='g')
            plt.plot(TC.chunkEnds[iChunk]*np.array([1.0, 1.0]), \
                         [0.5, 1.5], color='r', ls='--')


        plt.plot(TC.chunkStarts, yStarts, 'gs')
        plt.plot(TC.chunkEnds, yEnds, 'rs')
        plt.ylim(0.5,1.5)
                           
        # now plot the test array with and without the boolean
        plt.plot(tRand[bRand], yRand[bRand], 'ko')
        plt.plot(tRand[~bRand], yRand[~bRand]-0.005, marker='o', color='0.8', \
                     ms=3, ls='None')

        # also plot the joined dataset
        plt.plot(tJoin[bJoin], yJoin[bJoin], 'ko')
        plt.plot(tJoin[~bJoin], yJoin[~bJoin]-0.005, marker='o', color='0.8', \
                     ms=3, ls='None')

        return

    else:
        FLC.genSampleTimes(2000, mjdMax=6)

    # Fake data to simiulate an "old" lightcurve
    # FLC.genTemplateLightcurve()

    # generate blank template for sampling with BPL, using the input times
    FLC.genBlankLC()
    FLC.blankLCfromArrays()

    # now create a single BPL sample
    FLC.pickRNLfactor()

    print("Sampling...")
    FLC.sampleNoiseModel()

    # do LS on the noise model
    print("Calculating LS...")
    FLC.setupLombScargle(scaleLimits=True)
    FLC.lsNoiseModel()
    
    # show the results
    print("Plotting...")
    FLC.showLC()

    print FLC.LCsample.flux - FLC.LCblank.flux


def testDirect(filIn='lc92raw.txt', yStd=0.05, tBin=1e4):

    """More direct method - try using the DELCgen methods to
    characterize the LC as well"""

    
    tOld = Table.read(filIn, comment='%', format='ascii')

    tData = tOld['col1']
    yData = tOld['col2']
    eData = np.repeat(yStd, np.size(tData))

    lcData = DELCgen.Lightcurve(tData, yData, tBin, eData)

    print "Fitting PSD..."
    lcData.Fit_PSD(model=DELCgen.BendingPL, verbose=True)

    print "Fitting PDF..."
    lcData.Fit_PDF()

    print lcData.psdFit
    print lcData.psdModel


def testCombinedSample():

    """Tests building a combined sample with observations and a wider
    time baseline to distinguish observations from the wider
    sample."""

    FLC = FakeLC('DIA2017.csv')
    FLC.lcObsFromFile()
    FLC.findObsChunks()
    FLC.buildOvertimes()

    # since we now have a larger sample, we don't need such a large
    # RNLfactor. Use a smaller range.
    #FLC.rnlMin=10
    #FLC.rnlMax=20
    FLC.pickRNLfactor()
    
    # Sample the noise model...
    FLC.sampleNoiseModel()
    
    # ... and show the lightcurve
    FLC.showLC()

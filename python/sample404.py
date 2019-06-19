#
# sample404.py
#

# Started 2019-06-19 by WIC

# Purpose: sample fake lightcurve for V404 Cyg to test the hypothesis
# that our 201x datasets are drawn from the same underlying
# distribution as the Zurita et al. 2004 datasets to which we have
# access.

import os, sys, time
import numpy as np
import matplotlib.pylab as plt

from astropy.stats import LombScargle # for characterization
from sympy import factorint
import DELCgen

class FakeLC(object):

    """Object to hold fake lightcurve parameters and samples"""

    def __init__(self, psdModel='BendingPL'):

        """INIT"""

        # power spectrum model
        self.modelChoice = psdModel[:]
        self.methPSD = DELCgen.BendingPL
        self.parseModelChoice()

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
        self.PSDpars = [100., 1000, 1.5, 1.5, 0.]

        #self.PSDpars = [5.01345204e-03,\
        #                    1.93937705e-02,\
        #                    -2.39479873e-02, \
        #                    2.21425182e+00, \
        #                    1.04358075e-03]


        # rates and uncertainties for the template lightcurve
        self.tTemplate = np.array([])
        self.yTemplate = np.array([])
        self.eTemplate = np.array([])

        # control parameters for sampled lightcurve
        self.sampleMean = 16.0
        self.sampleStd = 0.3
        self.sampleLen = 0
        self.sampleTbin = 100000 # make this large for fine sampling

        # DELCgen objects
        self.LCtemplate = None  # template from previous data
        self.LCblank = None  # blank lc for passing to Simulate_TK
        self.LCsample = None # the sampled lightcurve

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

        if len(self.tSample) < 1:
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
            tkLength = (thisRNL * np.size(self.tSample)) +1
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
        

    def showLC(self):

        """Debug routine - shows the simulated lightcurve"""

        try:
            tSampl = self.LCsample.time
            ySampl = self.LCsample.flux
            eSampl = self.LCsample.errors
        except:
            if self.Verbose:
                print("FakeLC.showLC WARN - problem getting the sampled arrays")
            return

        fig1 = plt.figure(1)
        fig1.clf()
        ax1 = fig1.add_subplot(211)

        dumScatt = ax1.errorbar(tSampl, ySampl, yerr=eSampl, ls='none', \
                                    ms=2, alpha=0.5, marker='o', zorder=2)
        dumPlot = ax1.plot(tSampl, ySampl, alpha=0.3, zorder=1, c='b')
        ax1.set_xlabel('MJD')
        ax1.set_ylabel('Flux')

        # now show the LS of this sample
        ax2 = fig1.add_subplot(212)
        dumLS = ax2.loglog(self.lsPer, self.lsPow)
        ax2.set_xlabel('Period (d)')
        ax2.set_ylabel('LS power')
        ax2.set_xlim(self.lsPmin, self.lsPmax)

def testFakeLC():

    """Tests the functionality of the fake lc generator"""

    FLC = FakeLC()

    # Fake data to simiulate an "old" lightcurve
    FLC.genTemplateLightcurve()

    # generate blank template for sampling with BPL
    FLC.genSampleTimes(2000, mjdMax=8)
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

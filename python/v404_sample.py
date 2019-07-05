#
# v404_sample.py
#

# Started 2019-06-19 by WIC

# Purpose: sample fake lightcurve for V404 Cyg to test the hypothesis
# that our 201x datasets are drawn from the same underlying
# distribution as the Zurita et al. 2004 datasets.

# 2019-06-22: testCompact() currently performs the simulation using
# observation data to set the observation times and a template
# lightcurve to poll for the PSD characteristics. This is the one I'm
# developing at the moment so is a good place to start.

# 2019-06-22 to add:

# (i) Apply the identical FoM to the template lightcurve and/or the
# observation lightcurve, propagating the results through into the
# output file (so that everything is in the same location);
#
# (ii) A simple class to perform the output comparison. This would do
# the nice plots for the paper. AMB will likely have thoughts on this!

import os, sys, time
import copy
import numpy as np
import matplotlib.pylab as plt
plt.style.use('seaborn-whitegrid') # could put this into the object

from astropy.table import Table
from astropy.stats import LombScargle # for characterization
from sympy import factorint
import DELCgen

# for dumping the array of trial statistics to disk
from astropy.io import fits 

# For interpolation of existing data into a sample
from scipy import interpolate

class FakeLC(object):

    """Object to hold fake lightcurve parameters and samples"""

    def __init__(self, filSamples='DIA2017.csv', \
                     filTemplate='92Binned.fits', \
                     keyObsTime='time', \
                     keyObsFlux='flux', \
                     keyObsUnct='error', \
                     tSamplesFactor=1.0/1440.):

        # Samples file default keywords and samples factor are set
        # appropriately for CJF's DIA analysis
        
        # power spectrum model
        self.modelChoice = 'BendingPL'
        self.methPSD = DELCgen.BendingPL
        self.parseModelChoice()

        # perturb the samples with measurement uncertainty after the
        # generation of the power law noise
        self.pertByMeasureUncty = True

        # optional file containing prior observations we want to
        # reproduce
        self.filTemplate = filTemplate[:]

        # optional file containing at least the times of observation
        self.filSamples = filSamples[:]
        self.keyObsTime = keyObsTime[:] # 'time' # from CJF's dia
        self.keyObsFlux = keyObsFlux[:] # 'flux' # placeholder
        self.keyObsUnct = keyObsUnct[:] # 'error' # placeholder
        self.useObsFlux = False
        self.useObsUncty = False
        self.tSamplesFactor = np.float(tSamplesFactor)

        # 2019-07-05 (WIC) - I remembered we think some of the Zurita
        # files list 3-sigma for the uncertainties. This factor allows
        # that correction. The errors will be DIVIDED by these
        # quantities on import.
        self.errDivisorObs = 1.
        self.errDivisorTemplate = 3. # WATCHOUT - DEFAULT CHANGED.

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

        # default reference mag and flux for conversion. (2019-06-22
        # I'm undecided about whether to set this at the FakeLC level
        # or to set this object-by-object.)
        self.refMag = 16.2
        self.refFlux = 10000.

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

    def lcTemplateFromFile(self, refMag=-99., refFlux=-99., isMag=True):

        """Populates template arrays and lightcurve object from
        file."""

        # inherit the reference flux and magnitude from the instance
        # if not set as arguments
        if refMag < -90:
            refMag = np.copy(self.refMag)

        if refFlux < -90:
            refFlux = np.copy(self.refFlux)

        # Also sets the reference magnitude and flux at both the
        # instance level for FakeLC and for the lightcurve object.

        if not os.access(self.filTemplate, os.R_OK):
            if self.Verbose:
                print("FakeLC.lcTemplateFromFile WARN - cannot read path %s" \
                          % (self.filTemplate))
            return

        tTempl = Table.read(self.filTemplate)

        # transfer the data to the separate arrays and to the
        # LCtemplate object
        self.tTemplate = tTempl['tBin']
        self.yTemplate = tTempl['fBin']
        self.eTemplate = tTempl['uBin'] / self.errDivisorTemplate

        self.LCtemplate = DELCgen.Lightcurve(\
            self.tTemplate, self.yTemplate, self.sampleTbin, self.eTemplate)    

        # Set the reference quantities
        self.LCtemplate.isFlux = np.logical_not(isMag)
        self.LCtemplate.refFlux = refFlux
        self.LCtemplate.refMag = refMag

        # Ensure the FakeLC instance has these values too
        self.refFlux = refFlux
        self.refMag = refMag

    def lcObsFromFile(self, refMag=-99, refFlux=-99, isMag=True):

        """Imports observation lightcurve from file"""

        # inherit the reference values from the instance
        if refMag < -90:
            refMag = np.copy(self.refMag)

        if refFlux < -90:
            refFlux = np.copy(self.refFlux)

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
        if self.useObsUncty and self.keyObsUnct in tablObs.colnames:
            eObs = tablObs[self.keyObsUnct] / self.errDivisorObs

        yObs = np.random.normal(size=np.size(tObs))*eObs + self.defaultMean
        if self.useObsFlux and self.keyObsFlux in tablObs.colnames:
            yObs = tablObs[self.keyObsFlux]

        # now we have the time, flux, error, create the "Obs"
        # lightcurve object
        self.LCobs = DELCgen.Lightcurve(tObs, yObs, self.sampleTbin, eObs)    

        # Set the reference quantities
        self.LCobs.isFlux = np.logical_not(isMag)
        self.LCobs.refFlux = np.float(refFlux)
        self.LCobs.refMag = np.float(refMag)

    def findObsChunks(self):

        """Finds the chunks in the LCobs object"""

        try:
            timesObs = self.LCobs.time
        except:
            if self.Verbose:
                print("FakeLC.findObsChunks WARN - problem with LCobs.time")
            return

        self.ChunksObs = TimeChunks(timesObs, self.minDtChunk, runOnInit=True)

    def lcsAsFlux(self):

        """Ensures the lightcurve objects are expressed as flux"""

        # ensure the same reference flux is set for everything. (Note
        # that it's really the uncertainties we care about here in the
        # conversions)
        refFlux = self.LCtemplate.refFlux

        # This can be looped, but we write out here for clarity.
        if not self.LCobs.isFlux:
            self.LCobs = self.magToFlux(self.LCobs, refFlux)

        if not self.LCblank.isFlux:
            self.LCblank = self.magToFlux(self.LCblank, refFlux)

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

        # Take the median flux and uncertainty from the template
        # object, unless otherwise set
        if self.useObsUncty:
            LCforUnct = self.LCobs
        else:
            LCforUnct = self.LCtemplate
            
        eMed = np.median(LCforUnct.errors)
        eStd = np.std(LCforUnct.errors)

        print("INFO - eMed, eStd, isFlux:", eMed, eStd, \
                  LCforUnct.isFlux, LCforUnct.refMag)

        if self.useObsFlux:
            LCforFlux = self.LCobs
        else:
            LCforFlux = self.LCtemplate
        yMed = np.median(LCforFlux.flux)

        #eMed = np.median(self.LCtemplate.errors)
        #yMed = np.median(self.LCtemplate.
        #eMed = np.median(unctyObs)
        #yMed = np.median(fluxObs)

        # populate the larger sample including the gaps
        tLarge = np.linspace(tMin, tMax, nSim, endpoint=True)
        nLarge = np.size(tLarge)
        eLarge = np.repeat(eMed, nLarge)
        eLarge = np.random.normal(size=nLarge)*eStd + eMed
        yLarge = np.random.normal(size=nLarge)*eLarge + yMed
                           
        # What we use for the observation intervals depends on whether
        # we are using the observed yvalues and errors. If we are
        # using them, simply copy them in. Otherwise, generate new
        # ones using the parameters we copied from the template
        # lightcurve
        nObs = np.size(timesObs)
        if self.useObsUncty:
            obsErr = np.copy(unctyObs)
        else:
            obsErr = np.random.normal(size=nObs)*eStd + eMed

        if self.useObsFlux:
            obsFlux = np.copy(fluxObs)
        else:
            obsFlux = np.random.normal(size=nObs)*obsErr + yMed

        # now fuse the large with the observations. We do this for the
        # time, flux, error for both objects. 
        bLargeInChunk = self.ChunksObs.timesInChunks(tLarge)

        # now build the combined arrays and argsort by times
        # tCombo = np.hstack(( tLarge[~bLargeInChunk], timesObs ))
        # yCombo = np.hstack(( yLarge[~bLargeInChunk], fluxObs ))
        # eCombo = np.hstack(( eLarge[~bLargeInChunk], unctyObs ))

        tCombo = np.hstack(( tLarge[~bLargeInChunk], timesObs ))
        yCombo = np.hstack(( yLarge[~bLargeInChunk], obsFlux ))
        eCombo = np.hstack(( eLarge[~bLargeInChunk], obsErr ))

        
        lSor = np.argsort(tCombo)
        
        # pass this up to the blank LC object and set the chunk
        # indices
        self.LCblank = DELCgen.Lightcurve(tCombo[lSor], \
                                              yCombo[lSor], \
                                              tbin=self.sampleTbin, \
                                              errors=eCombo[lSor])

        # the LCblank should by this point be in flux (2019-06-22:
        # should add some syntax to ensure this!)
        self.LCblank.isFlux = True
        self.LCblank.refFlux = self.refFlux # MAY WANT TO CHECK
        self.LCblank.refMag = self.refMag 

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

    def getTemplateStats(self, unctyCorr = True):

        """Get the statistics for the sample lightcurve from the
        template object. If unctyCorr is true, corrects for the (mean)
        measurement uncertainty."""

        # When simulating in sampleNoiseModel(), self.sampleMean and
        # self.sampleStd are used.

        # (Note that if we're dealing with red noise, the stddev
        # should be much larger than the average measurement
        # uncertainty. But we want to get this right.)

        # Update 2019-06-23 - this recomputes the statistics rather
        # than simply reporting them across, since DELCgen.Lightcurve
        # only computes mean and std on initialization.
        try:
            yFlux = self.LCtemplate.flux
        except:
            if self.Verbose:
                print("FakeLC.getTemplateStats WARN - LCtemplate problem")
            return

        # compute the mean and stddev
        self.sampleMean = np.mean(yFlux)
        self.sampleStd = np.std(yFlux)

        if unctyCorr:
            errMean = np.mean(self.LCtemplate.errors)
            self.sampleStd = np.sqrt(np.var(yFlux) - errMean**2)

        # report the stats found to terminal
        if self.Verbose:
            print("FakeLC.getTemplateStats INFO - template mean %.2f, stddev %.3f" % (self.sampleMean, self.sampleStd))

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

        # perturb the results by measurement uncertainty
        if self.pertByMeasureUncty:
            ePert = np.random.normal(size=np.size(self.LCsample.errors)) \
                * self.LCsample.errors
            self.LCsample.flux += ePert

        # copy the units-quantities from blank
        self.LCsample.isFlux = self.LCblank.isFlux
        self.LCsample.refMag = self.LCblank.refMag
        self.LCsample.refFlux = self.LCblank.refFlux

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
            self.lsPmax = 2.0*(np.max(self.LCblank.time) \
                                   - np.min(self.LCblank.time))
            self.lsPmin = self.lsPmax / 1000.

        self.lsPer = np.logspace(np.log10(self.lsPmin), \
                                     np.log10(self.lsPmax), \
                                     self.lsNper, endpoint=True)
        
    def magToFlux(self, LCobj=None, refFlux=-99, asCopy=False):

        """Convenience-method to convert a lightcurve to flux. Uses
        the LC object's own attributes to get the reference
        magnitude"""

        # 2019-06-23 note to self: "Flux" is really "counts" or "count
        # rate". We could use astropy's unit conversion functionality
        # so that we don't have to have an arbitrary input

        try:
            yMag = LCobj.flux
            eMag = LCobj.errors
        except:
            if self.Verbose:
                print("FakeLC.magToFlux WARN - problem with LC object")
            return None

        # return the altered original or a copy?
        if asCopy:
            LCret = copy.deepcopy(LCobj)
        else:
            LCret = LCobj

        # if the output is already in flux, do nothing
        if hasattr(LCret, 'isFlux'):
            if LCret.isFlux:
                return LCret

        # Ensure the ref flux is set somewhere
        if refFlux < 0:
            refFlux = np.copy(LCobj.refFlux)
        else:
            refFlux = np.copy(refFlux)

        # for the moment, we take the refMag and refFlux from the LC
        # object for internal consistency.
        deltaMag = LCobj.refMag - yMag

        # print("magToFlux INFO - refMag, refFlux", LCobj.refFlux, LCobj.refMag)

        # convert the magnitude to flux, use the Taylor approximation
        # for the uncertainties
        yFlux = refFlux * 100.0**(deltaMag/5.0)
        eFlux = yFlux * eMag / 1.086

        # Update the attributes
        LCret.flux = yFlux
        LCret.errors = eFlux

        # we recompute the statistics, since LightCurve seems to
        # compute them once on initialization. We do the same ops here
        # as the first few lines of Lightcurve.__init__ in DELCgen.py.
        LCret.mean = np.mean(LCobj.flux)
        LCret.std = np.std(LCobj.flux) 

        # add attributes to the LC object (note that the default LC
        # object doesn't have these by default. We might consider
        # adding them whenever we generate these objects.)
        LCret.isFlux=True
        LCret.refFlux=np.float(refFlux)

        return LCret

    def fluxToMag(self, LCobj=None, refMag=-99, asCopy=False):

        """Converts an LC object in flux to magnitude"""

        try:
            yFlux = LCobj.flux
            eFlux = LCobj.errors
        except:
            if self.Verbose:
                print("FakeLC.fluxToMag WARN - problem with LC object")
            return None

        # returning as a copy?
        if asCopy:
            LCret = copy.deepcopy(LCobj)
        else:
            LCret = LCobj

        # If the output is ALREADY in magnitudes, do nothing
        if hasattr(LCret, 'isFlux'):
            if not LCret.isFlux:
                return LCret

        # ensure the ref mag is set somewhere
        if refMag < 0:
            refMag = np.copy(LCobj.refMag)
        else:
            refMag = np.copy(refMag)

        refFlux = np.copy(LCobj.refFlux)
        yMag = refMag - 2.5*np.log10(yFlux / refFlux)
        eMag = 1.086 * eFlux / refFlux 

        LCret.flux = yMag
        LCret.errors = eMag
        
        # recalculate (rather than transform) the mean and std
        LCret.mean = np.mean(LCret.flux)
        LCret.std = np.std(LCret.flux)

        LCret.isFlux = False
        LCret.refMag = np.float(refMag)

        return LCret
        
            
    def showLC(self, figname='testLC.png', showMag=True):

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

        # log-bin them
        nBin = 10
        logShoPer, logShoPow, _ = binData(np.log10(self.lsPer), \
                                              np.log10(lsSho), nBin, \
                                              useMedian=True)
        logAllPer, logAllPow, _ = binData(np.log10(self.lsPer), \
                                              np.log10(lsAll), nBin, \
                                              useMedian=True)
        
        # set colors upfront
        colorSho = 'k'
        colorOut = '0.5'

        # plot quantities
        nRows = 2
        nCols = 2
        figSz = (7,6)

        if self.LCsample.isFlux and showMag:
            nRows = 3
            figSz = (7,8)

        fig1 = plt.figure(1)
        fig1.set_size_inches(*figSz, forward=True)
        fig1.clf()

        fig1.subplots_adjust(hspace=0.35)

        ax1 = fig1.add_subplot(nRows, 1, 1)
        
        dumScatt = ax1.errorbar(tSampl[bSho], ySampl[bSho], \
                                    yerr=eSampl[bSho], ls='none', \
                                    ms=1, alpha=0.5, marker='o', zorder=2, \
                                    color=colorSho, ecolor=colorSho)

        dumOut = ax1.errorbar(tSampl[~bSho], ySampl[~bSho], \
                                  yerr=eSampl[~bSho], ls='none', \
                                  ms=0.5, alpha=0.5, marker='o', zorder=2, \
                                  color=colorOut, ecolor=colorOut)


        #dumPlot = ax1.plot(tSampl, ySampl, alpha=0.15, zorder=1, c='b')
        ax1.set_xlabel('MJD (d)')
        ax1.set_ylabel('Flux')

        # now show the LS of this sample. We might or might not use
        # two panels for this...
        if np.sum(~bSho) < 1:
            ax2 = fig1.add_subplot(nRows, 1, nRows)
        else:
            ax2 = fig1.add_subplot(nRows, nCols, nRows*nCols - 1)

        # do the LS for the "observations"
        dumLS = ax2.loglog(self.lsPer, lsSho, color=colorSho, \
                               lw=1)
        ax2.set_xlabel('Period (d)')
        ax2.set_ylabel('LS power')
        ax2.set_xlim(self.lsPmin, self.lsPmax)

        # overplot the binned LS
        dumBin = ax2.loglog(10.0**logShoPer, 10.0**logShoPow, 'r-', \
                                zorder=5)

        # if the entire observation set is different, plot that too.
        if np.size(lsAll) > 0:
            ax3 = fig1.add_subplot(nRows, nCols, nRows*nCols, sharey=ax2)
            dumLS3 = ax3.loglog(self.lsPer, lsAll, color=colorOut, \
                                    lw=1)
            ax3.set_xlabel('Period (d)')
            #ax3.set_ylabel('LS power')
            ax3.set_xlim(self.lsPmin, self.lsPmax)

            # overplot the binned LS
            dumBin = ax3.loglog(10.0**logAllPer, 10.0**logAllPow, 'r-', \
                                    zorder=5)


        # Finally, if asked, show the flux converted to magnitude
        # using the reference flux and magnitude assigned to the
        # LCsample object
        if showMag and nRows > 2:
            mags = self.LCsample.refMag - \
                2.5*np.log10(ySampl/self.LCsample.refFlux)  
            eMag = 1.086 * eSampl/ySampl

            axM = fig1.add_subplot(nRows, 1, 2, sharex=ax1)
            dumScatMag = axM.errorbar(\
                tSampl[bSho], mags[bSho], \
                    yerr=eMag[bSho], ls='none', \
                    ms=1, alpha=0.5, marker='o', zorder=2, \
                    color=colorSho, ecolor=colorSho)

            dumOutMag = axM.errorbar(\
                tSampl[~bSho], mags[~bSho], \
                    yerr=eMag[~bSho], ls='none', \
                    ms=1, alpha=0.5, marker='o', zorder=2, \
                    color=colorOut, ecolor=colorOut)

            axM.set_xlabel('MJD (d)')
            axM.set_ylabel('Mag')

        # save the figure to disk
        fig1.savefig(figname, rasterized=False)

    def wrapLoadLightcurves(self, Debug=False):

        """Convenience-wrapper to load lightcurves for simulations and
        convert magnitude to flux"""

        # Customization arguments to be added!

        # template lightcurve (to extract the stddev for the
        # simulation)
        self.lcTemplateFromFile(refMag=self.refMag, refFlux=self.refFlux)

        if Debug:
            self.getTemplateStats()

            # for the moment, be verbose about the magnitude conversion
            print("FakeLC.wrapPrepSims INFO - converting template mag to flux")
            print("FakeLC.wrapPrepSims INFO - ref flux %.2f, ref mag %.2f" \
                      % (self.LCtemplate.refFlux, self.LCtemplate.refMag))

        self.LCtemplate = self.magToFlux(self.LCtemplate)
        self.getTemplateStats()

        # output sampling including gaps
        self.lcObsFromFile(refMag=self.refMag, refFlux=self.refFlux)
        self.LCobs = self.magToFlux(self.LCobs)

    def wrapPrepRednoise(self):

        self.findObsChunks()
        self.buildOvertimes()
    
        # pick red noise factor
        self.pickRNLfactor()

    def wrapPrepSims(self):

        self.wrapLoadLightcurves()
        self.wrapPrepRednoise()

class TrialSet(object):

    """Convenience objet to hold a series of trials"""

    def __init__(self, nTrials=1, FakeTrial=None, \
                     filStats='tmp_trialStats.fits', \
                     nReport = 5, gapMin=999, \
                     doDetrend=False, detDeg=0, \
                     actOnFlux=True, \
                     choiceFom='simpleStd'):

        self.nTrials = nTrials
        
        # array of trial statistics
        self.aStats = np.array([])
        self.aStatsAll = np.array([])

        # output file for statistics (fits is convenient because we
        # can attach a header)
        self.filStatsAll = filStats[:]
        lSplit = os.path.splitext(self.filStatsAll)
        self.filStatsChunks = '%s_chunks%s' % (lSplit[0], lSplit[1])
        
        # the object holding the fake trial and methods
        self.FakeTrial = FakeTrial

        # Which method in class FoM() do we use?
        self.choiceFom = choiceFom[:]

        # report to screen every nReport trials
        self.nReport = nReport

        # min gap for chunking each fake set
        self.gapMin = gapMin
        self.gapAll = 999 # default gap for the entire dataset

        # Ensure the fake lightcurve is in flux units before
        # performing the figure of merit on it?
        self.actOnFlux = actOnFlux

        # Detrending information
        self.doDetrend = doDetrend
        self.detDeg = detDeg

        # *name* of the method used as figure of merit
        self.nameFomUsed = 'BLANK' # easily-recognizable default

        # time range covered by the requested sampling
        self.timeRangeSample = 0.

        # header for output file
        self.hdr = []

    def findTimeRange(self):

        """Finds the time range covered by the requested sampling"""

        tSampl = self.FakeTrial.LCblank.time
        self.timeRangeSample = np.max(tSampl) - np.min(tSampl)

    def doTrials(self):

        """Does the trials"""
        
        # initialize the outputs
        self.aStatsAll = np.array([])
        self.aStats = np.array([])

        # We only need to do the chunks if the min gap is smaller than
        # the time range covered by the data.
        self.findTimeRange()
        doChunks = self.timeRangeSample > self.gapMin

        for iTrial in range(self.nTrials):
            if iTrial % self.nReport < 1 and iTrial > 0:
                sys.stdout.write("\r TrialSet.doTrials INFO - trial %i of %i" \
                                     % (iTrial, self.nTrials))
                sys.stdout.flush()

            self.FakeTrial.sampleNoiseModel()
                
            if iTrial < 1:
                print("TrialSet.doTrials INFO - plotting %i..." \
                              % (iTrial))
                self.FakeTrial.showLC()
                print("... done.")

            # Convert the fake lightcurve to flux before acting on
            # it. Note that the magToFlux and fluxToMag methods do
            # nothing if the inputs are already in the desired output
            # unit.
            if self.actOnFlux:
                self.FakeTrial.LCsample = self.FakeTrial.magToFlux(\
                    self.FakeTrial.LCsample)
            else:
                self.FakeTrial.LCsample = self.FakeTrial.fluxToMag(\
                    self.FakeTrial.LCsample)
            

            FSall = FoMSet(self.FakeTrial.LCsample, \
                               self.FakeTrial.bObs, \
                               gapMin=self.gapAll, \
                               doDetrend=self.doDetrend, \
                               detDeg=self.detDeg, \
                               choiceFom=self.choiceFom, \
                               runOnInit=True)
            self.aStatsAll = self.accumArrays(self.aStatsAll, FSall.aFoms)

            if doChunks:
                FS = FoMSet(self.FakeTrial.LCsample, \
                                self.FakeTrial.bObs, \
                                gapMin=self.gapMin, \
                                doDetrend=self.doDetrend, \
                                detDeg=self.detDeg, \
                                choiceFom=self.choiceFom, \
                                runOnInit=True)
                self.aStats = self.accumArrays(self.aStats, FS.aFoms)
                
        # get the name of the method used for the fom (we might
        # normally specify this here too, but this way we read even if
        # the default was used)
        #self.nameFomUsed = FS.lFoms[0].methFom.__name__
        self.nameFomUsed = FSall.nameFoM
        self.doDetrend = FSall.lFoms[0].detrend
        self.detDeg = FSall.lFoms[0].detDeg

    def parsToHeader(self):

        """Constructs FITS header from the simulation parameters"""

        self.hdr = fits.Header()
        self.hdr['nTrials'] = (self.nTrials, \
                                   "Number of trial sets run")
        self.hdr['fomUsed'] = (self.nameFomUsed, \
                                   "Method used to calculate array contents")
        self.hdr['gapMin'] = ( self.gapAll, \
                                   "Minimum gap between time chunks in output")

        # Now for some more detailed parameters on the
        # simulation. 
        self.hdr['rnMeth'] = ( self.FakeTrial.methPSD.__name__, \
                                   "PSD simulation method")
        self.hdr['rnStd'] = ( self.FakeTrial.sampleStd, \
                                  "PSD simulation std dev")
        self.hdr['rnMean'] = ( self.FakeTrial.sampleMean, \
                                   "PSD simulation mean")
        self.hdr['rnlFac'] = ( self.FakeTrial.rnlFactor, \
                                   "DELCgen RNLfactor")
        self.hdr['tBin'] = ( self.FakeTrial.LCblank.tbin, \
                                 "DELGgen TBIN factor")

        # I don't think fits headers can handle array arguments, so we
        # parcel the PSD params into scalars
        psdPars = self.FakeTrial.PSDpars
        sStem = 'rnPar_'
        cStem = 'PSD parameter'
        for iPar in range(np.size(self.FakeTrial.PSDpars)):
            keyName = '%s%i' % (sStem, iPar)
            commen = '%s %i' % (cStem, iPar)
            self.hdr[keyName] = (psdPars[iPar], \
                                     commen)

        # some simulation control variables
        self.hdr['pertUnc'] = (self.FakeTrial.pertByMeasureUncty, \
                                   "Simulation perturbed by obs uncertainty")
        self.hdr['filSampl'] = (self.FakeTrial.filSamples, \
                                    "Obsn file for sample and/or flux, uncty")
        self.hdr['filTempl'] = (self.FakeTrial.filTemplate, \
                                    "Template lc for simulation")
        self.hdr['uObsFl'] = (self.FakeTrial.useObsFlux, \
                                  "Obsn used for flux")
        self.hdr['uObsErr'] = (self.FakeTrial.useObsUncty, \
                                   "Obsn used for uncertainties")

        # just so that we can trace the provenence of the mags and
        # fluxes:
        self.hdr['isFlux'] = (self.FakeTrial.LCsample.isFlux, \
                                  "Statistic evaluated on flux lightcurve" )
        self.hdr['refMag'] = (self.FakeTrial.LCsample.refMag, \
                                  "Reference apparent mag")
        self.hdr['refFlux'] = (self.FakeTrial.LCsample.refFlux, \
                                   "Reference flux")

        # Was detrending done?
        self.hdr['detrend'] = (self.doDetrend, \
                                   "Samples detrended before evaluating FoM")
        if self.doDetrend:
            self.hdr['detDeg'] = (self.detDeg, \
                                      "Detrending factor for samples")

    def writeTrials(self):

        """Writes the trialset to disk"""

        if len(self.filStatsAll) < 2:
            return

        # make a more sophisticated header
        self.parsToHeader()

        # populate the header with some keyword arguments
#        hdr = fits.Header()
#        hdr['nTrials'] = self.nTrials
#        hdr['fomUsed'] = self.nameFomUsed
#        hdr['gapMin'] = self.gapAll

        # could write various values of the FakeTrial object here too...
        if os.access(self.filStatsAll, os.W_OK):
            os.remove(self.filStatsAll)
            
        fits.writeto(self.filStatsAll, self.aStatsAll, header=self.hdr)

        # write the separate chunks stats file only if we populated it
        if np.size(self.aStats) < 1:
            return

        self.hdr['gapMin'] = self.gapMin
        if os.access(self.filStatsChunks, os.W_OK):
            os.remove(self.filStatsChunks)

        fits.writeto(self.filStatsChunks, self.aStats, header=self.hdr)

        

    def accumArrays(self, aMaster = np.array([]), aNew = np.array([])):

        """Convenience-method to copy or vstack aNew onto aMaster as
        appropriate"""

        if np.size(aMaster) < 1:
            aMaster = np.copy(aNew)
        else:
            aMaster = np.vstack(( aMaster, aNew ))

        return aMaster

class TimeChunks(object):

    """Convenience object to identify chunks in a time-series"""

    def __init__(self, times=np.array([]), chunkMin=0.01, \
                     runOnInit=True):

        self.chunkMin=chunkMin
        self.times = np.copy(times)

        # the chunk ID for each datapoint
        self.chunkIDs = np.array([])

        # the list of unique chunk IDs
        self.uniqueIDs = np.array([])

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

            # we grow the unique IDs array even though we probably
            # could just build it once and be done with it...
            self.uniqueIDs = np.hstack(( self.uniqueIDs, iChunk ))

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
     
class FoMSet(object):

    """Properties and methods for determining the figure of merit on
    an input simulation, including the ability to break up the dataset
    into nights. 

    """

    def __init__(self, LCobj=None, bObs=np.array([]), gapMin=999., \
                     parseOnInit=True, runOnInit=True, \
                     doDetrend=False, detDeg=0, \
                     choiceFom = 'simpleStd'):

        self.LCobj = LCobj
        self.bObs = np.copy(bObs)
        
        # make the time, flux, uncty into separate arrays (from the
        # LCobj) so that we can populate them directly from a calling
        # array if necessary.

        # First the full set...
        self.tFull = np.array([])
        self.yFull = np.array([])
        self.eFull = np.array([])
        
        # then the subset that will pass the bObs boolean
        self.tObs = np.array([])
        self.yObs = np.array([])
        self.eObs = np.array([])

        # Data Chunks object, minimum time-gap between chunks
        self.chunks = None
        self.gapMin = gapMin

        # list of figure of merit objects, array of fom values
        self.aFoms = np.array([])
        self.lFoms = []

        # control variable
        self.Verbose=True

        # record of the FoM name used (watchout - confusing namespace?)
        self.nameFoM=choiceFom[:]

        # detrending information for FoM
        self.doDetrend = doDetrend
        self.detDeg = detDeg

        # setup operations based on input arguments
        if parseOnInit:
            self.unpackLCobj()
            self.checkObsBoolean()
            self.extractObservedData()
            self.partitionIntoChunks()
            self.initFomHolders()

        if runOnInit:
            self.evaluateFoMs()

    def unpackLCobj(self):

        """Unpacks the LC object into time, flux, error arrays"""

        try:
            self.tFull = np.copy(self.LCobj.time)
            self.yFull = np.copy(self.LCobj.flux)
            self.eFull = np.copy(self.LCobj.errors)
        except:
            if self.Verbose:
                print("FomSet.unpackLCobj WARN - problem unpacking LC object")
            return

    def checkObsBoolean(self):

        """Ensure that the observation boolean is compatible with the
        dataset"""
        
        # I'm on the fence about whether this conditional is
        # necessary. If the tObs has not yet been set, then we can
        # easily generate the bObs to match it and deal with the
        # zero-size data elsewhere.
        if np.size(self.tFull) < 1:
            if self.Verbose:
                print("FomSet.checkObsBoolean WARN - times array not populated")
            return
        
        # This one, however, is the key function of this method.
        if np.size(self.bObs) <> np.size(self.tFull):
            if self.Verbose:
                print("FomSet.checkObsBoolean WARN - boolean size mismatch to times.")
                print("FomSet.checkObsBoolean WARN - will use all input datapoints.")
            self.bObs = np.repeat(True, np.size(self.tFull))

    def extractObservedData(self):

        """Selects only the 'observed' data using supplied boolean"""

        self.tObs = self.tFull[self.bObs]
        self.yObs = self.yFull[self.bObs]
        self.eObs = self.eFull[self.bObs]

    def partitionIntoChunks(self):

        """Partitions the dataset into chunks"""
        
        self.chunks = TimeChunks(self.tObs, self.gapMin, runOnInit=True)

    def initFomHolders(self):

        """Sets up a list of figure of merit objects, and an array of
        statistics"""

        nChunks = np.size(self.chunks.uniqueIDs)
        self.aFoms = np.repeat(-99., nChunks)
        self.lFoms = [None for i in range(nChunks)]

    def evaluateFoMs(self):

        """Evaluates the figure of merit for each chunk in turn."""

        # WATCHOUT for future index debugging - notice that this requires
        # the tObs, yObs, eObs arrays, NOT the tFull, yFull, eFull
        # arrays.

        for iChunk in range(np.size(self.chunks.uniqueIDs)):
            thisID = self.chunks.uniqueIDs[iChunk]
            bChunk = self.chunks.chunkIDs == thisID

            # we'll do this in pieces to clarify the operation of the
            # FoM object
            thisFom = FoM(self.tObs[bChunk], \
                              self.yObs[bChunk], \
                              self.eObs[bChunk], \
                              detrend=self.doDetrend, \
                              detDeg = self.detDeg, \
                              runOnInit=False, \
                              choiceFom = self.nameFoM)
            thisFom.chunkID = thisID
            thisFom.calcFoM()

            # 2019-06-20 - debug statement to fix a method-location problem
            # print("DEBUG:", iChunk, thisID, np.sum(bChunk), thisFom.fomStat)

            # now we update the list of FoM objects with this FoM and
            # pass the statistic to the holding array
            self.lFoms[iChunk] = thisFom
            self.aFoms[iChunk] = thisFom.fomStat

        # record the name of the FoM used
        self.nameFoM = self.lFoms[0].methFom.__name__

class FoM(object):

    """Figure of merit object for a single data chunk"""
                      
    def __init__(self, tObs=np.array([]), yObs=np.array([]), \
                     eObs=np.array([]), choiceFom='simpleStd', \
                     runOnInit=False, chunkID=0, \
                     detrend=False, detDeg=0):

        """Figure of merit object for a single dataset, optionally
        with measurement uncertainties"""

        # for convenience, we allow an identifier for which chunk this is
        self.chunkID = chunkID

        # Make copies of the input arrays to avoid circular references
        # or accidentally altering the input arrays themselves.
        self.tObs = np.copy(tObs)
        self.yObs = np.copy(yObs)
        self.eObs = np.copy(eObs)
        
        self.fomStat = -99.

        # doing the detrending?
        self.detrend = detrend
        self.detDeg = detDeg

        # choice of FoM (check that it actually works on np arrays)
        self.choiceFom = choiceFom
        self.setMethFom()  

        # verbose?
        self.Verbose = True

        if runOnInit:
            self.detrendY()
            self.setMethFom()
            self.calcFom()

    def setMethFom(self):

        """Sets the figure of merit to calculate"""

        choiceFom = self.choiceFom[:]
        choiceDefault = 'simpleStd'
        if not hasattr(self, self.choiceFom):
            if self.Verbose:
                print("FoM.checkMethFom WARN - requested FoM method not found")
                print("FoM.checkMethFom WARN - defaulting to %s" \
                          % (choiceDefault))
            choiceFom = choiceDefault[:]

        # set the method and record the name
        self.methFom = getattr(self, choiceFom)

    def detrendY(self):

        """Do the detrending of the datapoints"""

        if not self.detrend:
            return

        wts = np.ones(np.size(self.yObs))
        if np.size(self.eObs) <> np.size(self.yObs):
            wts = 1.0/self.eObs**2
                      
        detPars = np.polyfit(self.tObs, self.yObs, self.detDeg, \
                                 w=wts)

        yTrend = np.polyval(detPars, self.tObs)
        self.yObs = self.yObs - yTrend

    def calcFoM(self):

        """Calculates the figure of merit on the input dataset.
        """

        self.methFom()
        
        # 2019-06-20 WIC - I'm wondering if this even needs a method,
        # since we already set the reference in setMethFom!

    def simpleStd(self):

        """Computes the standard deviation of the input data"""

        self.fomStat = np.std(self.yObs)

    def quadDiff(self):

        """Find the quadrature difference between the stddev and the
        mean uncertainty"""

        varDiff = np.std(self.yObs)**2 - np.median(self.eObs)**2
        self.fomStat = np.sqrt(varDiff)
        

class Domino(object):

    """Straight resampling of v404 lightcurve using historical
    data."""

    def __init__(self, LCobj=None, doShuffleTiles=False, \
                     tZero=50000., \
                     doFlipToMatch=False, \
                     padFactor=0, \
                     runOnInit=True):

        """Accepts a lightcurve object, breaks into chunks, and offers
        various resampling options."""

        self.LC = LCobj  # required input
        self.doShuffleTiles = doShuffleTiles # control variable
        self.doFlipToMatch = doFlipToMatch # flip tiles to match?
        self.tZero = np.float(tZero) # time value of first chain
                                     # datapoint
        self.padFactor = padFactor
        self.rollStep = 2 # for padding. Experiment with this to best
                          # distribute the source tiles most
                          # evenly. In principle, 1 or 2 should work
                          # well.

        # some internal variables
        self.objTiles = None
        self.tileLen = 0.7  # days, in most applications will take
                            # this value

        # housekeeping variables for the original tile IDs and times
        self.origIDs = np.array([])
        self.origOrder = np.array([])
        self.origTimes = np.array([])  # might be redundant

        self.chainIDs = np.array([])  # which tile each obs belongs to
        self.chainTimes = np.array([])
        self.chainFLux = np.array([])
        self.chainUncty = np.array([])

        # DELCgen "Lightcurve" object populated by the chain
        self.chainLC = None

        # set up the chain object in preparation for construction
        self.findTiles()
        self.initChainIDs()

        if runOnInit:
            self.buildChain()
            #if self.doShuffleTiles:
            #    self.shuffleTiles()
            #self.replicateChain()
            #self.populateChainValues()
            #self.chainLCfromArrays()

    def findTiles(self):

        """Breaks the time-series into nights"""

        # (called 'tiles' to avoid namespace collision with the chunks
        # in other objects)
        self.objTiles = TimeChunks(self.LC.time, self.tileLen, \
                                       runOnInit=True)
        
        # set some useful instance-level properties
        self.origIDs = np.copy(self.objTiles.chunkIDs)
        self.origOrder = np.asarray(\
            np.unique(np.sort(self.objTiles.chunkIDs)), 'int')
        self.origTimes = np.copy(self.objTiles.times)

    def initChainIDs(self):

        """Initialises chain IDs from the original IDs"""

        self.chainOrder = np.copy(self.origOrder)

        # (This was originally longer since it did more...)

    def shuffleTiles(self):

        """Randomly reorders the time tiles in the chain"""

        xDum = np.random.uniform(size=np.size(self.origOrder))
        lDum = np.argsort(xDum)
        self.chainOrder = self.origOrder[lDum]

    def replicateChain(self):

        """Replicates the chain by padding factor self.nPad, rolling
        the chain by 1 for each replication. Useful in the case when
        the observations only cover part of a 24h cycle."""

        if self.padFactor < 1:
            return

        if len(self.chainOrder) < 1:
            return

        firstChain = np.copy(self.chainOrder)
        for iRep in range(self.padFactor-1):
            chainToAdd = np.roll(firstChain, 0 - (self.rollStep*iRep+1))
            self.chainOrder = np.hstack(( self.chainOrder, chainToAdd ))


    def initChainValues(self):

        """Initializes the chain times, flux, uncties"""

        self.chainTimes = np.array([])
        self.chainFlux = np.array([])
        self.chainUncty = np.array([])
        self.chainIDs = np.array([])

    def populateChainValues(self):

        """Populates the time, flux, uncty values for the chain of
        tiles"""

        tMax = np.copy(self.tZero)
        self.initChainValues()

        # we draw from the original set by tile ID, with the ordering
        # dictated by the reordered set.
        for iTile in range(len(self.chainOrder)):
            bThis = self.origIDs == self.chainOrder[iTile]

            if np.sum(bThis) < 1:
                continue

            # Since we're abutting samples together, we need to
            # re-compute the times. We subtract off the smallest time
            # in this tile, and add the median time interval between
            # datapoints in THIS tile. 
            tThis = self.origTimes[bThis] - np.min(self.origTimes[bThis])
            tSort = np.sort(tThis)
            tThis += np.median(np.abs(tSort - np.roll(tSort, -1)))
            tThis += tMax

            tMax = np.max(tThis)

            # To help the domino tiles match -- if doFlipToMatch is
            # set -- find which end of the domino best matches the
            # tail of the seequence, and abut the closest ends onto
            # each other.
            if self.doFlipToMatch and np.size(self.chainFlux) > 1:
                avgTail = np.mean(self.chainFlux[-6:-1])
                avgFron = np.mean(self.LC.flux[bThis][0:6])
                avgBack = np.mean(self.LC.flux[bThis][-6::])
                if np.abs(avgFron - avgTail) > np.abs(avgBack - avgTail):
                    tThis = tThis[::-1]

            # we now have our time, flux and uncertainties to slot in.
            self.chainTimes = np.hstack(( self.chainTimes, tThis ))
            self.chainFlux = np.hstack(( self.chainFlux, \
                                             self.LC.flux[bThis] ))
            self.chainUncty = np.hstack(( self.chainUncty, \
                                              self.LC.errors[bThis] ))

            theseIDs = np.repeat(self.chainOrder[iTile], np.size(tThis))
            self.chainIDs = np.hstack(( self.chainIDs, theseIDs ))

    def chainLCfromArrays(self):

        """Populates a DELCgen Lightcurve object with the chain
        values"""

        # note we do NOT want to simply copy the original tile because
        # the chain will often have a different number of datapoints.
        if not hasattr(self.chainLC, 'time'):
            tBin = np.copy(self.LC.tbin)
            self.chainLC = DELCgen.Lightcurve( self.chainTimes, \
                                                   self.chainFlux, \
                                                   tBin, \
                                                   self.chainUncty)

            # now we copy across the extra attributes we added
            self.chainLC.isFlux = np.copy(self.LC.isFlux)
            self.chainLC.refFlux = np.copy(self.LC.refFlux)
            self.chainLC.refMag = np.copy(self.LC.refMag)

        else:
            self.chainLC.time = self.chainTimes
            self.chainLC.flux = self.chainFlux
            self.chainLC.errors = self.chainUncty

            # If not initializing, we need to recalculate the
            # statistics
            recomputeLCstats(self.chainLC)

        # we will overload the Lightcurve object yet again, this time
        # with the chainIDs (which will come in handy when debugging
        # to learn from which tile each output point was drawn)
        self.chainLC.chainIDs = np.copy(self.chainIDs)

    def buildChain(self):

        """Wrapper to create a new chain once the input data has been
        appropriately set up"""

        if self.doShuffleTiles:
            self.shuffleTiles()
        self.replicateChain()
        self.populateChainValues()
        self.chainLCfromArrays()

class Counterfeiter(object):

    """Makes an imperfect copy of a lightcurve by sampling at
    particular points. Uses the DELCgen.Lightcurve object for
    convenience"""

    def __init__(self, LCtemplate=None, LCobs=None, \
                     interpKind='nearest'):

        """Initialize"""
    
        # the template and "observation" lightcurve objects
        self.LCtemplate = LCtemplate
        self.LCobs = LCobs

        # Lightcurve containing the sample
        self.LCsample = None

        ## Housekeeping attributes
        self.tDiffsObs = np.array([]) # deltas from start point
        self.tStartMin = -99.
        self.tStartMax = -99.

        # the interpolation objects
        self.interpFlux = None
        self.interpUncty = None
        self.interpID = None # useful to trace tile provenance
        self.interpKind = interpKind[:] # trust the user to know what
                                        # this is

        # should we perturb the interpolate by the uncertainty after
        # drawing it? (Not always warranted, e.g. if we're selecting
        # the nearest point only)
        self.doApplyUncty = False

        # now for the time samples
        self.sampleTimes = np.array([])
        self.sampleFlux = np.array([])
        self.sampleUncty = np.array([])
        self.sampleIDs = np.array([])

        # Common initialization tasks
        self.findObsTimediffs()
        self.setTimeBounds()
        self.setupInterpolation()

    def setTimeBounds(self):

        """The time lengths of the template and observation
        lightcurves set limits on how close to the edge we can set our
        start and stop times for the sample lightcurve."""
        
        obsDuration = np.max(self.LCobs.time) - np.min(self.LCobs.time)
        self.tStartMin = np.min(self.LCtemplate.time)
        self.tStartMax = np.max(self.LCtemplate.time) - obsDuration

        # build a half-hour buffer into the max start
        self.tStartMax -= 0.5/1440.

    def findObsTimediffs(self):

        """Turns all the observation times into time-diffs from the
        first datapoint"""

        self.tDiffsObs = self.LCobs.time - np.min(self.LCobs.time)

    def initSampleLC(self):

        """Initialises the sampling lightcurve"""

        # fortunately, this is easy!
        self.LCsample = copy.deepcopy(self.LCobs)

    def plantSampleTimes(self, t0=-99.):

        """Places a set of sample times using the observational time
        differences but starting at t0. If t0 unset, generates a new
        one from the bounds"""

        # I envisage this being called from a method that has already
        # generated a set of random start times. However, let's give
        # this the option to set its own (this is also a good time to
        # apply quality control to the input start time).
        if t0 < 0 or \
                t0 < self.tStartMin or \
                t0 > self.tStartMax:
            t0 = np.random.uniform(self.tStartMin, self.tStartMax)
            
        self.sampleTimes = self.tDiffsObs + t0

    def drawSampleValues(self):

        """Uses the interpolators to draw samples of the lightcurve at
        the given times in the sample."""
        
        self.sampleFlux = self.interpFlux(self.sampleTimes)
        self.sampleUncty = self.interpUncty(self.sampleTimes)
        self.sampleIDs = self.interpID(self.sampleTimes)

        if not self.doApplyUncty:
            return

        dY = np.random.normal(size=np.size(self.sampleUncty)) \
            * self.sampleUncty

        self.sampleFlux += dY

    def passSampleToLC(self):

        """Passes the current sample to the sample lightcurve
        object"""

        if not hasattr(self.LCsample, 'time'):
            self.initSampleLC()

        self.LCsample.time = self.sampleTimes
        self.LCsample.flux = self.sampleFlux
        self.LCsample.errors = self.sampleUncty

        # we recompute the statistics...
        recomputeLCstats(self.LCsample)

        #... and overload with our interpolated ID value
        self.LCsample.chainID = self.sampleIDs

    def makeSample(self, t0=-99):

        """Wrapper - plants sample times, draws sample values and
        passes the results up to the sample lightcurve object"""

        self.plantSampleTimes(t0)
        self.drawSampleValues()
        self.passSampleToLC()

    def setupInterpolation(self):

        """Sets up the interpolation objects for flux and uncertainty
        once the template lightcurve is populated """

        self.interpFlux = interpolate.interp1d(\
            self.LCtemplate.time, self.LCtemplate.flux, \
                kind=self.interpKind, \
                copy=False)

        self.interpUncty = interpolate.interp1d(\
            self.LCtemplate.time, self.LCtemplate.errors, \
                kind=self.interpKind, \
                copy=False)

        # Also useful to interpolate the original ID (so that we can
        # learn what coverage we are achieving.
        if not hasattr(self.LCtemplate, 'chainIDs'):
            return

        self.interpID = interpolate.interp1d(\
            self.LCtemplate.time, self.LCtemplate.chainIDs, \
                kind='nearest', \
                copy=False)

def recomputeLCstats(LC=None):

    """General method to re-perform the statistics for a DELCgen
    lightcurve. Modifies the Lightcurve object in-place."""

    LC.mean = np.mean(LC.flux)
    LC.std = np.std(LC.flux)
    LC.length = len(LC.flux)
    LC.freq = np.arange(1, LC.length/2.0 + 1)/(LC.length*LC.tbin)

    # this --should-- modify in place. Let's see if it does.

def binData(tIn=np.array([]), yIn=np.array([]), \
                nBins=100, nMin=2, \
                eIn=np.array([]), \
                useMedian=False):
#, \
#                logBins=False):

    """Bin time, rate, optionally uncertainty. Slow but robust. Uses
    nBins to specify the binning so that log-spaced binning is
    straightforward to implement if needed.

    useMedian = use the median (and std for uncertainty)

    """

    # lower and upper ranges for bins
    tMin = np.min(tIn)
    tMax = np.max(tIn)

    xAll = np.linspace(tMin, tMax, nBins+1, endpoint=True)
    #if logBins:
    #    xAll = np.logspace(np.log10(tMin), np.log10(tMax), nBins, \
    #                           endpoint=True)
    tLos = xAll[0:-1]
    tHis = xAll[1::]

    # weights if input errors set
    wgts = yIn*0.0 + 1.0
    if np.size(eIn) > 0:
        wgts = 1.0/eIn**2

    # Initialize the output bins
    tBin = np.zeros(nBins)
    yBin = np.copy(tBin)
    eBin = np.copy(tBin)
    bBin = np.repeat(False, nBins)

    for iBin in range(nBins):
        
        #tMin = tLos[iBin]
        #tMax = tHis[iBin]

        bThis = (tIn >= tLos[iBin]) & (tIn < tHis[iBin])
        nThis = np.sum(bThis)
        if nThis < nMin:
            continue

        bBin[iBin] = True

        # If uncertainties are given, use inverse-variance weights and
        # find the optimal average and its variance. Otherwise, use
        # the straight average and stddev. Only the uncertainty
        # differs from the two cases (if errors not given, every point
        # is given equal weight).


        if useMedian:
            tBin[iBin] = np.median(tIn[bThis])
            yBin[iBin] = np.median(yIn[bThis])
        else:
            sumWgts = np.sum(wgts[bThis])
            tBin[iBin] = np.sum(tIn[bThis]*wgts[bThis])/sumWgts
            yBin[iBin] = np.sum(yIn[bThis]*wgts[bThis])/sumWgts

        # if no errors provided OR using the median, return the stddev
        if np.size(eIn) <> np.size(tIn) or useMedian:
            eBin[iBin] = np.std(yIn[bThis])
        else:
            # if we have uncertainties, use the variance of the
            # inverse-variance-weighted average.
            eBin[iBin] = np.sqrt(1.0/sumWgts)

    # only return the complete bins
    return tBin[bBin], yBin[bBin], eBin[bBin]

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


def testCombinedSample(nTrials=1, gapMin=0.7):

    """Tests building a combined sample with observations and a wider
    time baseline to distinguish observations from the wider
    sample.

    nTrials = number of trials to run

    gapMin = minimum gap between "nights" in the observed dataset. Set
    to a large number to use all the data in a single chunk.
    """

    FLC = FakeLC('DIA2017.csv')

    # Use the previous file (if readable) for the broad statistics for
    # the simulation
    FLC.lcTemplateFromFile()
    FLC.getTemplateStats()

    # set up the sampling, including the gaps
    FLC.lcObsFromFile()
    FLC.findObsChunks()
    FLC.buildOvertimes()

    # since we now have a larger sample, we don't need such a large
    # RNLfactor. Use a smaller range.
    #FLC.rnlMin=10
    #FLC.rnlMax=20
    FLC.pickRNLfactor()
    
    ### 2019-06-20 all the material below could be refactored into a
    ### new class that just holds the simulation set... what do we
    ### think?

    # At this point we're ready to run our simulation set. Now we
    # generate an array for holding the big set of figures of merit

    # reporting interval (since I'm an impatient user)
    nReport = 5

    aFoms = np.array([])
    for iTrial in range(nTrials):

        # report progress to screen, since the simulations might be
        # slow-ish...
        if iTrial % nReport < 1 and iTrial > 0:
            sys.stdout.write("\r testCombinedSample INFO - on trial %i of %i" \
                                 % (iTrial, nTrials))
            sys.stdout.flush()

        # Sample the noise model...
        FLC.sampleNoiseModel()
    
        # ... and show the lightcurve (just do this once)
        if iTrial < 1:
            print("testCombinedSample INFO - plotting trial %i..." \
                      % (iTrial))
            FLC.showLC()
            print("testCombinedSample INFO - done.")

        # test the FomSet object for the figures of merit
        FS = FoMSet(FLC.LCsample, FLC.bObs, gapMin=gapMin)
        FS.evaluateFoMs()
    
        # if we haven't built the fom array yet, copy it from the
        # first instance
        if np.size(aFoms) < 1:
            aFoms = np.copy(FS.aFoms)
        else:
            aFoms = np.vstack(( aFoms, FS.aFoms ))

    print np.shape(aFoms)
    print aFoms[0]

def testCompact(nTrials=1, gapMin=0.7, \
                    filObstimes='DIA2017.csv', \
                    filTemplate='92Binned.fits', \
                    doDetrend=False, detDeg=0, \
                    actOnFlux=True, \
                    choiceFom='simpleStd', \
                    testFits = False):

    """Performs the simulation. Arguments:

    nTrials = number of trial sets to generate

    gapMin = minimum gap between chunks in the output data. The
    default 0.7 breaks the output into night-by-night chunks. Note
    that the figure of merit is ALWAYS calculated for the entire
    output dataset; setting gapMin to a value < the total time
    baseline ensures that the FoM is also calculated for the chunks
    (and stored in a separate output file).

    filObstimes = file used for the output observation times

    filTemplate = file to be used for the mean and standard deviation
    when simulating the red noise datasets.

    doDetrend -- detrend the output data before evaluating the figure of merit?

    detDeg = polynomial degree for detrending (if doDetrend is True)

    actOnFlux -- ensure the simulation is in "flux" units before
    calculating the fiture of merit? (Otherwise the simulation will be
    in "mag" units)

    choiceFom = name of the method in class FoM() used to compute the
    "figure of merit" for each simulated dataset

    testFits -- test the use of uncertainties from 2017 binned
    lightcurve? (Under development)

    --
    
    Example call to generate 4 trials, with constant-level
    subtraction, transforming the simulation from flux to apparent
    magnitude and using the error-corrected standard deviation as a
    figure of merit:

    v404_sample.testCompact(4, doDetrend=True, detDeg=0, actOnFlux=False, choiceFom='quadDiff')

    """

    FLC = FakeLC(filSamples=filObstimes, filTemplate=filTemplate)

    # adjust the keywords for AMB's binned lightcurves
    if testFits:
        FLC.filSamples = '17Binned.fits'
        FLC.keyObsTime = 'tBin'
        FLC.tSamplesFactor = 1.0
        FLC.keyObsUnct = 'uBin'
        FLC.keyObsFlux = 'fBin'

        FLC.useObsUncty=False
        FLC.useObsFlux=False

    # prepare for simulations
    FLC.wrapPrepSims()

    # now use an instance of our TrialSet class to run the iterations
    TS = TrialSet(nTrials, FakeTrial=FLC, gapMin=gapMin, \
                      doDetrend=doDetrend, detDeg=detDeg, \
                      actOnFlux = actOnFlux, \
                      choiceFom = choiceFom[:])
    TS.doTrials()
    TS.writeTrials()


def testFluxConversion(filRef='92Binned.fits'):

    """Tests conversion from magnitude to flux"""

    FF = FakeLC()
    FF.lcTemplateFromFile(refMag = 16.)
    FF.getTemplateStats()

    FF.LCtemplate = FF.magToFlux(FF.LCtemplate)
    FF.getTemplateStats()

    print FF.LCtemplate.isFlux

def testDominoes(filHistoric='92Binned.fits', \
                     filSamples='17Binned.fits', \
                     padFac=2, interpKind='nearest', \
                     doApplyUncty = False):

    """Tester for the direct resampling using 'dominoes' and
    interpolation using the 'Counterfeiter'"""

    # This includes diagnostic plots and timing lines so is somewhat
    # messy.

    # set up our methods for lightcurve object handling, passing
    # arguments appropriate for our 'historic' file as input arguments
    FL = FakeLC(filTemplate=filHistoric, \
                    filSamples=filSamples, \
                    keyObsTime='tBin', \
                    keyObsFlux='fBin', \
                    keyObsUnct='uBin', \
                    tSamplesFactor = 1.0)

    # 2019-07-05 (WIC) - I think this should be the default.
    FL.errDivisorTemplate = 3.
    
    FL.wrapLoadLightcurves()

    # Convert to flux and set up the template 
    #FL.lcTemplateFromFile(isMag=True)
    #FL.LCtemplate = FL.magToFlux(FL.LCtemplate)  # converts in-place
    #FL.getTemplateStats()

    ## load the observations file containing the times (or more
    ## importantly the time-differences) that will correspond to our
    ## samples.
    #FL.lcObsFromFile()
    #FL.LCobs = FL.magToFlux(FL.LCobs)

    # now set up the dominoes object. For testing, do one set with and
    # the other without shuffling.
    DOM = Domino(FL.LCtemplate, runOnInit=True, doShuffleTiles=False)

    tZer = time.time()
    DOMSHUF = Domino(FL.LCtemplate, runOnInit=True, doShuffleTiles=True, \
                         doFlipToMatch=True, padFactor=padFac)

    print("testDominoes INFO - time to set the chain to copy: %.3e s" \
              % (time.time() - tZer))

    # now let's see how long this takes to rebuild a chain after
    # initialization
    tZer = time.time()
    DOMSHUF.buildChain()
    print("testDominoes INFO - time to rebuild chain: %.3e" \
              % (time.time() - tZer))

    tPreSet = time.time()
    # now we try setting up a resampling.
    CF = Counterfeiter(DOMSHUF.chainLC, FL.LCobs, interpKind)

    # I don't trust nudging if we're just finding the
    # nearest. However, for testing purposes, let's give ourselves the
    # chance to toggle it as input arguments
    #CF.doApplyUncty = interpKind.find('nearest') < 0
    CF.doApplyUncty = doApplyUncty

    CF.initSampleLC()

    tZer = time.time()
    print("testDominoes INFO - time to set up Counterfeiter: %.3e s" \
              % (tZer - tPreSet))
    CF.makeSample()
    print("testDominoes INFO - time to draw a single sample: %.3e s" \
              % (time.time() - tZer))

    # try the second time through.
    # CF.makeSample()

    #CF.plantSampleTimes()
    #CF.drawSampleValues()
    #CF.passSampleToLC()

    # debug - print the list of tile-orderings
    ##print DOM.chainOrder
    ##print DOMSHUF.chainOrder

    # build the figure size programmatically
    figSz = (10,4)
    figSz2 = np.copy(figSz)
    figSz2[0] /= 2.
    cmap = plt.cm.prism

    fig1 = plt.figure(1, figsize=figSz)
    fig1.clf()
    ax1 = fig1.add_subplot(121)
    dum1 = ax1.errorbar(DOM.chainTimes, DOM.chainFlux, \
                            DOM.chainUncty, \
                            ms=.1, zorder=3, ecolor='0.6', \
                            ls='none')
    dum12 = ax1.scatter(DOM.chainTimes, DOM.chainFlux, \
                            c=DOM.chainIDs, zorder=10, \
                            cmap=cmap, s=9, \
                            edgecolor='0.5')

    ax2 = fig1.add_subplot(122)
    dum2 = ax2.errorbar(DOMSHUF.chainTimes, DOMSHUF.chainFlux, \
                            DOMSHUF.chainUncty, \
                            ms=.1, zorder=3, ecolor='0.6', \
                            ls='none', alpha=0.3)
    dum22 = ax2.scatter(DOMSHUF.chainTimes, DOMSHUF.chainFlux, \
                            c=DOMSHUF.chainIDs, zorder=10, \
                            cmap=cmap, s=9, \
                            edgecolor='0.5', alpha=0.3)


    # if we have interpolates, draw them!
    if np.size(CF.sampleUncty) < 1:
        return

    dumOver = ax2.errorbar(CF.sampleTimes, CF.sampleFlux, \
                               CF.sampleUncty, \
                               color='0.1', alpha=0.7, \
                               zorder=15, \
                               linestyle='None', \
                               ecolor='0.5')
    dumOverScatt = ax2.scatter(CF.sampleTimes, CF.sampleFlux, \
                                   c='k', edgecolor='0.5', \
                                   s=16, zorder=20)

    # now let's make another figure, this time showing the sampled
    # lightcurve ONLY.
    fig2 = plt.figure(2, figsize=figSz)
    fig2.clf()
    ax3 = fig2.add_subplot(121)
    dumSample = ax3.errorbar(CF.sampleTimes, CF.sampleFlux, \
                               CF.sampleUncty, \
                               color='0.5', alpha=0.5, \
                               zorder=15, \
                               linestyle='None', \
                               ecolor='0.5')
    dumScatt = ax3.scatter(CF.sampleTimes, CF.sampleFlux, \
                               c=CF.sampleIDs, edgecolor='0.5', \
                               s=36, zorder=20, cmap=cmap)

    ax4 = fig2.add_subplot(122)
    dumHist = ax4.hist(CF.sampleIDs, range=(0,20))

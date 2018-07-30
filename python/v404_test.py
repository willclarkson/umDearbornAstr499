# All imports from terminal first

import matplotlib.pylab as plt
plt.ion()
plt.style.use('ggplot')
from astropy.table import Table
import numpy as np
from scipy import optimize
import sys
from astropy.stats import sigma_clip
import loadOld #v404_test must be located in the same directory as loadOld for this to work
from astroML.time_series import lomb_scargle, lomb_scargle_bootstrap

# for testing existence of files
import os

# Reading the table

try:
	t92 = loadOld.loadPhot(fileIn="../../../../Shared/Data/v404Cyg_zurita04/199208elip_R.dat")
	t98 = loadOld.loadPhot(fileIn="../../../../Shared/Data/v404Cyg_zurita04/elipR1998.dat")
	t99A = loadOld.loadPhot(fileIn="../../../../Shared/Data/v404Cyg_zurita04/mv4_990706_apcorr.dat")
	t99B = loadOld.loadPhot(fileIn="../../../../Shared/Data/v404Cyg_zurita04/mv4_990707_apcorr.dat")
	t03 = loadOld.loadPhot(fileIn="../../../../Shared/Data/v404Cyg_zurita04/28_July_2003_iac80.qdp")
	t17 = loadOld.loadPhot(fileIn="../../../Desktop/v404_2017.txt")
	t18A = loadOld.loadPhot(fileIn="../../../Desktop/v404_2018A.txt")
except UnboundLocalError:
	print "One or more files is in the wrong directory. Please edit the script to ensure the code references the correct directory."

## t1 = Table.read('2018-03-09_v404CygPhotom_flag1.csv', format='ascii.csv') #Assuming the file is run from the Desktop, since that is where the script is saved
## t2 = Table.read('2018-03-09_v404CygPhotom_flag2.csv', format='ascii.csv')
## t3 = Table.read('2018-03-09_v404CygPhotom_flag3.csv', format='ascii.csv')
## t12 = Table.read('onesandtwos.csv', format='ascii.csv')
## t123 = Table.read('allflags.csv', format='ascii.csv')
## t7 = Table.read('onesandtwos_7days.csv', format='ascii.csv')
## t121 = Table.read('ones_a12.csv', format='ascii.csv')
## t122 = Table.read('twos_a12.csv', format='ascii.csv')
#t1212 = Table.read('onesandtwos_a12.csv', format='ascii.csv')

#Test

#useFlag = 12
## if useFlag == 1:
## 		tbl = t1
## 		flag = '1'
## if useFlag == 2:
## 		tbl = t2
## 		flag = '2'
## if useFlag == 3:
## 		tbl = t3
## 		flag = '3'
## if useFlag == 12:
## 		tbl = t12
## 		flag = tbl['Flag']
## if useFlag == 123:
## 		tbl = t123
## 		flag = tbl['Flag']

#dy = t1212['rel_flux_err_T1']
#jd = t1212['J.D.-2400000']
#mag = t1212['Flux_V404']
#flag = t1212['Flag']

def go(pctile=10., iCheck=1, useMags=True, \
	       clipOutliers=False, \
	       #oldAperture=False, plotRawCounts=False, \
	       useFlag=12, \
	       tStart=-1e9, \
	       tEnd = 57993., \
	       #tEnd=57992.0, \
	       binTime=0.0034722, bootstrapping=False, \
	       errorbars=True, \
	       magContam=17.5, magCompar=16.07, \
	       showPhase=True, \
	       writeOnly=False, \
	       data='2017', \
	       noCorrection=False, showBinned=False, \
	       plotBinnedData=False, \
	       plotBinnedLS=False, plotLS=False, \
	       plotNoiseData=False, \
	       plotNoiseLS=False, plotHistogram=False, \
	       plotBinnedNoiseData=False, plotBinnedNoiseLS=False, \
	       plotSubtractedData=False, plotEllipsoidal=False):
	# WIC - put the table reading back into go, to avoid scope
	# confusion

	warn = "WARNING: Ellipsoidal modulation will not work for 1999 and 2003!"
	warn2 = "WARNING: Binned and subtracted data only plots if showBinned, plotEllipsoidal, and plotSubtractedData are all set to True."

	if plotBinnedData or plotBinnedLS or plotBinnedNoiseData or plotBinnedNoiseLS or plotSubtractedData:
		print warn2


	if plotBinnedLS or plotBinnedNoiseLS or plotLS:
		print "WARNING: LS plots only work if the data the LS is being performed on is also plotted."

	if data == '1992':
		tbl = t92
	elif data == '1998':
		tbl = t98
		print "WARNING: Due to 1998 data already being in phase, binning and subtraction currently does not work."
	elif data == '1999A':
		tbl = t99A
		print warn
	elif data == '1999B':
		tbl = t99B
		print warn
	elif data == '2003':
		tbl = t03
		print warn
	elif data == '2017':
		tbl = t17
	elif data == '2018A':
		tbl = t18A
		print "WARNING: The current datafile does not contain the last night of data!"
	else:
		print "The value for 'data' must be either 1992, 1998, 1999A, 1999B, 2003, 2017, or 2018A."
		return

	# t1212 = Table.read(inPath, format='ascii.csv')
	# dy = t1212['rel_flux_err_T1']
	# jd = t1212['J.D.-2400000']
	# # mag = t1212['Flux_V404']
	# if 'Flag' in t1212.colnames:
	# 	flag = t1212['Flag']
	# else:
	flag = np.ones(len(tbl), 'int')
	tbl['Flag'] = flag

	# if oldAperture:
	# 	# WARN - this won't work without importing all the
	# 	# tables. Consider replacing the view with the loading
	# 	# from disk.
	# 	if useFlag == 1:
	# 		tbl = t1
	# 		flag = '1'
	# 	if useFlag == 2:
	# 		tbl = t2
	# 		flag = '2'
	# 	if useFlag == 3:
	# 		tbl = t3
	# 	formatlag = '3'
	# 	if useFlag == 12:
	# 		tbl = t12
	# 		flag = tbl['Flag']
	# 	if useFlag == 123:
	# 		tbl = t123
	# 		flag = tbl['Flag']
	# 	if useFlag == 7:
	# 		tbl = t7
	# 		flag = tbl['Flag']		
	#else:
		#if useFlag == 1:
			#tbl = t121
			#flag = '1'
		#if useFlag == 2:
			#tbl = t122
			#flag = '2'
		#if useFlag == 12:
			#tbl = t1212
			#flag = tbl['Flag']
	#if useFlag < 4:
		#flagColor = 'red'
	#else:
	flagColor = flag

	# Trying, but failing, to plot values by flag...
	##if flag == 1:
		##flagColor = 'green'
	##elif flag == 2:
		##flagColor = 'orange'
	##elif flag == 3:
		##flagColor = 'red'

	# 2018-04-07 WIC - went back to the output from AIJ. This is
	# the flux from v404 + contaminant as a multiple of the star
	# "C4" of Casares et al. (1993).
	if 'flux' in tbl.colnames:
		relFlux = tbl['flux']
		if 'fluxErr' in tbl.colnames:
			errFlux = tbl['fluxErr']

		# Since now we know the contaminant and the comparison in
		# apparent magnitude space, we convert flux to magnitude.
		magBoth = magCompar - 2.5*np.log10(relFlux)
		if 'fluxErr' in tbl.colnames:
			errMag = 1.086 * errFlux
	else:
		magBoth = tbl['mag']
		if 'magErr' in tbl.colnames:
			errMag = tbl['magErr']
	
	# so, get the apparent magnitude of v404 cyg itself.
	if not noCorrection:
		mag = correctMagForContaminant(magBoth, magContam)
		print "WARNING: if using data from Z04, it will be INCORRECT unless noCorrection=True!"
	else:
		mag = np.copy(magBoth)
	dy = np.copy(errMag)

	if 'time' in tbl.colnames:
		jd = tbl['time']
		# compute the phase on the orbital ephemeris.
		phase, u_phase = phaseFromJD(jd)
	else:
		phase = tbl['phase']
		jd = tbl['phase'] #so the rest of the script doesn't return an error.........

	# 2018-04-07 useful to plot the RAW data. Rather than
	# interrupt the flow, we port this off to another method.
	
	#if plotRawCounts:
		#showRawCounts(tbl)

	#mag = tbl['Average Mag(V404)']  # superseded by the material above.
	#dy = 1.086 * tbl['rel_flux_err_T1']
	num = len(mag)


	if not useMags:
		mag = np.copy(relFlux)
		dy = np.copy(errFlux)
		#mag = tbl['rel_flux_T1'] 
		#dy = tbl['rel_flux_err_T1']

	### 2018-06-01 WIC - throw in an argument to just plot for one
	### night only
	if writeOnly:
		plt.style.use('ggplot')
		fig2 = plt.figure(2, figsize=(10,4))
		fig2.clf()

		#print jd

		#print np.shape(jd)
		#print np.shape(mag)

		tBin, fBin, uBin, nBin = \
		    BinData(jd, mag, dy, tStart=tStart, tEnd=tEnd, \
				    BinTime=binTime, plotDBG=True)


		xSho = np.copy(jd)
		ySho = np.copy(mag)
		eSho = np.copy(dy)

		if showBinned:
			xSho = np.copy(tBin)
			ySho = np.copy(fBin)
			eSho = np.copy(uBin)

		ax1 = fig2.add_subplot(111)
		#dum = ax1.plot(jd, mag, 'bo')
		dum2 = ax1.errorbar(xSho, ySho, eSho, fmt='bo', ls='None', \
					    ms=1, alpha=0.5)

	#	dum2 = ax1.scatter(tBin, fBin, color='b', alpha=0.5, s=4)


		#dum2 = ax1.errorbar(tBin, fBin, uBin, fmt='b.', ls='None')

		yMin = np.percentile(fBin, 1)

		yRang = np.copy(ax1.get_ylim())
		ax1.set_ylim(yMin+0.4, yMin)
		
		### 2018-06-02 WIC - thrown-together estimate for the
		### chisq
		bLo = xSho < 58270.953695
		fig3 = plt.figure(3, figsize=(5,5))
		fig3.clf()
		ax3 = fig3.add_subplot(111)

		chi = (ySho - np.median(ySho))/eSho

		dum = ax3.hist(chi, 100)
		print "INFO - std(chi) = %.3f" % (np.std(chi))

		# also print out the uncertainties
		print "INFO: dy median, std: %.3e, %.3e" \
		    % (np.median(eSho), np.std(eSho))

		ax1.plot(xSho[bLo], ySho[bLo], 'r.', zorder=5)

		return


	# This was my (failed) attempt to use the lambda function for
	# y like in the example


	# y = lambda p, x: ((p[0]*np.sin(((2*np.pi*x)+p[1])/p[2]))+(p[3]*np.sin(((2*np.pi*x)+p[1]/0.5*p[2]))))+p[4]

	# Initial Guess for Parameters

	if plotEllipsoidal:
		a1 = 0.1 # First Amplitude
		phi = -4.0 # sin(2*pi*t/P) + phi <-This is phi. Offset; horizontal shift
		orbital_period = 6.4714 # According to Pavlenko et al (1996)
		a2 = 0.2 # Second Amplitude
		diff = 16.6 # Shift due to average magnitude

		p = [a1, phi, orbital_period, diff, a2] #This is for my attempt

		##p = [a1, phi, orbital_period, diff] #for the equation in line 33

		pGuess = np.copy(p)

	# make bounds for the chunks
	tBounds = makeBounds(jd)
	chunks = classifyChunks(jd, tBounds)

	tLow, yLow = assignLowerEnvelope(jd, mag, chunks, pctile, \
						 useMags=useMags, \
						 clipOutliers=clipOutliers)

	if iCheck > -1:
		bCheck = chunks == iCheck
		if plotHistogram:	
			showHist(mag[bCheck])


 	if plotEllipsoidal:
	 	p0 = np.array([a1, phi, orbital_period, a2, diff])
	 	p1, success = \
		    optimize.leastsq(errFunc, pGuess[:], args=(jd, mag), \
					     maxfev=int(1e6), ftol=1e-10)
	 	pLow, successLow = \
		    optimize.leastsq(errFunc, pGuess[:], args=(tLow, yLow), \
					     maxfev=int(1e6), ftol=1e-10)

	 	# by this point we have the ellipsoidal modulation fit to the
	 	# dataset
	 	ySub = mag - twoSine(pLow, jd)

		# write out the ellipsoidal modulation model to disk for
		# uniform characterization by v404ellipsoidal.py
		tModel = Table()
		tFine = np.linspace(np.min(jd), np.max(jd)+8., 400)
		pFine, _ = phaseFromJD(tFine)
		yFine = twoSine(pLow, tFine)
		ll = np.argsort(pFine)

		tModel['x'] = pFine[ll]
		tModel['y'] = yFine[ll]
	##filEll = 'mdm_fig2_2017.txt'  # use the fig2_yyyy.txt convention
	##if os.access(filEll, os.R_OK):
		##os.remove(filEll)
	##tModel.write(filEll, format='ascii')

 	# by THIS point we have our ellipsoidal-subtracted
 	# lightcurve. Now let's try binning this...
 	if showBinned:
	 	tBin, fBin, uBin, nBin = \
		    BinData(jd, mag, dy, tStart=tStart, tEnd=tEnd, \
				    BinTime=binTime, plotDBG=True)

	 	print "binned arrays shape, nMin, nMax:", \
		    np.shape(tBin), np.min(nBin), np.max(nBin)

 	# if you want to subtract the ellipsoidal modulation from the
 	# binned data, you might do:
 	if plotSubtractedData:
	 	fBinSub = fBin - twoSine(pLow, tBin)

	 	# or just the raw data...
	 	fSub = mag - twoSine(pLow, jd)

	 	# let's plot this...


 		fig11 = plt.figure(11)
 		fig11.clf()
 		ax11 = fig11.add_subplot(111)
 		ax11.scatter(jd, fSub)

		# 2018-04-07 WIC - write the binned subtracted lightcurve to
		# file to plot with other methods.
		tExp = Table()
		tExp['tBin'] = tBin
		tExp['fBin'] = fBin
		tExp['uBin'] = uBin
		tExp['nBin'] = nBin
		tExp['fBinSub'] = fBinSub
		tExp['phsBin'], _  = phaseFromJD(tBin)

		# export this to disk
		filExp='v404_binSub.fits'
		if os.access(filExp, os.R_OK):
			os.remove(filExp)
		tExp.write(filExp)

 	orbital_period = 6.4714
 	period = np.linspace(0.0005, orbital_period, 1000)
 	omega = 2 * np.pi / period
 	if plotBinnedLS:
 		binnedLS = lomb_scargle(tBin, fBinSub, uBin, omega, generalized=False)
 	
 	if bootstrapping == True:
 		Dbin = lomb_scargle_bootstrap(tBin, fBinSub, uBin, omega, generalized=False, N_bootstraps=1000, random_state=0)
 		sig1b, sig5b = np.percentile(Dbin, [99,95])

 	if plotBinnedData:
 		fig4 = plt.figure(4)
 		fig4.clf()
 		ax4 = fig4.add_subplot(111)
 		ax4.plot(tBin, fBinSub, 'ko')
 		plt.errorbar(tBin, fBinSub, yerr=uBin, fmt='o', ms=4, ecolor='0.5', alpha=0.5)

 	
 	if plotBinnedLS:
 		fig5 = plt.figure(5)
 		fig5.clf()
		#ax5 = fig5.add_subplot(111, xscale='log', yscale='log')
		ax5 = fig5.add_subplot(111)
		ax5.plot(period, binnedLS, 'ko', ls='-', ms=4)
		if bootstrapping:
	   		ax5.plot([period[0], period[-1]], [sig1b, sig1b], ':', c='red') #sig1 means percentile is '99'
	   		ax5.plot([period[0], period[-1]], [sig5b, sig5b], ':', c='green') #sig5 means percentile is '95'

 	##period = orbital_period
 	period = np.linspace(0.0005, orbital_period, 1000)

 	# let's try logspace. Here we specify the limits in terms of a power of ten.
 	if plotLS:
	 	period = np.logspace(-3., 0.5, 1000)
	 	#dy = tbl['fluxErr']
	 	omega = 2 * np.pi / period
		# 	flicker = lomb_scargle(period, ySub, dy, omega, generalized=False)
	 	flicker = lomb_scargle(jd, ySub, dy, omega, generalized=False)
	 	num2 = len(flicker)
	 	if bootstrapping == True:
	 		D = lomb_scargle_bootstrap(jd, ySub, dy, omega, generalized=False, N_bootstraps=1000, random_state=0)
	 		sig1, sig5 = np.percentile(D, [99,95])

		 		#Attempting to plot the periodogram


		fig3 = plt.figure(3)
		fig3.clf()
		#ax3 = fig3.add_subplot(111, xscale='log', yscale='log')
		ax3 = fig3.add_subplot(111)

		ax3.plot(period, flicker, 'ko', ls='-', ms=4)
		if bootstrapping == True:
			ax3.plot([period[0], period[-1]], [sig1, sig1], ':', c='red') #sig1 means percentile is '99'
			ax3.plot([period[0], period[-1]], [sig5, sig5], ':', c='green') #sig5 means percentile is '95'

 	if plotNoiseData:
	 	magsNoise = np.random.normal(size=np.size(jd))*dy
	 	if plotNoiseLS:
		 	flickNoise = lomb_scargle(jd, magsNoise, dy, omega, generalized=False)

		    #Plot the LS
		
			fig6 = plt.figure(6)
			fig6.clf()
			#ax5 = fig5.add_subplot(111, xscale='log', yscale='log')
			ax6 = fig6.add_subplot(111)
		 	# Get Significance via bootstrap
	 		if bootstrapping == True:
	 			DNoise = lomb_scargle_bootstrap(jd, magsNoise, dy, omega, generalized=False, N_bootstraps=1000, random_state=0)
				sig1N, sig5N = np.percentile(DNoise, [99,95])


			ax6.plot(period, flickNoise, 'ko', ls='-', ms=4)
			if bootstrapping == True:
				ax6.plot([period[0], period[-1]], [sig1N, sig1N], ':', c='red') #sig1 means percentile is '99'
				ax6.plot([period[0], period[-1]], [sig5N, sig5N], ':', c='green') #sig5 means percentile is '95'

		# Plot the random noise

		fig7 = plt.figure(7)
		fig7.clf()
		ax7 = fig7.add_subplot(111)

		ax7.scatter(jd, magsNoise)


		#Let's try binning the random noise
		if plotBinnedNoiseData:
			tBinNoise, fBinNoise, uBinNoise, nBinNoise = BinData(jd, magsNoise, dy, tStart=tStart, tEnd=tEnd, BinTime=binTime, plotDBG=True)
		 	print "binned arrays shape, nMinNoise, nMaxNoise:", np.shape(tBinNoise), np.min(nBinNoise), np.max(nBinNoise)

		 	#Attempting to plot the binned noise

	 		fig8 = plt.figure(8)
	 		fig8.clf()
	 		ax8 = fig8.add_subplot(111)
	 		ax8.plot(tBinNoise, fBinNoise, 'ko')
	 		plt.errorbar(tBin, fBinNoise, yerr=uBinNoise, fmt='o', ms=4, ecolor='0.5', alpha=0.5)

	 	#Lomb-Scargle for binned noise

	 	if plotBinnedNoiseLS:
		 	orbital_period = 6.4714
		 	period = np.linspace(0.0005, orbital_period, 1000)
		 	omega = 2 * np.pi / period
		 	binnedLSNoise = lomb_scargle(tBinNoise, fBinNoise, uBinNoise, omega, generalized=False)
		 	if bootstrapping == True:
		 		Dnb = lomb_scargle_bootstrap(tBinNoise, fBinNoise, uBinNoise, omega, generalized=False, N_bootstraps=1000, random_state=0)
		 		sig1nb, sig5nb = np.percentile(Dnb, [99,95])


		 	#Plot this

	 		fig9 = plt.figure(9)
	 		fig9.clf()
	 		ax9 = fig9.add_subplot(111)

	 		ax9.plot(period, binnedLSNoise, 'ko', ls='-', ms=4)
	 		if bootstrapping == True:
	 			ax9.plot([period[0], period[-1]], [sig1nb, sig1nb], ':', c='red') #sig1 means percentile is '99'
				ax9.plot([period[0], period[-1]], [sig5nb, sig5nb], ':', c='green') #sig5 means percentile is '95'


    #blah = optimize.leastsq(errFunc, p0[:], args=(Tx, tX), full_output=True)
    #print blah

 	##print p1, success
 	##print pLow, successLow



	if plotEllipsoidal:
		yPred = twoSine(p1, jd)
		yPredLow = twoSine(pLow, jd)

		#Generating random points throughout the best fit line
		#num_points = 1401
		#randfunc = 0.0393104 * np.sin(2.0*np.pi*jd/6.47155519  + -1.68599974)+16.5961013 + 0.15487321 * np.sin(4.0*np.pi*jd/6.47155519 + -1.68599974)*np.random.rand(num_points)

		xGrid = np.linspace(57985., 57991., 1000)
		#randfunc = twoSine(pLow,xGrid)

		# Plotting the v404 data and the best fit line

		##print tBounds
		#	print np.shape(jd)
		#print np.shape(chu#nks)
		#return

		# 2018-04-07 WIC - modified Austin's plotting stanza to show
		# the phase curve.

		# create a few phase bins
		phaseLow, _  = phaseFromJD(tLow)
		phaseGrid, _  = phaseFromJD(xGrid)

		# Because the phasing is not the same order as the jd, we need
		# to sort the line-drawn objects so that we don't get
		# crossover lines
		lG = np.argsort(phaseGrid)

		# plot by jd or by phase?
		tSho = np.copy(jd)
		tLo = np.copy(tLow)
		tGrid = np.copy(xGrid)
		cScatt = flagColor
		if showPhase:
			tSho = np.copy(phase)
			tLo = np.copy(phaseLow)
			tGrid = np.copy(phaseGrid)
			cScatt = jd - np.min(jd)

		plt.figure(1)
		plt.clf()
		#plt.scatter(jd, mag, alpha=0.5, color='darkmagenta')
		dum = plt.scatter(tSho, mag, \
					  alpha=1., c=cScatt, s=16, \
					  cmap='inferno', zorder=25, \
					  edgecolor='0.4')
		plt.plot(tLo, yLow, 'ko', ms=7, zorder=25)
		plt.plot(tGrid[lG], twoSine(pLow,xGrid[lG]), c='k')
		plt.plot(tGrid[lG], twoSine(p1,xGrid[lG]), c='g', ls='--')
		# plt.scatter(jd, ySub, c='violet', s=16)
		if errorbars:
			plt.errorbar(tSho, mag, yerr=dy, fmt='o', \
					     ms=1, ecolor='0.3', alpha=0.5, \
					     zorder=10)
		
		#plt.plot(phase, yPred, 'gx', lw=2, alpha=0.5, ls='-')
		#plt.plot(phase, yPredLow, 'k+', lw=2, alpha=0.5, ls='--')
		if useFlag > 4 or usePhase:
			plt.colorbar(dum)

		##plt.plot(jd, 0.0393104 * np.sin(2.0*np.pi*jd/6.47155519  + -1.68599974)+16.5961013 + 0.15487321 * np.sin(4.0*np.pi*jd/6.47155519 + -1.68599974), '-b') #attempt to plot the function itself

		plt.show(block=False) # block=False to ensure the terminal doesn't hang
		#return

		yLims = np.copy(plt.gca().get_ylim())

		if useMags:
			plt.ylim([yLims[1], yLims[0]])

		# let's show the bounds
		if not showPhase:
			for iBound in range(np.size(tBounds)):
				tThis = tBounds[iBound]
				plt.plot([tThis, tThis], yLims, 'k--')

		# And actually label the axes
		if showPhase:
			plt.xlabel('Phase (6.4714 d period)')
		else:
			plt.xlabel('JD - 2 400 000.0 d')
		
		if useMags:
			plt.ylabel('R (mag)')
		else:
			plt.ylabel('Flux / F(ref)')
		

		# now set the plot limits
		if useMags:
			plt.ylim(17.0, 16.5)
		if showPhase:
			plt.xlim(0., 1.)
		
		# 2018-04-07 - save the plot as png, naming convention after
		# the figure numbers in the routine.
		#plt.savefig('v404_fig1.png')

		##ax3.show(block=False)

	plt.figure(0)
	plt.clf()
	dum0 = plt.scatter(jd, mag,\
					  alpha=1., c=flagColor, s=16, \
					  cmap='inferno', zorder=25, \
					  edgecolor='0.4')
	if errorbars:
			plt.errorbar(jd, mag, yerr=dy, fmt='o', \
					     ms=1, ecolor='0.3', alpha=0.5, \
					     zorder=10)
	plt.show(block=False)
	
	yLims = np.copy(plt.gca().get_ylim())

	if useMags:
		plt.ylim([yLims[1], yLims[0]])	


def oneSine(p,x):

	"""Single sine wave"""

	return p[0] * np.sin(2.0*np.pi*x/p[2] + p[1]) + p[3]

def twoSine(p,x):

	"""'Hardcoded' double sine"""

	sineOne = p[0] * np.sin(2.0*np.pi*x/p[2] + p[1]) + p[3]
	sineTwo = p[4] * np.sin(4.0*np.pi*x/p[2] + p[1])

	return sineOne + sineTwo

def twiceSine(p,x, debug=False):

	"""Constructs a sum of two sines from a single sine"""

	pOne = p[0:4]
	pTwo = np.copy(pOne)
	pTwo[0] = p[-1]
	pTwo[2] = pOne[2] * 0.5  # half the period
	pTwo[3] = 0. # don't want to add the diff twice!

	if debug:
		print "twoSine DEBUG:"
		print pOne
		print pTwo

	sineOne = oneSine(pOne, x)
	sineTwo = oneSine(pTwo, x)

	return sineOne + sineTwo

	##return p[0] * np.sin(2.0*np.pi*(x-p[1])/p[2]) + p[3] # Taken from the termninal archive

	##return ((p[0]*np.sin(((2*np.pi*x)+p[1])/p[2]))+(p[3]*np.sin(((2*np.pi*x)+p[1]/0.5*p[2]))))+p[4] # My attempt at the function on the board...

def twoSine(p,x):

	"""'Hardcoded' double sine"""

	sineOne = p[0] * np.sin(2.0*np.pi*x/p[2] + p[1]) + p[3]
	sineTwo = p[4] * np.sin(4.0*np.pi*x/p[2] + p[1])

	return sineOne + sineTwo

def errFunc(p, x, y):

    """The error function"""

    return twoSine(p, x) - y

def makeBounds(times=np.array([]), buf=0.3, interval=1.):

	"""Makes boundaries for a times array"""

	tMin = np.min(times) - np.abs(buf)
	tMax = np.max(times) + np.abs(buf) + interval
	borders = np.arange(tMin, tMax, interval)

	return borders

def classifyChunks(times=np.array([]), bounds=np.array([])):

	"""Classifies time points by which chunk they appear in"""

	whichChunk = np.zeros(np.size(times), 'int')

	for iChunk in range(np.size(bounds)-1):
		tMin = bounds[iChunk]
		tMax = bounds[iChunk+1]

		bThis = (times >= tMin) & (times < tMax)

		# now we assign the chunk ID to all the objects in this chunk
		whichChunk[bThis] = int(iChunk)

	return whichChunk

def assignLowerEnvelope(t=np.array([]), y=np.array([]), chunks=np.array([]), pct=10., useMags=True, \
		clipOutliers=True):

	"""Assigns the lower envelopes to the t,y array , partitioned by chunks, using the pct'th percentile """

	# we will need average time, mag arrays at least
	nChunks = np.max(chunks) - np.min(chunks)+1
	tEnv = np.zeros(nChunks)
	yEnv = np.copy(tEnv)

	# now we populate:
	for iChunk in range(nChunks):

		# which points are in this chunk?
		bChunk = chunks == iChunk

		if np.sum(bChunk) < 1:
			continue

		# find the envelope for THIS chunk...
		tMed, yMed = findLowerValue(t[bChunk], y[bChunk], pct, useMags=useMags, clipOutliers=clipOutliers)
		tEnv[iChunk] = tMed
		yEnv[iChunk] = yMed

	return tEnv, yEnv



def findLowerValue(t=np.array([]), y=np.array([]), pctile=10., useMags=True, restrictMedian=False, clipOutliers=True, clipSigma=3., clipIters=3):

	"""Return the pctile'th value of y"""

	pctUse = np.copy(pctile)
	if useMags:
		pctUse = 100.-pctUse

	yUse = np.copy(y)
	if clipOutliers:
		yClip = sigma_clip(y, sigma=clipSigma, iters=clipIters)
		bKeep = ~yClip.mask
		yUse = yUse[bKeep]

	# find the percentile value
	pctVal = np.percentile(yUse, pctUse)

	# find all points that pass this criterion
	if useMags:
		bUse = y > pctVal
	else:
		bUse = y < pctVal

	# now compute the median time
	if restrictMedian:
		tMed = np.median(t[bUse])
	else:
		tMed = np.median(t)

	return tMed, pctVal	

def BinData(vTime=np.array([]), vRate=np.array([]), vError=np.array([]), nMin=2, tStart=-1e9, tEnd=57992., BinTime=0.0034722, \
	Verbose=True, plotDBG=False):

    """Bin data given a time-series.

    RETURNS binned time, rate, uncertainty, and the number per bin"""

	# Initialize output arrays:
    vOutTime = np.array([])
    vOutRate = np.array([])
    vOutError = np.array([])
    vNPerBin = np.array([])

    if np.size(vTime) < 2:
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
        print "BinTime = %.2f" % (BinTime)
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

    print np.shape(vBinIDs)

    for iBin in range(0, np.size(vBinIDs)):
        
        # When does this bin start and end?
        ThisTStart = vBinTimes[iBin]
        ThisTEnd   = ThisTStart + BinTime

        # ignore data > tend
        if ThisTEnd >= tEnd:
            continue

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

    # At thsi point we've been through the loops. Now we send the
    # three vectors we've produced to the routine that called it. So:
    ##np.shape(vOutTime), np.shape(vOutRate), np.shape(vOutError)
    ##np.savetxt('BLAH.lis', np.transpose(np.vstack((vOutTime,vOutRate,vOutError))))
    ##Time, Rate, Error = np.loadtxt('BLAH.lis', unpack=True)
    ##return vOutTime, vOutRate, vOutError, vNPerBin

    # QQQ probably not a good idea to do the subtraction in the binning method... 

  #   tBounds = makeBounds(vOutTime)
  #   chunksBin = classifyChunks(vOutTime, tBounds)
  #   chunks = classifyChunks(vOutTime, tBounds)
  #   pctile = 10.
  #   useMags = True
  #   clipOutliers = True

  #   a1 = 0.1 # First Amplitude
  #   phi = -4.0 # sin(2*pi*t/P) + phi <-This is phi. Offset; horizontal shift
  #   orbital_period = 6.4714 # According to Pavlenko et al (1996)
  #   a2 = 0.2 # Second Amplitude
  #   diff = 16.6 # Shift due to average magnitude

  #   tLowb, yLowb = assignLowerEnvelope(vOutTime, vOutRate, chunks, pctile, useMags=useMags, clipOutliers=clipOutliers)

  #   p0b = np.array([a1, phi, orbital_period, a2, diff])
  #   p1b, success2 = optimize.leastsq(errFunc, pGuess[:], args=(vOutTime, vOutRate), maxfev=int(1e6), ftol=1e-10)
  #   pLow2, successLow2 = optimize.leastsq(errFunc, pGuess[:], args=(tLow, yLow), maxfev=int(1e6), ftol=1e-10)

 	# # by this point we have the ellipsoidal modulation fit to the dataset
  #   ySub2 = vOutRate - twoSine(pLow2, vOutTime)



    # return vOutTime, vOutRate, vOutError, vNPerBin



def showHist(x=np.array([]), bins=20, checkClipping=True):

	"""Utility - show a histogram"""

	if np.size(x) < 2:
		return

	# enforce the same range
	xRange = [np.min(x), np.max(x)]


	fig2 = plt.figure(2)
	fig2.clf()
	ax = fig2.add_subplot(111)
	dum = ax.hist(x, bins, alpha=0.5, range=xRange)

	# now we sigma-clip
	xClip = sigma_clip(x, sigma=3, iters=2)
	dum2 = ax.hist(xClip[~xClip.mask], bins, alpha=0.5, range=xRange)

	return

def genData(nPoints=1000, nNights=6, Bootstrapping=False):

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

    #Lomb-Scargle for this data



    # to generate the random noise magnitudes:
    # magsNoise = np.random.normal(size=np.size(time))*unctys

    # orbital_period = 6.4714
    # period = np.linspace(0.0005, orbital_period, 1000)
    # omega = 2 * np.pi / period
    # randLS = lomb_scargle(time, magsNoise, unctys, omega, generalized=False)

    return times[lSor], unctys[lSor], mags[lSor]

def correctMagForContaminant(magObs=20., magContam=17.2):

	"""Utility for correcting v404 cyg's apparent magnitude for
	contamination."""
	
	# Moved out into a separate method so that we can use it
	# elsewhere (e.g. when plotting ellipsoidal modulation),
	# possibly from other modules or the IPython command line..

	# Notice that the flux of whatever reference object we used
	# divides from both sides of the expression (flux obs = flux
	# obj + flux contam), so we can ignore it. I break the
	# expression into steps here to aid debugging if needed.

	fObs = 10.0**(-0.4*magObs)
	fCon = 10.0**(-0.4*magContam)

	# so, what's the relative flux of the object itself?
	return -2.5*np.log10(fObs - fCon)

def phaseFromJD(jdShort=np.array([]), per=6.4714, tZer=48813.873, \
			u_per = 0.0001, u_tZer=0.004):

	"""Compute the phase from given times. Note that the times are
	provided as JD - 2 400 000.0"""

	# We'll do this in pieces so that debugging is easier:
	dt = jdShort - tZer

	# compute the phase...
	nOrbs = dt/per
	phs = nOrbs - np.floor(nOrbs)

	# ... the fractional uncertainty in nOrbs. If it's greater
	# than one orbit then we really don't know the phase very well
	# at all! 
	frac_var_nOrbs = (u_tZer/dt)**2 + (u_per/per)**2
	u_phs = nOrbs * np.sqrt(frac_var_nOrbs)

	return phs, u_phs

# def showRawCounts(tPhot=Table(), figNam='test_rawCounts.png'):

# 	"""Show the raw counts"""

# 	if len(tPhot) < 1:
# 		return

# 	jd = tPhot['J.D.-2400000']
# 	countsObj = tPhot['Source-Sky_T1']
# 	countsCom = tPhot['Source-Sky_C2']
# 	flag = tPhot['Flag']


# 	fig1 = plt.figure(10)
# 	fig1.clf()
# 	ax1=fig1.add_subplot(211)
# 	ax2=fig1.add_subplot(212, sharex=ax1, sharey=ax1)

# 	fig1.subplots_adjust(hspace=0.05)

# 	pSyms = ['x', '+', 'x', 'v', 's']
# 	pCols = ['b', 'g', 'b', 'g', '0.4']

# 	# loop through the flag quality:
# 	for fl in range(np.max(flag)+1):
# 		bThis = np.abs(fl-flag) < 1.0e-2
# 		if np.sum(bThis) < 1:
# 			continue
	
# 		# which plot symbol?
# 		pSym = pSyms[fl % len(pSyms)]
# 		colo = pCols[fl % len(pCols)]

# 		dum1 = ax1.scatter(jd[bThis], countsObj[bThis], \
# 				   marker=pSym, \
# 				   color=colo, \
# 				   alpha=0.5, \
# 				   s=25)
# 		dum2 = ax2.scatter(jd, countsCom, \
# 				   marker=pSym, \
# 				   color=colo, \
# 				   alpha=0.5, \
# 				   s=25)

# 	ax2.set_xlabel('JD - 2 400 000.0 d')
# 	ax1.set_ylabel('Source minus sky, counts')
# 	ax2.set_ylabel('Source minus sky, counts')

# 	# standardize the y range
# 	ax1.set_ylim(0., 15000.)
# 	#ax2.set_ylim(0., 15000.)

# 	ax1.annotate('v404 Cyg + contaminant', (0.95, 0.95), \
# 		     ha='right', va='top', \
# 		     xycoords='axes fraction', \
# 		     color='k', fontsize=14)

# 	ax2.annotate('Comparison star "C4" (R = 16.07)', (0.95, 0.95), \
# 		     ha='right', va='top', \
# 		     xycoords='axes fraction', \
# 			     color='k', fontsize=14)


# 	#save to disk
# 	fig1.savefig(figNam)

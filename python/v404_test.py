# All imports from terminal first

import matplotlib.pylab as plt
plt.ion()
plt.style.use('ggplot')
from astropy.table import Table
import numpy as np
from scipy import optimize
import sys
from astropy.stats import sigma_clip
import loadOld
from astroML.time_series import lomb_scargle, lomb_scargle_bootstrap

# for testing existence of files
import os

# Reading the tables

try:
	t92 = loadOld.loadPhot(fileIn="/Users/Shared/Data/v404Cyg_zurita04/199208elip_R.dat")
	t98 = loadOld.loadPhot(fileIn="/Users/Shared/Data/v404Cyg_zurita04/elipR1998.dat")
	t99A = loadOld.loadPhot(fileIn="/Users/Shared/Data/v404Cyg_zurita04/mv4_990706_apcorr.dat")
	t99B = loadOld.loadPhot(fileIn="/Users/Shared/Data/v404Cyg_zurita04/mv4_990707_apcorr.dat")
	t03 = loadOld.loadPhot(fileIn="/Users/Shared/Data/v404Cyg_zurita04/28_July_2003_iac80.qdp")
	t17 = loadOld.loadPhot(fileIn="/Users/amblevin/Desktop/G4_2017.txt")
	t18A = loadOld.loadPhot(fileIn="/Users/amblevin/Desktop/G4_2018A.txt")
	t18B = loadOld.loadPhot(fileIn="/Users/amblevin/Desktop/2018B_Sorted.txt")


	#t17 = loadOld.loadPhot(fileIn="../../../MDM/austinsFiles/v404_2017.txt")
	#t18A = loadOld.loadPhot(fileIn="../../../MDM/austinsFiles/v404_2018A.txt")

except UnboundLocalError:
	print "One or more files not found. Check the directories indicated in the source code."

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

# 2019-02-11 HARDCODE
#tZer = 57964.4326

def go(pctile=10., iCheck=1, useMags=True, \
	       clipOutliers=False, \
	       #oldAperture=False, plotRawCounts=False, \
	       useFlag=12, \
	       tStart=-1e9, \
	       tEnd = 57993., \
	       #tEnd=57992.0, \
	       binTime=0.0034722, bootstrapping=False, \
	       errorbars=True, \
	       magContam=17.52, magCompar=16.07, \
	       showPhase=True, \
	       writeOnly=False, \
	       data='2017', \
	       noCorrection=False, showBinned=False, \
	       plotBinnedData=False, \
	       plotBinnedLS=False, plotLS=False, \
	       plotNoiseData=False, \
	       plotNoiseLS=False, plotHistogram=False, \
	       plotBinnedNoiseData=False, plotBinnedNoiseLS=False, \
	       plotSubtractedData=False, plotEllipsoidal=False, limit=0.1, \
	       plotBinnedOnSubtracted=False, writeEllipsoidal=False, \
	       overlayEllipsoidal=False, compareEllipsoidals=False, \
	       genOS=False, moreTimes=False, genTS=False, amp1=-0.04, amp2=0.108, \
	       plot2Phases=False, \
	       filBinned='v404_binSub.fits', \
	       filExp='v404_unbinned_sub.fits'):
	# WIC - put the table reading back into go, to avoid scope
	# confusion

	warn = "WARNING: Ellipsoidal modulation will not work for 1999 and 2003!"
	warn2 = "WARNING: Binned and subtracted data only plots if plotEllipsoidal is set to True."
	#warn3 = "2017 and 2018A data will not work for this version of v404_test on Christian's account."

	if not plotEllipsoidal:
		if showBinned or plotSubtractedData:
			print warn2
		if compareEllipsoidals:
				print "plotEllipsoidal must be set to True in order to compare ellipsoidals."


	if plotBinnedLS or plotBinnedNoiseLS or plotLS:
		print "WARNING: LS plots only work if the data the LS is being performed on is also plotted."

	if data == '1992':
		tbl = t92
		#print tbl
		t0 = projectEphemeris()
	elif data == '1998':
		tbl = t98
		t0 = projectEphemeris()
		print "WARNING: Due to 1998 data already being in phase, binning currently does not work."
	elif data == '1999A':
		tbl = t99A
		t0 = projectEphemeris()
		print warn
	elif data == '1999B':
		tbl = t99B
		t0 = projectEphemeris()
		print warn
	elif data == '2003':
		tbl = t03
		t0 = projectEphemeris()
		print warn
	elif data == '2017':
		tbl = t17
		oldT0 = np.min(tbl['time'])
		t0 = projectEphemeris(minJD=oldT0)
		#print warn3
	elif data == '2018A':
		tbl = t18A
		oldT0 = np.min(tbl['time'])
		t0 = projectEphemeris(minJD=oldT0)
		#print warn3
	elif data == '2018B':
		tbl = t18B
		oldT0 = np.min(tbl['time'])
		t0 = projectEphemeris(minJD=oldT0)
	else:
		print "The value for 'data' must be either 1992, 1998, 1999A, 1999B, 2003, 2017, 2018A or 2018B."
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

	if data == '2018B':
		ucty = tbl['fluxErr']
		goodErr = ucty < limit



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
		if data == '1992' or data == '1998' or data == '1999A' or data == '1999B' or data == '2003':
			print "WARNING: if using data from Z04, it will be INCORRECT unless noCorrection=True!"
	else:
		mag = np.copy(magBoth)
	dy = np.copy(errMag)

	if 'time' in tbl.colnames:
		jd = tbl['time']
		# compute the phase on the orbital ephemeris.
		phase, u_phase = phaseFromJD(jd, tZer=t0)
	else:
		phase = tbl['phase']
		jd = tbl['phase'] #so the rest of the script doesn't return an error.........

	# 2018-04-07 useful to plot the RAW data. Rather than
	# interrupt the flow, we port this off to another method.
	
	#if plotRawCounts:
		#showRawCounts(tbl)

	#mag = tbl['Average Mag(V404)']  # superseded by the material above.
	#dy = 1.086 * tbl['rel_flux_err_T1']
	num = len(jd)


	if not useMags:
		mag = np.copy(relFlux)
		dy = np.copy(errFlux)
		#mag = tbl['rel_flux_T1'] 
		#dy = tbl['rel_flux_err_T1']

	### 2018-06-01 WIC - throw in an argument to just plot for one
	### night only
	# if writeOnly:
	# 	plt.style.use('ggplot')
	# 	fig2 = plt.figure(2, figsize=(10,4))
	# 	fig2.clf()

	# 	#print jd

	# 	#print np.shape(jd)
	# 	#print np.shape(mag)

	# 	tBin, fBin, uBin, nBin = \
	# 	    BinData(jd, mag, dy, tStart=tStart, tEnd=tEnd, \
	# 			    BinTime=binTime, plotDBG=True)


	# 	xSho = np.copy(jd)
	# 	ySho = np.copy(mag)
	# 	eSho = np.copy(dy)

	# 	if showBinned:
	# 		xSho = np.copy(tBin)
	# 		ySho = np.copy(fBin)
	# 		eSho = np.copy(uBin)

	# 	ax1 = fig2.add_subplot(111)
	# 	#dum = ax1.plot(jd, mag, 'bo')
	# 	dum2 = ax1.errorbar(xSho, ySho, eSho, fmt='bo', ls='None', \
	# 				    ms=1, alpha=0.5)

	# #	dum2 = ax1.scatter(tBin, fBin, color='b', alpha=0.5, s=4)


	# 	#dum2 = ax1.errorbar(tBin, fBin, uBin, fmt='b.', ls='None')

	# 	yMin = np.percentile(fBin, 1)

	# 	yRang = np.copy(ax1.get_ylim())
	# 	ax1.set_ylim(yMin+0.4, yMin)
		
	# 	### 2018-06-02 WIC - thrown-together estimate for the
	# 	### chisq
	# 	bLo = xSho < 58270.953695
	# 	fig3 = plt.figure(3, figsize=(5,5))
	# 	fig3.clf()
	# 	ax3 = fig3.add_subplot(111)

	# 	chi = (ySho - np.median(ySho))/eSho

	# 	dum = ax3.hist(chi, 100)
	# 	print "INFO - std(chi) = %.3f" % (np.std(chi))

	# 	# also print out the uncertainties
	# 	print "INFO: dy median, std: %.3e, %.3e" \
	# 	    % (np.median(eSho), np.std(eSho))

	# 	ax1.plot(xSho[bLo], ySho[bLo], 'r.', zorder=5)

	# 	return


	# This was my (failed) attempt to use the lambda function for
	# y like in the example


	# y = lambda p, x: ((p[0]*np.sin(((2*np.pi*x)+p[1])/p[2]))+(p[3]*np.sin(((2*np.pi*x)+p[1]/0.5*p[2]))))+p[4]

	# Initial Guess for Parameters

	if plotEllipsoidal:
		if not useMags:
			print "WARNING: Lower-envelope ellipsoidals will not plot while useMags=False."
		a1 = 0.1# 0.02 # First Amplitude -- 'DEFAULT' VALUE IS 0.1
		phi = -4.0 # sin(2*pi*t/P) + phi <-This is phi. Offset; horizontal shift-- Default is -4.0
		orbital_period = 6.4714 # According to Pavlenko et al (1996)
		diff = 16.6 # Shift due to average magnitude- Default is 16.6
		a2 = 0.2 # Second Amplitude -- 'DEFAULT' VALUE IS 0.2 
		

		p = [a1, phi, orbital_period, diff, a2] #This is for my attempt
		#p0, p1, p2, p3, p4 = a1, phi, orbital_period, diff, a2

		##p = [a1, phi, orbital_period, diff] #for the equation in line 33

		pGuess = np.copy(p)

		# make bounds for the chunks
		tBounds = makeBounds(jd)
		chunks = classifyChunks(jd, tBounds)

		#tLow, yLow = assignLowerEnvelope(jd, mag, chunks, pctile, \
		#					 useMags=useMags, \
		#					 clipOutliers=clipOutliers)

		tLow, yLow = assignLowerEnvelope(jd, mag, chunks, data=data)

		if iCheck > -1:
			bCheck = chunks == iCheck
			if plotHistogram:	
				showHist(mag[bCheck])


	 	p0 = np.array([a1, phi, orbital_period, diff, a2])
	 	guess = np.array([a1, phi, diff, a2])
	 	bounds = ([-np.inf, -np.inf, -np.inf, -np.inf], [np.inf, np.inf, np.inf, np.inf])
	 	#print "Guess[a1, phi, o_p, diff, a2]: ", p0
	 	if overlayEllipsoidal:
	 		p1 = np.loadtxt("/Users/amblevin/Desktop/p12018A.txt", unpack=True)
	 		pLow = np.loadtxt("/Users/amblevin/Desktop/pLow2018A.txt", unpack=True)
	 	else:
	 		ph, _ = phaseFromJD(jd, tZer=t0)
	 		phLow, _ = phaseFromJD(tLow, tZer=t0)
	 		print "phLow: ", phLow
	 		print "tLow: ", tLow
	 		print "yLow: ", yLow
			p1, pcov_1 = optimize.curve_fit(twoSinePhase, ph, mag, p0=guess, method='trf', \
		    	bounds=bounds, sigma=dy, absolute_sigma=True) 
			
			#pLow, successLow = optimize.leastsq(errFunc, pGuess[:], args=(tLow, yLow), maxfev=int(1e6), ftol=1e-10)

			pLow, pcov_low = optimize.curve_fit(twoSinePhase, phLow, yLow, p0=guess, method='trf', \
		    	bounds=bounds) # Sigma CANNOT be used for the lower envelope!

		    
		if writeEllipsoidal:
			np.savetxt('p1'+str(data)+'.txt', p1)
			np.savetxt('pLow'+str(data)+'.txt', pLow)
			print "Saved p1: ", p1
			print "Saved pLow: ", pLow

		
		print "Fit[a1, phi, o_p, diff, a2]: ", p1
		#print "Lower Envelope: ", pLow

	 	# by this point we have the ellipsoidal modulation fit to the dataset
	 	# 2019-03-19 altered twoSineTime --> twoSine2019
	 	ySub = mag - twoSine2019(jd, *pLow)

		# write out the ellipsoidal modulation model to disk for
		# uniform characterization by v404ellipsoidal.py
		tModel = Table()
		tFine = np.linspace(np.min(jd), np.max(jd)+8., 400)
		pFine, _ = phaseFromJD(tFine, tZer=t0)
		yFine = twoSine2019(tFine, *pLow, t0=t0)
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

 	# if you want to subtract the ellipsoidal modulation from the data, you might do:
 	if plotSubtractedData:
	 	if showBinned:
	 		fBinSub = fBin - twoSine2019(tBin, *pLow)
	 	else:
	 		fSub = mag - twoSine2019(jd, *pLow)

	 	# let's plot this...

 		fig11 = plt.figure(11)
 		fig11.clf()
 		ax11 = fig11.add_subplot(111)
 		if showBinned:
 			ax11.scatter(tBin, fBinSub)
 			if errorbars:
 				ax11.errorbar(tBin, fBinSub, yerr=uBin, fmt='o', ms=4, ecolor='0.5', alpha=0.5)
 		else:
 			ax11.scatter(jd, fSub)
 			if errorbars:
 				ax11.errorbar(jd, fSub, yerr=dy, fmt='o', ms=4, ecolor='0.5', alpha=0.5)
 		if plotBinnedOnSubtracted:
 			ax11.plot(tBin, fBinSub, 'bo')
 			ax11.errorbar(tBin, fBinSub, yerr=uBin, fmt='o', ms=4, ecolor='0.5', alpha=0.5)

 		yLims = np.copy(plt.gca().get_ylim())
 		plt.ylim([yLims[1], yLims[0]])



		# 2018-04-07 WIC - write the binned subtracted lightcurve to
		# file to plot with other methods.

		# 2019-03-19 - recommend including extra syntax like this to also output the un-binned data table
		# (lcPlot.py should work just fine as long as the column header names are the same).
		if writeOnly:
			if showBinned:
				tExp = Table()
				tExp['tBin'] = tBin
				tExp['fBin'] = fBin
				tExp['uBin'] = uBin
				tExp['nBin'] = nBin
				tExp['fBinSub'] = fBinSub
				tExp['phsBin'], _  = phaseFromJD(tBin, tZer=t0)
			else:
				tReg = Table()
				tReg['tBin'] = jd
				tReg['uBin'] = dy
				tReg['fBinSub'] = fSub
				tReg['phsBin'], _ = phaseFromJD(jd, tZer=t0)


			# export this to disk
			# filBinned='v404_binSub.fits'  ### 2019-03-18 promoted to argument to the method
			if showBinned:
				if os.access(filBinned, os.R_OK):
					os.remove(filBinned)
				tExp.write(filBinned)
			else:
				if os.access(filExp, os.R_OK):
					os.remove(filExp)
				tReg.write(filExp)

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
		yPred = twoSine2019(jd, *p1, t0=t0)
		yPredLow = twoSine2019(jd, *pLow, t0=t0)

		#Generating random points throughout the best fit line
		#num_points = 1401
		#randfunc = 0.0393104 * np.sin(2.0*np.pi*jd/6.47155519  + -1.68599974)+16.5961013 + 0.15487321 * np.sin(4.0*np.pi*jd/6.47155519 + -1.68599974)*np.random.rand(num_points)

		if data == '2017':
			xGrid = np.linspace(57985., 57991., 1000)
		elif data =='2018A':
			xGrid = np.linspace(58269., 58275., 1000)
		elif data == '2018B':
			xGrid = np.linspace(58337., 58346., 1000)
		elif data == '1992':
			xGrid = np.linspace(8801., 8817., 1000)
		elif data == '1998':
			xGrid = np.linspace(51004., 51015., 1000)
		#randfunc = twoSine(pLow,xGrid)

		# Plotting the v404 data and the best fit line

		##print tBounds
		#	print np.shape(jd)
		#print np.shape(chu#nks)
		#return

		# 2018-04-07 WIC - modified Austin's plotting stanza to show
		# the phase curve.

		# create a few phase bins
		phaseLow, _  = phaseFromJD(tLow, tZer=t0)
		phaseGrid, _  = phaseFromJD(xGrid, tZer=t0)
		#per=p1[2]) # 2019-02-01 CHECK
		#print("ONe last debug:", p1)

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


		if compareEllipsoidals:
	 		p1_17 = np.loadtxt("/Users/amblevin/Desktop/p12017.txt")
	 		pLow_17 = np.loadtxt("/Users/amblevin/Desktop/pLow2017.txt")
	 		p1_18A = np.loadtxt("/Users/amblevin/Desktop/p12018A.txt")
	 		pLow_18A = np.loadtxt("/Users/amblevin/Desktop/pLow2018A.txt")
	 		p1_18B = np.loadtxt("/Users/amblevin/Desktop/p12018B.txt") 
	 		pLow_18B = np.loadtxt("/Users/amblevin/Desktop/pLow2018B.txt") # Currently incorrect ellipsoidal
	 		p1_92 = np.loadtxt("/Users/amblevin/Desktop/p11992.txt")
	 		pLow_92 = np.loadtxt("/Users/amblevin/Desktop/pLow1992.txt")
	 		p1_98 = np.loadtxt("/Users/amblevin/Desktop/p11998.txt")
	 		pLow_98 = np.loadtxt("/Users/amblevin/Desktop/pLow1998.txt") # Currently incorrect ellipsoidal

	 		Y17e = twoSine2019(xGrid[lG], *p1_17, t0=57983.8468)
	 		Y17Le = twoSine2019(xGrid[lG], *pLow_17, t0=57983.8468)
	 		Y18Ae = twoSine2019(xGrid[lG], *p1_18A, t0=58268.5884)
	 		Y18ALe = twoSine2019(xGrid[lG], *pLow_18A, t0=58268.5884)
	 		Y18Be = twoSine2019(xGrid[lG], *p1_18B, t0=58333.3024)
	 		Y18BLe = twoSine2019(xGrid[lG], *pLow_18B, t0=58333.3024)
	 		Y92e = twoSine2019(xGrid[lG], *p1_92, t0=48813.873)
	 		Y92Le = twoSine2019(xGrid[lG], *pLow_92, t0=48813.873)
	 		Y98e = twoSine2019(xGrid[lG], *p1_98, t0=48813.873)
	 		Y98Le = twoSine2019(xGrid[lG], *pLow_98, t0=48813.873)

	 		Y17 = Y17e
	 		Y17L = Y17Le
	 		Y18A = Y18Ae
	 		Y18AL = Y18ALe
	 		Y18B = Y18Be
	 		Y18BL = Y18BLe
	 		Y92 = Y92e
	 		Y92L = Y92Le
	 		Y98 = Y98e
	 		Y98L = Y98Le
	 		XC = tGrid[lG]

	 		if plot2Phases:
	 			XC, Y17 = plotAnotherPhase(phase=tGrid[lG], y=Y17e)
	 			XC, Y17L = plotAnotherPhase(phase=tGrid[lG], y=Y17Le)
	 			XC, Y18A = plotAnotherPhase(phase=tGrid[lG], y=Y18Ae)
	 			XC, Y18AL = plotAnotherPhase(phase=tGrid[lG], y=Y18ALe)
	 			XC, Y18B = plotAnotherPhase(phase=tGrid[lG], y=Y18Be)
	 			XC, Y18BL = plotAnotherPhase(phase=tGrid[lG], y=Y18BLe)
	 			XC, Y92 = plotAnotherPhase(phase=tGrid[lG], y=Y92e)
	 			XC, Y92L = plotAnotherPhase(phase=tGrid[lG], y=Y92Le)
	 			XC, Y98 = plotAnotherPhase(phase=tGrid[lG], y=Y98e)
	 			XC, Y98L = plotAnotherPhase(phase=tGrid[lG], y=Y98Le)

	 		# Plot of Median ellipsoidals
	 			
	 		fig10 = plt.figure(10)
	 		fig10.clf()
	 		ax10 = fig10.add_subplot(111)
	 		ax10.plot(XC, Y92, c='blue', label='1992')
	 		ax10.plot(XC, Y98, c='violet', label='1998')
	 		ax10.plot(XC, Y17, c='red', label='2017')
	 		ax10.plot(XC, Y18A, c='orange', label='2018A')
	 		ax10.plot(XC, Y18B, c='yellow', label='2018B (suspect)')
	 		plt.legend()
	 		plt.title('Median Ellipsoidals')
	 		yLims = np.copy(plt.gca().get_ylim())
	 		plt.ylim([yLims[1], yLims[0]])
	 		plt.xlabel('Phase (6.4714 d period)')
	 		plt.ylabel('R (mag)')

	 		# Plot of Lower Envelope Ellipsoidals (obviously incorrect ellipsoidals are not shown)

	 		fig17 = plt.figure(17)
	 		fig17.clf()
	 		ax17 = fig17.add_subplot(111)
	 		ax17.plot(XC, Y92L, c='blue', label='1992')
	 		ax17.plot(XC, Y17L, c='red', label='2017')
	 		ax17.plot(XC, Y18AL, c='orange', label='2018A')
	 		plt.legend()
	 		plt.title('Lower Envelope Ellipsoidals')
	 		yLims = np.copy(plt.gca().get_ylim())
	 		plt.ylim([yLims[1], yLims[0]])
	 		plt.xlabel('Phase (6.4714 d period)')
	 		plt.ylabel('R (mag)')

	 		# Plot of Median Ellipsoidals in the stlye of Z04 (Once 1998 lower envelope is fixed, use those)

	 		fig18 = plt.figure(18)
	 		fig18.clf()
	 		ax18 = fig18.add_subplot(111)
	 		ax18.plot(XC, Y92, c='k', label='1992')
	 		ax18.plot(XC, Y98, c='k', ls='--', label='1998')
	 		plt.legend()
	 		plt.xlabel('Phase (6.4714 d period)')
	 		plt.ylabel('R (mag)')
	 		plt.ylim([17., 16.4])

		plt.figure(1)
		plt.clf()
		#plt.scatter(jd, mag, alpha=0.5, color='darkmagenta')
		dum = plt.scatter(tSho, mag, \
					  alpha=1., c=cScatt, s=16, \
					  cmap='inferno', zorder=25, \
					  edgecolor='0.4')
		plt.plot(tLo, yLow, 'ko', ms=7, zorder=25)
		elLowX = tGrid[lG]
		lowEllip = twoSine2019(xGrid[lG], *pLow, t0=t0)
		elLowY = lowEllip
		if plot2Phases:
			if not showPhase:
				print "showPhase must be set to True to plot ellipsoidals over 2 phases."
				return
			else:
				elLowX, elLowY = plotAnotherPhase(phase=tGrid[lG], y=lowEllip)

		plt.plot(elLowX, elLowY, c='k')
		medEllip = twoSine2019(xGrid[lG], *p1, t0=t0)
		elMidX = tGrid[lG]
		elMidY = medEllip
		if plot2Phases:
			elMidX, elMidY = plotAnotherPhase(phase=tGrid[lG], y=medEllip)
		plt.plot(elMidX, elMidY, c='g', ls='--')
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
			#plt.xlabel('Phase (%.6f d period)' % (p1[2]))			
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
			if plot2Phases:
				plt.xlim(0., 2.)
			else:
				plt.xlim(0., 1.)
		
		# 2018-04-07 - save the plot as png, naming convention after
		# the figure numbers in the routine.
		#plt.savefig('v404_fig1.png')

		##ax3.show(block=False)

	plt.figure(0)
	plt.clf()
	if data == '2018B':
		dum0 = plt.scatter(jd, mag,\
						  alpha=1., c=tbl['telescope'], s=16, \
						  cmap='RdYlGn_r', zorder=25, \
						  edgecolor='0.4')
		plt.colorbar(dum0)
		if errorbars:
				plt.errorbar(jd, mag, yerr=dy, fmt='o', \
						     ms=1, ecolor='0.3', alpha=0.5, \
						     zorder=10)
	else:
		dum0 = plt.scatter(jd, mag,\
						  alpha=1., c=flagColor, s=16, \
						  cmap='inferno', zorder=25, \
						  edgecolor='0.4')
		if errorbars:
			plt.errorbar(jd, mag, yerr=dy, fmt='o', \
						     ms=1, ecolor='0.3', alpha=0.5, \
						     zorder=10)

		# let's try overplotting the model
		#if plotEllipsoidal:
		#	plt.plot(jd, twoSine2019(jd, *pLow), 'k-')
		#	plt.xlim(8804.4, 8804.8)



	plt.show(block=False)
	
	yLims = np.copy(plt.gca().get_ylim())

	if useMags:
		plt.ylim([yLims[1], yLims[0]])	

	# Generate Data along a single sine-wave and plot it

	if genOS:
		if not plotEllipsoidal:
			print "WARNING: Generation will not work unless plotEllipsoidal=True"

		# set up time array to use for generations
		tGen = np.copy(jd)
		dy1 = np.copy(dy)

		if moreTimes:
			tRange = np.max(jd) - np.min(jd)
			tGen = np.hstack(( jd, jd+tRange+1.0))
			dy1 = np.hstack(( dy, dy ))

		parsTrue = [0.1, -4.0, 16.6]

		xDum = oneSine2019(tGen, *parsTrue)

		print "True Values: ", parsTrue
		#xDum = (tDum + mag)
		yDum = xDum + np.random.normal(size=np.size(xDum))*dy1
		#np.savetxt('yDum1.txt', yDum)

		# p1s, success1s = \
		#     	optimize.leastsq(osEF, pGuess[:], args=(tGen, yDum), \
		# 			     maxfev=int(1e6), ftol=1e-10)
		# predicted at the measurement dates
		#yPredict = oneSine(p1s, jd)

		# 2018-12-07 testing curve_fit on oneSine!

		#boundsOne = ([0.09, -4.01, 6.4713, 16.3, -np.inf], [0.11, -3.99, 6.4715, 16.9, np.inf])

		#print("Double-checking oneSine Guess:", pGuess)

		#popt1, pcov1 = optimize.curve_fit(oneSine4fit, xDum, yDum, bounds=boundsOne, p0=p0, method='trf', check_finite=True, maxfev=1e9)
		#print "popt1: ", popt1

		# 2019-01-11 Testing initial guess for multiple trials; constraining PHI

		aGuess = 43543058.
		nSets = 16
		phiGuess = np.linspace(-8.0, 8.0, num=nSets)
		oGuess = 16.
		parsFound = np.zeros((nSets, np.size(parsTrue)))

		iPlot = 0 #0 plots everything; 3 doesn't plot the first 3 plots[0,1, and 2], etc.

		# Generate fit using input guesses

		for iSet in range(nSets):
			p0OS = [aGuess, phiGuess[iSet], oGuess]

			boundsOS = ([0, -np.inf, -np.inf], [np.inf, np.inf, np.inf])

			try:
				paramsOS, pcovOS = optimize.curve_fit(oneSine2019, tGen, yDum, p0=p0OS, method='trf',\
					bounds=boundsOS, sigma=dy1, absolute_sigma=False)

				# slot the best-fit parameters into the results array
				parsFound[iSet] = np.copy(paramsOS)
			except:
				print("WARN- failed fit for trial %i, %.2f" %(iSet, p0OS[1]))
				iPlot += 1
				continue

			# fine-grained version for plotting
			tFine = np.linspace(np.min(tGen), np.max(tGen), endpoint=True, num=np.size(tGen))
			yFine = oneSine2019(tFine, *paramsOS)
			pFine1, _ = phaseFromJD(tFine, tZer=t0)
			ll2 = np.argsort(pFine1)

			# plot by JD or phase?

			tSho1 = np.copy(tGen)
			tGrid = np.copy(xGrid)
			phase1, u_phs1 = phaseFromJD(tGen, tZer=t0)
			#print np.size(phase2)
			#print np.size(tGen2)
			if showPhase:
				tSho1 = np.copy(phase1)
				tGrid = np.copy(phaseGrid)
			lG1 = np.argsort(tGrid)
				#print lG

			# Old Plots

			# plt.figure(12)
			# plt.clf()
			# plt.errorbar(tGen, yDum, dy1, ls='none', marker='o', color='b', zorder=10, alpha=0.5)
			# plt.plot(tGen, oneSine4fit(tGen, *popt1), color='r', ls='--')
			# plt.plot(tGen, oneSine4fit(tGen, *p0), color='g', ls=':')
			# #plt.plot(jd, (1.09084866e-01*np.sin((2.0*np.pi*jd/7.65892886e+00)+8.90158004e+03)+1.69284759e+01), c='orange')

			# yLims = np.copy(plt.gca().get_ylim())

			# if useMags:
			# 	plt.ylim([yLims[1], yLims[0]])

			# New Plots

			if iSet == iPlot:
				fig12 = plt.figure(12)
				fig12.clf()
				plt.scatter(tSho1, yDum, color='b', s=10, zorder=10)
				plt.errorbar(tSho1, yDum, yerr=dy1, ls='none', color='b', alpha=0.3, zorder=11)

				plt.title('phi guess: %.2f' % (p0OS[1]) )
				if not showPhase:
					plt.xlabel('time(jd - 2 400 000)')
				else:
					plt.xlabel('phase')
				plt.ylabel('mag')

			plt.plot(tGrid[lG1], oneSine2019(xGrid[lG1], *paramsOS))#, label='fit: a=%5.3f, phi=%5.3f, offset=%5.3f' % tuple(paramsOS), zorder=1)
		#plt.plot(tFine, yFine, label='attempt to plot fine-grained curve', color='k', zorder=1)
		plt.legend()
		plt.show()


		# Plot only the trials that didn't fail
		phiOKOS = np.abs(parsFound[:,0]) > 1e-3

		# Plot to show how brst-fit parameters depend on choice of initial-guess phi

		fig13 = plt.figure(13)
		fig13.clf()
		ax13a = fig13.add_subplot(311)
		dum13a = ax13a.scatter(phiGuess[phiOKOS], parsFound[phiOKOS,0], color='r', marker='v', zorder=2, edgecolor='0.5')
		ax13b = fig13.add_subplot(312)
		dum13b = ax13b.scatter(phiGuess[phiOKOS], parsFound[phiOKOS,1], color='r', marker='v', zorder=2, edgecolor='0.5')
		ax13c = fig13.add_subplot(313)
		dum13c = ax13c.scatter(phiGuess[phiOKOS], parsFound[phiOKOS,2], color='r', marker='v', zorder=2, edgecolor='0.5')

		ax13a.set_ylabel('Final a')
		ax13b.set_ylabel('Final phi')
		ax13c.set_ylabel('Final Offset')

		ax13a.set_ylim(parsTrue[0]-0.2, parsTrue[0]+0.2)
		ax13b.set_ylim(parsTrue[1]-4, parsTrue[1]+12)
		ax13c.set_ylim(parsTrue[2]-0.2, parsTrue[2]+0.2)

		# Overplotting the true value

		ax13a.plot(phiGuess, np.repeat(parsTrue[0], np.size(phiGuess)), color='g', lw='1', zorder=1)
		ax13b.plot(phiGuess, np.repeat(parsTrue[1], np.size(phiGuess)), color='g', lw='1', zorder=1)
		ax13c.plot(phiGuess, np.repeat(parsTrue[2], np.size(phiGuess)), color='g', lw='1', zorder=1)

		for ax in [ax13a, ax13b, ax13c]:
			ax.set_xlabel('Initial guess phi')

		#print "p0: ", p0
		#print "p1s: ", p1s
		#print "success1s: ", success1s
		#print "size(xDum) ", np.size(xDum)

	if genTS:
		if not plotEllipsoidal:
			print "WARNING: Generation will not work unless plotEllipsoidal=True"

		# set up time array to use for generations
		tGen2 = np.copy(jd)

		#2019-02-01 let's subtract off the smallest jd value to remove extreme sensitivity to the period
		# tGen2 -= np.min(tGen2)
		dy2 = np.copy(dy)

		if moreTimes:
			tRange2 = np.max(jd) - np.min(jd)
			tGen2 = np.hstack(( tGen2, tGen2+tRange2+1.0))
			dy2 = np.hstack(( dy, dy ))

		parsTrue2 = [amp1, -2., 16.73, amp2]

		xDum2 = twoSine2019(tGen2, *parsTrue2, t0=t0)

		print "True values: ", parsTrue2
		#xDum = (tDum + mag)
		yDum2 = xDum2 + np.random.normal(size=np.size(xDum2))*dy2
		#np.savetxt('yDum2.txt', yDum2)

		# 2018-12-03 let's try specifying bounds on the parameters
		#boundsTwo = ([0.,2.],[0., 2.0*np.pi],[6.46, 6.48],[16.0, 17.0],[0., 2.])
		#boundsTwo = ([0.,0.,6.46, 16.0, 0.], \
		#	[2., 2.*np.pi, 6.48, 17.0, 2.0])

		#boundsTwo = ([0, -2.0*np.pi, 6.4713, -np.inf, -np.inf], [0.5, 2.0*np.pi, 6.4715, np.inf, np.inf])
#			0.,0.,6.46, 16.0, 0.], \
#			[2., 2.*np.pi, 6.48, 17.0, 2.0])


		#print("DOUBLE-CHECKING THE GUESS:", pGuess)
			# NOTE - inserting the actual parameters seems to not reproduce the actual data...

		#popt2, pcov2 = optimize.curve_fit(twoSine4fit, xDum2, yDum2, bounds=boundsTwo, p0=p0, method='trf')# Guess[:])
		#print "popt2: ", popt2
		#print "shape(popt): ", np.shape(popt)
		#print "pcov: ", pcov
		#print "shape(pcov): ", np.shape(pcov)

		# 2018-12-03 let's try invoking a better minimizer that can handle bounds...
		#p2s, success2s = \
		#    	optimize.leastsq(errFunc, pGuess[:], args=(tGen2, yDum2), \
		#			     maxfev=int(1e6), ftol=1e-10)
		# let's keep the old version in a comment
		#p2s, success2s = \
		#bigfit = \
		#    	optimize.least_squares(errFunc, pGuess[:], args=(tGen2, yDum2), \
		#    		bounds=boundsTwo,\
		#			     max_nfev=int(1e6), ftol=1.0e-10)

		#let's lift out what we want
		# p2s = np.copy(returnedStuff)
		# success2s = returnedStuff.status
		# print("FITTED", np.shape(p2s))
		# #print(dum)
		# #return
		# # predicted at the measurement dates
		# yPred2 = twoSine(p2s, jd)


		# 2019-01-14 Testing initial twoSine guess for multiple trials; constraining phi

		a1Guess = 54.
		nSets2 = 8
		phiGuess2 = np.linspace(-4.0, 4.0, num=nSets2)
		#phiGuess2 = -3.
		oGuess2 = 16.
		a2Guess = -345.54
		parsFound2 = np.zeros((nSets2, np.size(parsTrue2)))

		iPlot2 = 0 #0 plots everything; 3 doesn't plot the first 3 plots [0, 1, and 2], etc.


		# DEBUG

		# p0TS = [a1Guess, phiGuess2, oGuess2, a2Guess]
		# boundsTS = ([-np.inf, -np.inf, -np.inf, -np.inf], [np.inf, np.inf, np.inf, np.inf])
		# paramsTS, pcovTS = optimize.curve_fit(twoSine2019, tGen2, yDum2, p0=p0TS, method='trf',\
		# 	bounds=boundsTS, sigma=dy2, absolute_sigma=False)

		# fig14 = plt.figure(14)
		# fig14.clf()
		# plt.scatter(tGen2, yDum2, color='b', s=10, zorder=10)
		# plt.errorbar(tGen2, yDum2, yerr=dy2, ls='none', color='b', alpha=0.3, zorder=11)

		# plt.title('phi guess: %.2f' % (p0TS[1]) )
		# plt.xlabel('time(jd - 2 400 000)')
		# plt.ylabel ('mag')

		# Generate fit using input guesses

		for iSet2 in range(nSets2):
			p0TS = [a1Guess, phiGuess2[iSet2], oGuess2, a2Guess]

			boundsTS = ([-np.inf, -np.inf, -np.inf, -np.inf], [np.inf, np.inf, np.inf, np.inf])

			pGen2, _ = phaseFromJD(tGen2, tZer=t0)

			# 2019-02-18 note to self: once this works, can put it back into the try/except clause
			#try:
			paramsTS, pcovTS = optimize.curve_fit(twoSinePhase, pGen2, yDum2, p0=p0TS, method='trf',\
					bounds=boundsTS, sigma=dy2, absolute_sigma=True)
				# QQQ 2019-02-18

			try:
				# slot the best-fit parameters into the results array
				parsFound2[iSet2] = np.copy(paramsTS)
			except:
				print("WARN- failed fit for trial %i, %.2f" %(iSet2, p0TS[1]))
				iPlot2 += 1
				continue

			if parsFound2[iSet2][3] < 0:
				parsFound2[iSet2][3] *= -1
				parsFound2[iSet2][0] *= -1
				parsFound2[iSet2][1] -= np.pi
			if parsFound2[iSet2][1] > np.pi:
				parsFound2[iSet2][1] -= 2 * np.pi
			if parsFound2[iSet2][1] < -np.pi:
				parsFound2[iSet2][1] += 2 * np.pi

			# fine-grained version for plotting
			tFine2 = np.linspace(np.min(tGen2), np.max(tGen2), endpoint=True, num=1000)
			pFine2, _ = phaseFromJD(tFine2, tZer=t0)

			yFine2 = twoSine2019(tFine2, *paramsTS, t0=t0)
			ll2 = np.argsort(pFine2)

			# plot by JD or phase?

			tSho2 = np.copy(tGen2)
			tGrid = np.copy(xGrid)
			phase2, u_phs2 = phaseFromJD(tGen2, tZer=t0)
			#print np.size(phase2)
			#print np.size(tGen2)
			if showPhase:
				tSho2 = np.copy(phase2)
				tGrid = np.copy(phaseGrid)
			lG2 = np.argsort(tGrid)
				#print lG
				#print tGrid

			# old plots

			# plt.figure(14)
			# plt.clf()
			# plt.errorbar(tGen2, yDum2, dy2, ls='none', marker='o', color='b', zorder=10, alpha=0.5)
			# plt.plot(tGen2, twoSine4fit(tGen2, *popt2), color='r', ls='--')
			# plt.plot(tGen2, twoSine4fit(tGen2, *p0), color='g', ls=':')

			# yLims = np.copy(plt.gca().get_ylim())

			# if useMags:
			# 	plt.ylim([yLims[1], yLims[0]])

			#print "p0: ", p0
			#print "p2s: ", p2s
			#print "success2s: ", success2s
			#print "xDum2: ", xDum2
			#print "size(xDum2) ", np.size(xDum2)

			if iSet2 == iPlot2:
				#print "x: ", np.size(tSho2)
				#print "y: ", np.size(yDum2)

				ellip = twoSine2019(xGrid[lG2], *paramsTS, t0=t0)
				elX = tGrid[lG2]
				elY = ellip

				#print("tGrid: ", np.size(tGrid[lG2]))
				#print("xGrid: ", np.size(xGrid[lG2]))

				if plot2Phases:
						elX, elY = plotAnotherPhase(phase=tGrid[lG2], y=ellip)

				fig14 = plt.figure(14)
				fig14.clf()
				plt.scatter(tSho2, yDum2, color='b', s=10, zorder=10)
				plt.errorbar(tSho2, yDum2, yerr=dy2, ls='none', color='b', alpha=0.3, zorder=11)

				plt.title('phi guess: %.2f' % (p0TS[1]) )
				if not showPhase:
					plt.xlabel('time(jd - 2 400 000)')
				else:
					plt.xlabel('phase')
				plt.ylabel ('mag')

			# 2019-02-01 WIC - fudge to force the same parameters

			plt.plot(elX, elY)#, label='fit: a1=%5.3f, phi=%5.3f, offset=%5.3f, a2=%5.3f' % tuple(paramsTS), zorder=1)

		# 2019-01-02 - make a copy of p1 so that we can put it through twosine:
		p1For2019 = np.array([p1[0], p1[1], p1[2], p1[3]])
		p1_X = tGrid[lG2]
		p1_X2 = tGrid[lG2]
		#oldp1 = twoSine(p1, xGrid[lG2])
		newerp1 = twoSine2019(xGrid[lG2], *p1For2019, t0=t0)
		p1_2019_Y = newerp1
		#p1_Y = oldp1
		if plot2Phases:
			p1_X, p1_2019_Y = plotAnotherPhase(phase=tGrid[lG2], y=newerp1)
		#	p1_X2, p1_Y = plotAnotherPhase(phase=tGrid[lG2], y=oldp1)
		plt.plot(p1_X, p1_2019_Y, 'g-.', label=("twoSine2019(p1)"))
		#plt.plot(p1_X2, p1_Y, 'g--', lw=2, label="twoSine(p1)")


		plt.legend()
		plt.ylim(16.95, 16.45)
		plt.show()

		print("DEBUG: 'TRUE:' twoSine2019:", parsTrue2)
		print("DEBUG: 'FITTED:' twoSine2019:", paramsTS)
		print("DEBUG: 'twoSine' green dashed:", p1)
		print("DEBUG - time range: %.2f, %.2f" %( np.min(tGen2), np.max(tGen2)))


		# Plot only the trials that didn't fail
		phiOKTS = np.abs(parsFound2[:,0]) > 1e-3

		#Plot to show how best-fit parameters depend on choice of initial guess phi

		fig15 = plt.figure(15)
		fig15.clf()
		ax15a = fig15.add_subplot(221)
		dum15a = ax15a.scatter(phiGuess2[phiOKTS], parsFound2[phiOKTS,0], color='r', marker='v', zorder=2, edgecolor='0.5')
		ax15b = fig15.add_subplot(222)
		dum15b = ax15b.scatter(phiGuess2[phiOKTS], parsFound2[phiOKTS,1], color='r', marker='v', zorder=2, edgecolor='0.5')
		ax15c = fig15.add_subplot(223)
		dum15c = ax15c.scatter(phiGuess2[phiOKTS], parsFound2[phiOKTS,2], color='r', marker='v', zorder=1, edgecolor='0.5')
		ax15d = fig15.add_subplot(224)
		dum15d = ax15d.scatter(phiGuess2[phiOKTS], parsFound2[phiOKTS,3], color='r', marker='v', zorder=1, edgecolor='0.5')


		ax15a.set_ylabel('Final a1')
		ax15b.set_ylabel('Final phi')
		ax15c.set_ylabel('Final offset')
		ax15d.set_ylabel('Final a2')

		ax15a.set_ylim(parsTrue2[0]-0.2, parsTrue2[0]+0.2)
		ax15b.set_ylim(parsTrue2[1]-2, parsTrue2[1]+6)
		ax15c.set_ylim(parsTrue2[2]-0.2, parsTrue2[2]+0.2)
		ax15d.set_ylim(parsTrue2[3]-0.2, parsTrue2[3]+0.2)

		# Overplotting the true value

		ax15a.plot(phiGuess2, np.repeat(parsTrue2[0], np.size(phiGuess2)), color='g', lw='1', zorder=1)
		ax15b.plot(phiGuess2, np.repeat(parsTrue2[1], np.size(phiGuess2)), color='g', lw='1', zorder=1)
		ax15c.plot(phiGuess2, np.repeat(parsTrue2[2], np.size(phiGuess2)), color='g', lw='1', zorder=1)
		ax15d.plot(phiGuess2, np.repeat(parsTrue2[3], np.size(phiGuess2)), color='g', lw='1', zorder=1)

		for ax2 in [ax15a, ax15b, ax15c, ax15d]:
			ax2.set_xlabel('initial guess phi')

		# Subtract the ellipsoidal and plot the subtracted data

		if plotSubtractedData:
			fSub2 = yDum2 - twoSine2019(tGen2, *paramsTS, t0=t0)

			fig16 = plt.figure(16)
			fig16.clf()
			ax16 = fig16.add_subplot(111)
			ax16.scatter(tGen2, fSub2)
			if errorbars:
				ax16.errorbar(tGen2, fSub2, yerr=dy2, fmt='o', ms=4, ecolor='0.5', alpha=0.5)
			ax16.set_xlabel('Time (jd - 2 400 000)')	


	#return tZer

def oneSine(p,x):

	"""Single sine wave"""

	return p[0] * np.sin(2.0*np.pi*x/p[2] + p[1]) + p[3]

def twoSine(p,x):

	"""'Hardcoded' double sine"""

	sineOne = p[0] * np.sin(2.0*np.pi*x/p[2] + p[1]) + p[3]
	sineTwo = p[4] * np.sin(4.0*np.pi*x/p[2] + p[1])

	return sineOne + sineTwo

def oneSine4fit(x, p0, p1, p2, p3, p4):
	"""OneSine for fit"""

	return p0 * np.sin(2.0*np.pi*x/p2 + p1) + p3

def twoSine4fit(x, p0, p1, p2, p3, p4):
	"""TwoSine for fit"""

	sineOne = p0 * np.sin(2.0*np.pi*x/p2 + p1) + p3
	sineTwo = p4 * np.sin(4.0*np.pi*x/p2 + p1)

	return sineOne + sineTwo

def oneSine2019(x, a, phi, offset):
	"""Basic single sine ellipsoidal to be fitted."""

	return a * np.sin(2.0*np.pi*x/6.4714 + phi) + offset

def twoSine2019(x, a1, phi, offset, a2, t0=48813.873):
	"""Double-amplitude ellipsoidal fit"""

	# 2019-02-01 - replaced 6.4714 with 6.471528

	# 2019-02-01 - calculate the phase in the exact same way as the data
	phase, _ = np.array(phaseFromJD(x, tZer=t0))

	result = twoSinePhase(phase, a1, phi, offset, a2)

#	#sineOne = a1 * np.sin(2.0*np.pi*x/6.471528 + phi) + offset
#	#sineTwo = a2 * np.sin(4.0*np.pi*x/6.471528 + phi)

#	sineOne = a1 * np.sin(2.0*np.pi*phase + phi) + offset
#	sineTwo = a2 * np.sin(4.0*np.pi*phase + phi)


#	return sineOne + sineTwo

	return result

def twoSinePhase(phase, a1, phi, offset, a2):

	"""2019-02-18 - return the double-sinusoid expecting phase as the independent variable"""

	sineOne = a1 * np.sin(2.0*np.pi*phase + phi) + offset
	sineTwo = a2 * np.sin(4.0*np.pi*phase + phi)

	return sineOne + sineTwo

def twoSineTime(time, a1, phi, offset, a2):

	"""2019-02-21 - return the double-sinusoid expecting time as the independent variable"""

	sineOne = a1 * np.sin(2.0*np.pi*time + phi) + offset
	sineTwo = a2 * np.sin(4.0*np.pi*time + phi)

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

# def twoSine(p,x):

# 	"""'Hardcoded' double sine"""

# 	sineOne = p[0] * np.sin(2.0*np.pi*x/p[2] + p[1]) + p[3]
# 	sineTwo = p[4] * np.sin(4.0*np.pi*x/p[2] + p[1])

# 	return sineOne + sineTwo

def errFunc(p, x, y):

    """The error function"""

    return twoSine(p, x) - y

def osEF(p, x, y):

	"""The error function for OneSine"""

	return oneSine(p, x) - y

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
		clipOutliers=True, data='2018'): #dy=np.array([])

	"""Assigns the lower envelopes to the t,y array , partitioned by chunks, using the pct'th percentile """

	# we will need average time, mag arrays at least
	nChunks = np.max(chunks) - np.min(chunks)+1
	tEnv = np.zeros(nChunks)
	yEnv = np.copy(tEnv)

	#print("CHUNK DEBUG 0:", np.shape(t), np.shape(y), np.shape(chunks))

	# now we populate:
	for iChunk in range(nChunks):

		# which points are in this chunk?
		bChunk = chunks == iChunk

		if np.sum(bChunk) < 1:
			#print("WARNING - no points in chunk %i: %.2f %.2f" % (iChunk, np.min(t), np.max(t) ))
			continue

		# find the envelope for THIS chunk...
		#tMed, yMed = findLowerValue(t[bChunk], y[bChunk], pct, useMags=useMags, clipOutliers=clipOutliers)

		#print("CHUNK DEBUG:", np.sum(bChunk))
		tMed, yMed = findLowPoints(t[bChunk], y[bChunk], data=data)

		print "tMed: ", tMed
		print "yMed: ", yMed

		tEnv[iChunk] = tMed
		yEnv[iChunk] = yMed

		print("Chunk %i: tMed=%.2f, yMed=%.2f" % (iChunk, tMed, yMed))

	#print("Chunk INFO:", tEnv)
	#print("Chunk INFO:", yEnv)


	# 2019-03-11 let's try trimming out the points with zero time:
	bGood = tEnv > 0

	return tEnv[bGood], yEnv[bGood]


def findLowPoints(t=np.array([]), y=np.array([]), restrictMedian=True, clipSigma=3., clipIters=3, data='2018'): #dy=np.array([])

	"""Finds a lower envelope using the method from Pavlenko et al. (1996)."""

	yUse = np.copy(y)
	yClip = sigma_clip(y, sigma=clipSigma, iters=clipIters)
	bKeep = ~yClip.mask
	#yUse = yUse[bKeep]

	#dyUse = dy[bKeep]

	if data == '1992':
		sig = 0.02
	else:
		sig = np.std(yUse[bKeep])

	#print "dy: ", dy[bKeep]
	print "std: ", sig

	# 2019-03-19 - let's apply both booleans separately to retain the dimensions
	bNotClipped = bKeep

	# Apply the Zurita criterion (remember larger numbers mean FAINTER!)
	bUse = ((np.max(yUse[bNotClipped])-yUse)/sig) < 2

	# now we take the subset of points not clipped AND passing the zurita criterion:
	bEnv = (bNotClipped) & (bUse)

	print np.sum(bNotClipped), np.sum(bUse), np.sum(bEnv)

	print "Number of points used in this chunk: ", np.size(yUse[bEnv]), "of ", np.size(yUse)

	#print ((np.max(yUse)-yUse)/sig)

	yAvg = np.average(yUse[bEnv])

	if restrictMedian:
		tMed = np.average(t[bUse])
	else:
		tMed = np.average(t)

	# fig19 = plt.figure(19)
	# fig19.clf()
	# ax19 = fig19.add_subplot(111)
	# ax19.scatter(t[bUse], y[bUse])
	# #ax19.errorbar(t[bUse], y[bUse], yerr=dy[bUse], fmt='o', ls='none', ms=4, alpha=0.5, ecolor='0.5')
	# yLims = np.copy(plt.gca().get_ylim())
	# plt.ylim([yLims[1], yLims[0]])


	return tMed, yAvg


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

def correctMagForContaminant(magObs=20., magContam=17.52):

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

def phaseFromJD(jdShort=np.array([]), per=6.4714, \
		tZer=48813.873, \
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

	print("tZer: ", tZer)

	return phs, u_phs

def projectEphemeris(minJD=48813.873, POrb=6.4714):

	"""Projects the Casares ephemeris forward based on the first date of an observing run."""

	NDays = minJD - 48813.873 #Casares Ephemeris is JD 48813.873
	NPOrb = NDays / POrb
	NInt = np.floor(NPOrb)
	newT0 = NInt * POrb + 48813.873

	return newT0

def plotAnotherPhase(phase=np.array([]), y=np.array([])):

	"""Enables plotting of ellipsoidal curves over two phases."""

	num = np.size(phase)
	#print("num: ", num)
	#print ("y: ", np.size(y))
	addition = np.ones(num)
	extra = phase + addition
	pNew = np.hstack((phase, extra))
	#print("pNew: ", np.size(pNew))
	yNew = np.hstack((y, y))
	#print ("yNew: ", np.size(yNew))

	return pNew, yNew

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

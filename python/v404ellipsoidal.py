#
# v404ellipsoidal.py
#

# moved the ellipsoidal modulation info across to a separate module to
# avoid confusion

import pickle
import os
import numpy as np
import matplotlib.pylab as plt
from astropy.table import Table
from scipy import optimize


# to fit the ellipsoidal modulation
from scipy.interpolate import UnivariateSpline

def fitPaperEllipsoidal(pathIn='zurita04_fig2_2003.txt'):

	"""Fits an ellipsoidal modulation to the curves screen-grabbed
	from Zurita et al."""

	if not os.access(pathIn, os.R_OK):		
		print "fitPaperEllipse FATAL - cannot read path %s" % (pathIn)
		aZer = np.array([])
		return aZer, aZer, aZer, aZer

	tCurv = Table.read(pathIn, format='ascii')
	phs = tCurv['x']
	mag = tCurv['y']

	# print np.shape(phs)

	phs = np.hstack([phs, phs + 1.0])
	mag = np.hstack([mag, mag])

	iSor = np.argsort(phs)
	phs = phs[iSor]
	mag = mag[iSor]

	# guess for the twosine
	pGuess = np.array([0.2, 1.32, 1.0, 0.05, 0.1])
	pEll, successEll = optimize.leastsq(errFunc, pGuess[:], args=(phs, mag))
#, \
#						    maxfev=int(1e8), ftol=1e-15)
	# print pGuess
	# print pEll
	# print successEll
	return pEll, successEll, phs, mag

def showPaperEllipsoidal(smo=1e-3, yStep=0.3, lineupVert=False, \
				 showOnlyDipRoots=False, \
				 reportMinima=True):

	"""Fit and show the ellipsoidal modulations from the
	literature. Set smo=0 for the spline to simply connect the
	dots. A few arguments:

	smo=1e-3 -- spline smoothing factor for ellipsoidal representation

	yStep = 0.3 -- progressive vertical step between epochs

	lineupVert -- line up the curves by the phase=0.25 peak

	showOnlyDipRoots: only highlight the locations of the minima

	reportMinima -- reports the minima difference to the screen"""

	# ensure light background
	plt.style.use('bmh')
	plt.style.use('ggplot')

	# construct the filename strings
	lTimes = '2003', '2001', '1998', '1992'

	lCurvs = ['zurita04_fig2_%s.txt' % (sTime) for sTime in lTimes]

	# reverse the list
	lCurvs.reverse()

	# add bernardini to the list
	lCurvs.append('bernardini_fig2_2015.txt')

	# now add our own dataset
	lCurvs.append('mdm_fig2_2017.txt')

	# set up an array for the fit parameters
	aPars = np.zeros((len(lCurvs), 5))

	# set up the figure
	fig1 = plt.figure(8)
	fig1.clf()
	ax1 = fig1.add_subplot(111)
	tFine = np.linspace(0., 2., 5000, endpoint=True)

	# plot symbols
	lSyms = ['o', '^', 'v', '+', 'x', 's']

	# filename stem for serialization
	DPickles = {}

	# count down through the ellipsoidal modulations
	for iCurv in range(len(lCurvs)):

		pars, succ, jd, mag = fitPaperEllipsoidal(lCurvs[iCurv])

		# we don't want those arbitrary offsets in our
		# magnitudes. So remove them. FUDGE WARNING.
		if lCurvs[iCurv].find('2017') < 1:
			mag = mag - np.median(mag)

		#print "INFO:: ", iCurv, np.median(mag)

		if np.size(pars) < 1:
			continue

		# slot these into the params array
		aPars[iCurv] = pars
		yFine = twoSine(pars, tFine)

		# Now try with the spline
		ss = UnivariateSpline(jd, mag, k=4, s=smo)

		# re-extract the date string
		sDate = lCurvs[iCurv].split('fig2_')[-1].split('.')[0]

		# pass this up to the dictionary
		DPickles[sDate] = {}
		DPickles[sDate]['spln'] = ss

		# plot symbol
		pSym = lSyms[iCurv % len(lSyms)]

		#magSho = mag - pars[3] + 0.2*iCurv
		#yFineSho = yFine - pars[3] + 0.2*iCurv

		# predicted y spline and derivative
		ySpline = ss(tFine)

		# Find the stationary points and evaluate the curve at those points
		ySplineDeriv = ss.derivative(n=1)	
		xRoots=ySplineDeriv.roots()
		lRoots = (xRoots > 0.1) & (xRoots < 1.1)
		xRoots = xRoots[lRoots]
		yRoots=ss(xRoots)

		# pass these up too
		DPickles[sDate]['xRoots'] = xRoots
		DPickles[sDate]['yRoots'] = yRoots

		# also the magnitude difference between the phase ~
		# 1.0 dip and the phase ~0.5 dip
		DPickles[sDate]['dDip'] = yRoots[3] - yRoots[1]

		# offset to add to all the magnitudes when overplotting
		yOff = yStep*iCurv - np.median(ySpline)

		if lineupVert:
			yOff = 0. - yRoots[0]

		# a few plot symbol variables. Change this to change
		# the behavior on the plot.
		lw = 1
		pCol='k'
		ls='-'
		if sDate.find('2017') > -1:
			lw = 2
			pCol='b'
			ls='-'

		if sDate.find('2015') > -1:
			pCol='0.3'
			ls='--'

		magSho = mag + yOff
		ySplineSho = ySpline + yOff
		yRoots += yOff

		dum1 = ax1.scatter(jd, magSho, zorder=10, \
					   label=sDate, marker=pSym, \
					   s=.1)

		#dum2 = ax1.plot(tFine, yFineSho, ls='-', color='0.3')
		
		dum3 = ax1.plot(tFine, ySplineSho, ls=ls, color=pCol, \
					lw=lw, zorder=15)

		# which roots do we show?
		ll = np.argsort(xRoots)
		if showOnlyDipRoots:
			ll = np.array([1,3], 'int')

		dum4 = ax1.plot(xRoots[ll], yRoots[ll], 'ks', zorder=10, ms=6)

		if sDate.find('2015') > -1:
			sDate = '2006-2015'

		if reportMinima:
			print "%s: min(0.5) - min(1.0) = %+.3f" \
			    % (sDate, yRoots[1]-yRoots[3])


		# annotate
		iAnno = np.argmax(jd)
		dumAnno = ax1.annotate(sDate, \
					       (jd[iAnno]+0.05, magSho[iAnno]), \
					       xycoords='data', ha='left', va='top', \
					       fontsize=12, color=pCol)

	# serialize the ellipsoidal curves to disk
	filEllip = 'test_ellip.pickle'
	pickle.dump(DPickles, open(filEllip, 'w'))

	# flip the axes
	yAxes = np.copy(ax1.get_ylim())
	ax1.set_ylim(yAxes[1], yAxes[0])

	# label the axes
	ax1.set_xlabel(r'Phase')
	ax1.set_ylabel(r'$\Delta R$, mag')

	# save the figure 
	fig1.savefig('zurita04_ellipsoidal.png')

def loadEllipsoidal(pathIn='test_ellip.pickle'):

	"""Loads ellipsoidal modulation characterization"""
	
	if not os.access(pathIn, os.R_OK):
		print "loadEllipsoidal WARN - cannot read path %s" \
		    % (pathIn)
		return {}

	DPickles = pickle.load(open(pathIn, 'r'))
	return DPickles

def twoSine(p,x):

	"""'Hardcoded' double sine"""

	sineOne = p[0] * np.sin(2.0*np.pi*x/p[2] + p[1]) + p[3]
	sineTwo = p[4] * np.sin(4.0*np.pi*x/p[2] + p[1])

	return sineOne + sineTwo

def errFunc(p, x, y):

    """The error function"""

    return y - twoSine(p, x)

def makeEllipLayers(dirOut='ellipPNGs', filEllip='test_ellip.pickle', \
			    showDips=True, lightbg=True, \
			    hideHorizAxis=False):

	"""Makes pngs with transparent background for overplotting in
	powerpoint or slides. Frames are put into their own
	subdirectory. Some arguments:

	showDips -- label the phase 0.5 and 1.0 points

	lightbg -- make frames with light background?

	hideHorizAxis -- don't draw the horizontal axis line."""

	# It's probably better for presentations and posters to allow
	# the author to do this by hand (e.g. the time interval
	# between datasets is >> the time interval along the
	# ellipsoidal modulation). So:

	if lightbg:
		plt.style.use('seaborn-white')
		#plt.style.use('seaborn-ticks')
		lineColorAll='k'
		sLight = 'lightBg'
	else:
		plt.style.use('dark_background')
		lineColorAll='w'
		sLight = 'darkBg'

	# load the set of ellipsoidal modulation representations
	if not os.access(filEllip, os.R_OK):
		print "makeEllipLayers WARN - cannot read input file %s" \
		    % (filEllip)
		return

	DD = pickle.load(open(filEllip, 'r'))
	
	# generate a phase array we'll use for all the curves
	phs = np.linspace(0., 2., 200)

	# standardize the vertical range
	xRange = [0., 2.]
	#yRange = [17.0, 16.4]
	yRange = [0.25, -0.3]

	# set up the output path
	dirOut = '%s_%s' % (dirOut, sLight)
	if not os.access(dirOut, os.R_OK):
		os.mkdir(dirOut)

	print "makeEllipLayers INFO - will put images into subdirectory %s" \
	    % (dirOut)

	# now loop through the characterizations we've made
	for sYear in DD.keys():
		func = DD[sYear]['spln']
		mag = func(phs)

		# we ensure that we're dealing with delta-mag from the
		# median in each case.
		yMed = np.median(mag)
		mag -= yMed

		# use our color scheme of black (or white) for
		# pre-2017, and blue (or yellow) for the 2017
		if sYear.find('2017') < 0:
			lineColor = lineColorAll
		else:
			if lightbg:
				lineColor='b'
			else:
				lineColor='y'

		fig1=plt.figure(1, figsize=(10,10))
		fig1.clf()
		ax1=fig1.add_subplot(111)
		dum = ax1.plot(phs, mag, lw=3, color=lineColor)
		ax1.set_xlim(xRange)
		ax1.set_ylim(yRange)

		lRoots = np.array([1,3], 'int')
		if showDips:

			# (re-) find the roots
			xRoots = func.derivative(n=1).roots()
			bRoots = xRoots > 0.1
			xRoots = xRoots[bRoots]

			ll = np.argsort(xRoots)
			xRoots = xRoots[ll]
			yRoots = func(xRoots) - yMed

			# we'll do one symbol per root

			ax1.scatter(xRoots[1], yRoots[1], \
					    marker='s', s=100, \
					    zorder=25, color=lineColor)

			ax1.scatter(xRoots[3], yRoots[3], \
					    marker='p', s=100, \
					    zorder=20, color=lineColor)
			ax1.scatter(xRoots[3], yRoots[3], \
					    marker='p', s=49, \
					    color='w', \
					    zorder=25)


		# some axis carpentry
		ax1.spines['top'].set_visible(False)
		ax1.spines['right'].set_visible(False)
		
		# give the range in the spines
		yLo = np.min(mag)
		yHi = np.max(mag)
		
		ax1.spines['left'].set_bounds(yLo, yHi)

		ax1.set_ylabel(r'$\Delta R$', fontsize=16, rotation=00)
		plt.yticks(fontsize=14)
		plt.xticks(fontsize=14)
		ax1.set_xlabel('Phase', fontsize=16)

		if not hideHorizAxis:
			#ax1.set_xlabel('Phase', fontsize=16)
			ax1.spines['bottom'].set_visible(True)
			ax1.spines['bottom'].set_bounds(0., 2.)
			ax1.xaxis.set_ticks_position('bottom')
		else:
			ax1.spines['bottom'].set_visible(False)

		# set up the output path
		pathOut = '%s/ellipFrame_%s_%s.png' % (dirOut, sYear, sLight)
		fig1.savefig(pathOut, transparent=True)

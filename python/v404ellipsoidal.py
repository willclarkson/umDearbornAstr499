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

def showPaperEllipsoidal(smo=1e-3):

	"""Fit and show the ellipsoidal modulations from the
	literature. Set smo=0 for the spline to simply connect the
	dots."""

	# construct the filename strings
	lTimes = '2003', '2001', '1998', '1992'

	lCurvs = ['zurita04_fig2_%s.txt' % (sTime) for sTime in lTimes]
	
	# add bernardini to the list
	lCurvs.append('bernardini_fig2_2015.txt')

	# set up an array for the fit parameters
	aPars = np.zeros((len(lCurvs), 5))

	# set up the figure
	fig1 = plt.figure(1)
	fig1.clf()
	ax1 = fig1.add_subplot(111)
	tFine = np.linspace(0., 2., 5000, endpoint=True)

	# plot symbols
	lSyms = ['o', '^', 'v', '+', 'x', 's']

	# filename stem for serialization
	DPickles = {}

	for iCurv in range(len(lCurvs)):
		pars, succ, jd, mag = fitPaperEllipsoidal(lCurvs[iCurv])

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

		# offset to add to all the magnitudes when overplotting
		yOff = 0.3*iCurv - np.median(ySpline)
		
		magSho = mag + yOff
		ySplineSho = ySpline + yOff
		yRoots += yOff

		dum1 = ax1.scatter(jd, magSho, zorder=10, label=sDate, marker=pSym)
		#dum2 = ax1.plot(tFine, yFineSho, ls='-', color='0.3')
		
		dum3 = ax1.plot(tFine, ySplineSho, ls='-', color='0.3')

		dum4 = ax1.plot(xRoots, yRoots, 'ks', zorder=10, ms=6)

		if sDate.find('2015') > -1:
			sDate = '2006-2015'

		# annotate
		iAnno = np.argmax(jd)
		dumAnno = ax1.annotate(sDate, \
					       (jd[iAnno]+0.05, magSho[iAnno]), \
					       xycoords='data', ha='left', va='top', \
					       fontsize=12)

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


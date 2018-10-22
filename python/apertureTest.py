#
#apertureTest.py
#
#2018-10-22 AMB
#

import matplotlib.pylab as plt
plt.ion()
plt.style.use('ggplot')
from astropy.table import Table, Column
import numpy as np

#loadTable def first so tables can be read next

def loadTable(fileIn='/Users/amblevin/Desktop/apT_2017.txt'):
	"""Converts AstroImageJ (txt) tables to an output AstroPy table in the following format:
	JD/Rel. Flux/ Rel. Flux Error/ FWHM(V404)/ FWHM (C4)"""

	aDum = np.genfromtxt(fileIn, unpack=False)

	tPho = Table()

	nCols = np.shape(aDum)[-1]
	time = aDum[:,0]
	rel_flux = aDum[:,1]
	rel_flux_err = aDum[:,2]
	widthT1 = aDum[:,3]
	widthC2 = aDum[:,4]

	tPho['JD - 2 400 000'] = Column(time)
	tPho['Relative Flux'] = Column(rel_flux)
	tPho['Relative Flux Error'] = Column(rel_flux_err)
	tPho['FWHM (T1)'] = Column(widthT1)
	tPho['FWHM (C2)'] = Column(widthC2)

	return tPho



# Reading the tables

try:
	t17 = loadTable(fileIn='/Users/amblevin/Desktop/apT_2017.txt')
	t18A = loadTable(fileIn='/Users/amblevin/Desktop/apT_2018A.txt')
	t18B = loadTable(fileIn='/Users/amblevin/Desktop/apT_2018B.txt')
except UnboundLocalError:
	print "One or more files not found. Check the directories indicated in the source code."

def go(data='2017', plotFlux=False, plotC4FWHM=True, plotV404FWHM=True, plotHist=True, nBins=100):
	if data == '2017':
		tbl = t17
	elif data == '2018A':
		tbl = t18A
	elif data == '2018B':
		tbl = t18B
	else:
		print "The value for 'data' must be either 2017, 2018A, or 2018B."

	#print tbl


	#Defining Variables

	jd = tbl['JD - 2 400 000']
	flux = tbl['Relative Flux']
	dy = tbl['Relative Flux Error']
	fwhmV404 = tbl['FWHM (T1)']
	fwhmC4 = tbl['FWHM (C2)']

	#Plot of JD vs. Flux

	if plotFlux:
		fig1 = plt.figure(1)
		fig1.clf()
		plt.scatter(jd, flux, alpha=1., s=16, zorder=25)
		plt.errorbar(jd, flux, yerr=dy, fmt='o', ms=1, ecolor='0.3', alpha=0.5, zorder=10)
		plt.xlabel('JD - 2 400 000 d')
		plt.ylabel('Flux (V404) / Flux (C4)')
		#plt.show(block=False)


	#Plot of rel_flux_err vs. FWHM(C4)

	if plotC4FWHM:
		fig2 = plt.figure(2)
		fig2.clf()
		dum2 = plt.scatter(dy, fwhmC4, alpha=1., s=16, zorder=25, c=jd, cmap='hsv')
		plt.xlabel('Relative Flux Uncertainty')
		plt.ylabel('FWHM (C4)')
		plt.colorbar(dum2)


	#Plot of rel_flux_err vs. FWHM(V404)

	if plotV404FWHM:
		fig3 = plt.figure(3)
		fig3.clf()
		dum3 = plt.scatter(dy, fwhmV404, alpha=1., s=16, zorder=25, c=jd, cmap='hsv')
		plt.xlabel('Relative Flux Uncertainty')
		plt.ylabel('FWHM (V404)')
		plt.colorbar(dum3)


	#Plot of Histogram(dy)

	if plotHist:

		Range = [np.min(dy), np.max(dy)]

		fig4 = plt.figure(4)
		fig4.clf()
		ax = fig4.add_subplot(111)
		ax.hist(dy, nBins, alpha=0.5, range=Range)
		plt.xlabel('dy')
		plt.ylabel('Frequency')







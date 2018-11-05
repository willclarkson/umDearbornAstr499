#
#v404_comparison.py
#
#2018-11-05 AMB
#

import matplotlib.pylab as plt
plt.ion()
plt.style.use('ggplot')
from astropy.table import Table, Column
import numpy as np

#loadTable def first so tables can be read next

def loadTable(fileIn='/Users/amblevin/Desktop/ApC_2017.txt'):
	"""Converts AstroImageJ (txt) tables to an output AstroPy table in the following format:
	JD/Rel.Flux(V404/C4)/Rel.Flux(V404/C5)/Rel.Flux(C4/C5)/Rel.FluxError(V404/C4)/Rel.FluxError(V404/C5)/Rel.FluxError(C4/C5)/FWHM(V404)/FWHM (C4)/FWHM(C5)"""

	aDum = np.genfromtxt(fileIn, unpack=False)

	tPho = Table()

	nCols = np.shape(aDum)[-1]
	time = aDum[:,0]
	rfVC4 = aDum[:,1]
	rfVC5 = aDum[:,2]
	rfC4C5 = aDum[:,3]
	rfeVC4 = aDum[:,4]
	rfeVC5 = aDum[:,5]
	rfeC4C5 = aDum[:,6]
	fwhmV = aDum[:,7]
	fwhmC4 = aDum[:,8]
	fwhmC5 = aDum[:,9]

	tPho['JD - 2 400 000'] = Column(time)
	tPho['Relative Flux(V404/C4)'] = Column(rfVC4)
	tPho['Relative Flux(V404/C5)'] = Column(rfVC5)
	tPho['Relative Flux(C4/C5)'] = Column(rfC4C5)
	tPho['Relative Flux Error(V404/C4)'] = Column(rfeVC4)
	tPho['Relative Flux Error(V404/C5)'] = Column(rfeVC5)
	tPho['Relative Flux Error(C4/C5)'] = Column(rfeC4C5)
	tPho['FWHM (V404)'] = Column(fwhmV)
	tPho['FWHM (C4)'] = Column(fwhmC4)
	tPho['FWHM (C5)'] = Column(fwhmC5)

	return tPho


# Reading the tables

try:
	t17 = loadTable(fileIn='/Users/amblevin/Desktop/ApC_2017.txt')
	t18A = loadTable(fileIn='/Users/amblevin/Desktop/ApC_2018A.txt')
	t18B = loadTable(fileIn='/Users/amblevin/Desktop/ApC_2018B.txt')
except UnboundLocalError:
	print "One or more files not found. Check the directories indicated in the source code."

def go(data='2017', plotFlux=True, plotErrorvsFWHM=True, plotJDvsFWHM=True, errorbars=True):

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
	fluxVC4 = tbl['Relative Flux(V404/C4)']
	fluxVC5 = tbl['Relative Flux(V404/C5)']
	fluxC4C5 = tbl['Relative Flux(C4/C5)']
	dyVC4 = tbl['Relative Flux Error(V404/C4)']
	dyVC5 = tbl['Relative Flux Error(V404/C5)']
	dyC4C5 = tbl['Relative Flux Error(C4/C5)']
	V404FWHM = tbl['FWHM (V404)']
	C4FWHM = tbl['FWHM (C4)']
	C5FWHM = tbl['FWHM (C5)']

	#Plot of JD vs. Flux(V404/C4)

	if plotFlux:
		fig1 = plt.figure(1)
		fig1.clf()
		plt.scatter(jd, fluxVC4, alpha=0.4, s=16, zorder=25)
		if errorbars:
			plt.errorbar(jd, fluxVC4, yerr=dyVC4, fmt='o', ms=1, ecolor='0.3', alpha=0.2, zorder=10)
		plt.xlabel('JD - 2 400 000 d')
		plt.ylabel('Flux (V404) / Flux (C4)')

		# Plot of JD vs. Flux(V404/C5)

		fig2 = plt.figure(2)
		fig2.clf()
		plt.scatter(jd, fluxVC5, alpha=0.4, s=16, zorder=25)
		if errorbars:
			plt.errorbar(jd, fluxVC5, yerr=dyVC5, fmt='o', ms=1, ecolor='0.3', alpha=0.2, zorder=10)
		plt.xlabel('JD - 2 400 000 d')
		plt.ylabel('Flux (V404) / Flux (C5)')

		# Plot of JD vs. Flux(C4/C5)

		fig3 = plt.figure(3)
		fig3.clf()
		plt.scatter(jd, fluxC4C5, alpha=0.4, s=16, zorder=25)
		if errorbars:
			plt.errorbar(jd, fluxC4C5, yerr=dyC4C5, fmt='o', ms=1, ecolor='0.3', alpha=0.2, zorder=10)
		plt.xlabel('JD - 2 400 000 d')
		plt.ylabel('Flux (C4) / Flux (C5)')

	# Plot of Error(V404/C4) vs. FWHM(V404)

	if plotErrorvsFWHM:
		fig4 = plt.figure(4)
		fig4.clf()
		dum4 = plt.scatter(dyVC4, V404FWHM, alpha=0.4, s=16, zorder=25, c=jd, cmap='hsv')
		plt.xlabel('Relative Flux Uncertainty (V404/C4)')
		plt.ylabel('FWHM (V404)')
		plt.colorbar(dum4)

		# Plot of Error(V404/C4) vs. FWHM(C4)

		fig5 = plt.figure(5)
		fig5.clf()
		dum5 = plt.scatter(dyVC4, C4FWHM, alpha=0.4, s=16, zorder=25, c=jd, cmap='hsv')
		plt.xlabel('Relative Flux Uncertainty (V404/C4)')
		plt.ylabel('FWHM (C4)')
		plt.colorbar(dum5)

		# Plot of Error(V404/C5) vs. FWHM(V404)

		fig6 = plt.figure(6)
		fig6.clf()
		dum6 = plt.scatter(dyVC5, V404FWHM, alpha=0.4, s=16, zorder=25, c=jd, cmap='hsv')
		plt.xlabel('Relative Flux Uncertainty (V404/C5)')
		plt.ylabel('FWHM (V404)')
		plt.colorbar(dum6)

		# Plot of Error(V404/C5) vs. FWHM(C5)

		fig7 = plt.figure(7)
		fig7.clf()
		dum7 = plt.scatter(dyVC5, C5FWHM, alpha=0.4, s=16, zorder=25, c=jd, cmap='hsv')
		plt.xlabel('Relative Flux Uncertainty (V404/C5)')
		plt.ylabel('FWHM (C5)')
		plt.colorbar(dum7)

		# Plot of Error(C4/C5) vs. FWHM(C4)

		fig8 = plt.figure(8)
		fig8.clf()
		dum8 = plt.scatter(dyC4C5, C4FWHM, alpha=0.4, s=16, zorder=25, c=jd, cmap='hsv')
		plt.xlabel('Relative Flux Uncertainty (C4/C5)')
		plt.ylabel('FWHM (C4)')
		plt.colorbar(dum8)

		# Plot of Error(C4/C5) vs. FWHM(C5)

		fig9 = plt.figure(9)
		fig9.clf()
		dum9 = plt.scatter(dyC4C5, C5FWHM, alpha=0.4, s=16, zorder=25, c=jd, cmap='hsv')
		plt.xlabel('Relative Flux Uncertainty (C4/C5)')
		plt.ylabel('FWHM (C5)')
		plt.colorbar(dum9)

	# Plot of JD vs. FWHM(V404)

	if plotJDvsFWHM:
		fig10 = plt.figure(10)
		fig10.clf()
		plt.scatter(jd, V404FWHM, alpha=0.4, s=16, zorder=25)
		plt.xlabel('JD - 2 400 000 d')
		plt.ylabel('FWHM (V404)')

		# Plot of JD vs. FWHM(C4)

		fig11 = plt.figure(11)
		fig11.clf()
		plt.scatter(jd, C4FWHM, alpha=0.4, s=16, zorder=25)
		plt.xlabel('JD - 2 400 000 d')
		plt.ylabel('FWHM (C4)')

		# Plot of JD vs. FWHM(C5)

		fig12 = plt.figure(12)
		fig12.clf()
		plt.scatter(jd, C5FWHM, alpha=0.4, s=16, zorder=25)
		plt.xlabel('JD - 2 400 000 d')
		plt.ylabel('FWHM (C5)')


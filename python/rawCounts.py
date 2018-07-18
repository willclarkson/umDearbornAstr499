import matplotlib.pylab as plt
plt.ion()
plt.style.use('ggplot')
from astropy.table import Table
import numpy as np

def go(inPath='../../../Desktop/onesandtwos_a12.csv'):
	tbl = Table.read(inPath, format='ascii.csv')
	showRawCounts(tbl)



def showRawCounts(tPhot=Table(), figNam='test_rawCounts.png'):

	"""Show the raw counts"""

	if len(tPhot) < 1:
		return

	jd = tPhot['J.D.-2400000']
	countsObj = tPhot['Source-Sky_T1']
	countsCom = tPhot['Source-Sky_C2']
	flag = tPhot['Flag']


	fig1 = plt.figure(10)
	fig1.clf()
	ax1=fig1.add_subplot(211)
	ax2=fig1.add_subplot(212, sharex=ax1, sharey=ax1)

	fig1.subplots_adjust(hspace=0.05)

	pSyms = ['x', '+', 'x', 'v', 's']
	pCols = ['b', 'g', 'b', 'g', '0.4']

	# loop through the flag quality:
	for fl in range(np.max(flag)+1):
		bThis = np.abs(fl-flag) < 1.0e-2
		if np.sum(bThis) < 1:
			continue
	
		# which plot symbol?
		pSym = pSyms[fl % len(pSyms)]
		colo = pCols[fl % len(pCols)]

		dum1 = ax1.scatter(jd[bThis], countsObj[bThis], \
				   marker=pSym, \
				   color=colo, \
				   alpha=0.5, \
				   s=25)
		dum2 = ax2.scatter(jd, countsCom, \
				   marker=pSym, \
				   color=colo, \
				   alpha=0.5, \
				   s=25)

	ax2.set_xlabel('JD - 2 400 000.0 d')
	ax1.set_ylabel('Source minus sky, counts')
	ax2.set_ylabel('Source minus sky, counts')

	# standardize the y range
	ax1.set_ylim(0., 15000.)
	#ax2.set_ylim(0., 15000.)

	ax1.annotate('v404 Cyg + contaminant', (0.95, 0.95), \
		     ha='right', va='top', \
		     xycoords='axes fraction', \
		     color='k', fontsize=14)

	ax2.annotate('Comparison star "C4" (R = 16.07)', (0.95, 0.95), \
		     ha='right', va='top', \
		     xycoords='axes fraction', \
			     color='k', fontsize=14)


	#save to disk
	fig1.savefig(figNam)
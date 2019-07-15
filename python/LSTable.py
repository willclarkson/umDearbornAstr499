from astropy.table import Table
import numpy as np
import v404_test
import matplotlib.pylab as plt
plt.ion()
plt.style.use('ggplot')

def go(year='2017', nBins=35, cannedRanges=True, old=False, showTrendline=True):

	if old:
		if year == '2017':
			tbl = Table.read('/Users/amblevin/Desktop/OLD/17lsAll.txt', format='ascii')
		elif year == '2018A':
			tbl = Table.read('/Users/amblevin/Desktop/OLD/18lsAll.txt', format='ascii')
		else:
			print "'year' must be either 2017 or 2018A if old=True!"
			return
	else:
		if year == '2017':
			tbl = Table.read('/Users/amblevin/Desktop/17lsAll.txt', format='ascii')
		elif year == '2018A':
			tbl = Table.read('/Users/amblevin/Desktop/18AlsAll.txt', format='ascii')
		elif year == '2018B':
			tbl = Table.read('/Users/amblevin/Desktop/18BlsAll.txt', format='ascii')
		elif year == '1992':
			tbl = Table.read('/Users/amblevin/Desktop/92lsAll.txt', format='ascii')
		elif year == '1998':
			tbl = Table.read('/Users/amblevin/Desktop/98lsAll.txt', format='ascii')
		elif year == '1999':
			tbl = Table.read('/Users/amblevin/Desktop/99lsAll.txt', format='ascii')
		else:
			print "'year' must be '2017', '2018A', '2018B', '1992', '1998', or '1999'."
			return
	
	fAll = np.asarray(tbl['col0'])
	lsAll = np.asarray(tbl['col1'])

	# Log of LS

	logFreq = np.log10(fAll)
	logPower = np.log10(lsAll)

	# Beta value for full LS

	pwrFit = np.polyfit(logFreq, logPower, 1)

	# plot this

	f1 = plt.figure(1, figsize=(10,5))
	f1.clf()
	ax1 = f1.add_subplot(111)
	#ax1.scatter(logFreq, logPower, color='b')
	ax1.step(logFreq, logPower, color='k', label=year, lw=2)
	if showTrendline:
		ax1.plot(logFreq, np.polyval(pwrFit, logFreq), c='r', label=r"$\beta$" +": %5.2f" % pwrFit[0], lw=4, alpha=0.7)
	ax1.legend()
	ax1.set_xlabel("log(Frequency)")
	ax1.set_ylabel("log(LS Power)")

	ax1.set_title('Lomb-Scargle periodogram, %s' % (year))

	# bin LS

	freqBin, pwrBin, errBin, nPerBin = v404_test.BinData(fAll, lsAll, np.ones(len(fAll)), \
		tStart=10**-3.9, tEnd=10**-2.8, BinTime=(0.00146 / nBins))

	logFreqBin = np.log10(freqBin)
	logPwrBin = np.log10(pwrBin)
	pwrFitBin, covarBin = np.polyfit(logFreqBin, logPwrBin, 1, cov=True)
	stdResid = np.std(logPwrBin - np.polyval(pwrFitBin, logFreqBin))

	# plot binned LS

	f2 = plt.figure(2)
	f2.clf()
	ax2 = f2.add_subplot(111)
	ax2.scatter(logFreqBin, logPwrBin, color='b', label=year)
	if showTrendline:
		ax2.plot(logFreqBin, np.polyval(pwrFitBin, logFreqBin), c='r', label=r"$\beta$" +": %5.2f" % pwrFitBin[0], lw=2)

	# set legend and labels for both plots
	ax2.set_title('Binned Lomb-Scargle periodogram, %s' % (year))

	for ax in [ax1, ax2]:
		ax.legend()
		ax.set_xlabel(r"log$_{10}$(Frequency)")
		ax.set_ylabel(r"log$_{10}$(LS Power)")


	# if canned ranges set, reset the axes
	if cannedRanges:
		for ax in [ax1, ax2]:
			ax.set_xlim([-4.0, -2.6])
			ax.set_ylim([-4.0,-1.0])

	# save the figure to disk
	f1.savefig('LS_unbinned_%s.png' % (year))

	# write out the binned LS so we can use it later

	fBinned = '%s_lsBinned_%i.dat' % (year, nBins)
	tBinned = Table()
	tBinned['logFreqBin'] = logFreqBin
	tBinned['logPwrBin'] = logPwrBin
	tBinned.write(fBinned, format='ascii',overwrite=True)

from astropy.table import Table
import numpy as np
import v404_test
import matplotlib.pylab as plt
plt.ion()
plt.style.use('ggplot')

def go(year='2017', nBins=35):

	if year == '2017':
		tbl = Table.read('/Users/amblevin/Desktop/17lsAll.txt', format='ascii')
	elif year == '2018':
		tbl = Table.read('/Users/amblevin/Desktop/18lsAll.txt', format='ascii')
	else:
		print "'year' must be either 2017 or 2018!"
		return
	
	fAll = np.asarray(tbl['col0'])
	lsAll = np.asarray(tbl['col1'])

	# Log of LS

	logFreq = np.log10(fAll)
	logPower = np.log10(lsAll)

	# Beta value for full LS

	pwrFit = np.polyfit(logFreq, logPower, 1)

	# plot this

	f1 = plt.figure(1)
	f1.clf()
	ax1 = f1.add_subplot(111)
	ax1.scatter(logFreq, logPower, color='b')
	ax1.plot(logFreq, logPower, color='b')
	ax1.plot(logFreq, np.polyval(pwrFit, logFreq), c='r', label=r"$\beta$" +": %5.2f" % pwrFit[0])
	ax1.legend()
	ax1.set_xlabel("log(Frequency)")
	ax1.set_ylabel("log(LS Power)")

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
	ax2.scatter(logFreqBin, logPwrBin, color='b')
	ax2.plot(logFreqBin, np.polyval(pwrFitBin, logFreqBin), c='r', label=r"$\beta$" +": %5.2f" % pwrFitBin[0])
	ax2.legend()
	ax2.set_xlabel("log(Frequency)")
	ax2.set_ylabel("log(LS Power)")

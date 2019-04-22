#
#v404_analysis.py
#
#2019-04-02 AMB
#
#Script for analyzing sigma-z values from written data from lcPlot.
#


# Import everything needed and set plot settings

from astropy.table import Table
import numpy as np
import matplotlib.pylab as plt
plt.ion()
plt.style.use('ggplot')

def go(binned=True, showOutburst=True, degPoly=0, offset=False, figPrep=False, \
       underlayNight=True, useMedian=True, showCannedMean=False, showScaledZ04=False, \
       annoOutburst=True, forProposal=False):

	"""Reads sigma-z data and produces plots similar to Figure 4 of Zurita et al (2004). Some control variables:

        underlayNight = underlay the night-by-night statistics beneath the run-by-run averages?
        useMedian = use the median to calculate averages; otherwise use the mean
        showCannedMean = show the Z04 mean value
        showScaledZ04 = rescale the Z04 2001-2003 data by the median scale from our determination to Z04 1992-1999
        annoOutburst: annotate the outburst on the figure (rather than adding the label to the legend)

        forProposal - use style sheet defaults appropriate for printing on paper"""

        # 2019-04-18 WIC - style sheet to use
        if not forProposal:
                #plt.style.use('classic')
                plt.style.use('ggplot')
        else:
                plt.style.use('seaborn-poster')
                plt.style.use('seaborn-white')

        # 2019-04-18 WIC 
        # set the averaging method
        methAverage = np.mean
        if useMedian:
                methAverage=np.median
        
	# Load the data in format [JD, sigma-z]

	if degPoly == 2:
		if binned:
			try:
				t92 = Table.read('/Users/amblevin/Desktop/sigma92Binned.txt', format='ascii')
				t98 = Table.read('/Users/amblevin/Desktop/sigma98.txt', format='ascii') # Not binned
				t99 = Table.read('/Users/amblevin/Desktop/sigma99Binned.txt', format='ascii')
				t17 = Table.read('/Users/amblevin/Desktop/sigma17Binned.txt', format='ascii')
				t18A = Table.read('/Users/amblevin/Desktop/sigma18ABinned.txt', format='ascii')
				t18B = Table.read('/Users/amblevin/Desktop/sigma18BBinned.txt', format='ascii')
			except UnboundLocalError:
				print "One or more files not found. Check the directories indicated in the source code."
		else:
			try:
				t92 = Table.read('/Users/amblevin/Desktop/sigma92.txt', format='ascii') # SUSPECT- binTime seems wrong
				t98 = Table.read('/Users/amblevin/Desktop/sigma98.txt', format='ascii')
				t99 = Table.read('/Users/amblevin/Desktop/sigma99.txt', format='ascii')
				t17 = Table.read('/Users/amblevin/Desktop/sigma17.txt', format='ascii')
				t18A = Table.read('/Users/amblevin/Desktop/sigma18A.txt', format='ascii')
				t18B = Table.read('/Users/amblevin/Desktop/sigma18B.txt', format='ascii')
			except UnboundLocalError:
				print "One or more files not found. Check the directories indicated in the source code."
	elif degPoly == 1:
		if binned:
			try:
				t92 = Table.read('/Users/amblevin/Desktop/sigma92Binnedd1.txt', format='ascii')
				t98 = Table.read('/Users/amblevin/Desktop/sigma98d1.txt', format='ascii') # Not binned
				t99 = Table.read('/Users/amblevin/Desktop/sigma99Binnedd1.txt', format='ascii')
				t17 = Table.read('/Users/amblevin/Desktop/sigma17Binnedd1.txt', format='ascii')
				t18A = Table.read('/Users/amblevin/Desktop/sigma18ABinnedd1.txt', format='ascii')
				t18B = Table.read('/Users/amblevin/Desktop/sigma18BBinnedd1.txt', format='ascii')
			except UnboundLocalError:
				print "One or more files not found. Check the directories indicated in the source code."
		else:
			try:
				t92 = Table.read('/Users/amblevin/Desktop/sigma92d1.txt', format='ascii') # SUSPECT- binTime seems wrong
				t98 = Table.read('/Users/amblevin/Desktop/sigma98d1.txt', format='ascii')
				t99 = Table.read('/Users/amblevin/Desktop/sigma99d1.txt', format='ascii')
				t17 = Table.read('/Users/amblevin/Desktop/sigma17d1.txt', format='ascii')
				t18A = Table.read('/Users/amblevin/Desktop/sigma18Ad1.txt', format='ascii')
				t18B = Table.read('/Users/amblevin/Desktop/sigma18Bd1.txt', format='ascii')
			except UnboundLocalError:
				print "One or more files not found. Check the directories indicated in the source code."
	elif degPoly == 0:
		if binned:
			try:
				t92 = Table.read('/Users/amblevin/Desktop/sigma92Binnedd0.txt', format='ascii')
				t98 = Table.read('/Users/amblevin/Desktop/sigma98d0.txt', format='ascii') # Not binned
				t99 = Table.read('/Users/amblevin/Desktop/sigma99Binnedd0.txt', format='ascii')
				t17 = Table.read('/Users/amblevin/Desktop/sigma17Binnedd0.txt', format='ascii')
				t18A = Table.read('/Users/amblevin/Desktop/sigma18ABinnedd0.txt', format='ascii')
				t18B = Table.read('/Users/amblevin/Desktop/sigma18BBinnedd0.txt', format='ascii')
			except UnboundLocalError:
				print "One or more files not found. Check the directories indicated in the source code."
		else:
			try:
				t92 = Table.read('/Users/amblevin/Desktop/sigma92d0.txt', format='ascii') # SUSPECT- binTime seems wrong
				t98 = Table.read('/Users/amblevin/Desktop/sigma98d0.txt', format='ascii')
				t99 = Table.read('/Users/amblevin/Desktop/sigma99d0.txt', format='ascii')
				t17 = Table.read('/Users/amblevin/Desktop/sigma17d0.txt', format='ascii')
				t18A = Table.read('/Users/amblevin/Desktop/sigma18Ad0.txt', format='ascii')
				t18B = Table.read('/Users/amblevin/Desktop/sigma18Bd0.txt', format='ascii')
			except UnboundLocalError:
				print "One or more files not found. Check the directories indicated in the source code."

	# Extract all columns into arrays

	jd92 = np.asarray(t92['JD'])
	sig92 = np.asarray(t92['std'])

	jd98 = np.asarray(t98['JD'])
	sig98 = np.asarray(t98['std'])

	jd99 = np.asarray(t99['JD'])
	sig99 = np.asarray(t99['std'])

	jd17 = np.asarray(t17['JD'])
	sig17 = np.asarray(t17['std'])

	jd18A = np.asarray(t18A['JD'])
	sig18A = np.asarray(t18A['std'])

	jd18B = np.asarray(t18B['JD'])
	sig18B = np.asarray(t18B['std'])

	# Find the mean sigma-z for Z04 data

	sigZ04 = np.hstack([sig92, sig98, sig99])
	mszZ04 = methAverage(sigZ04)

	# Find the mean sigma-z for MDM data

	sigMDM = np.hstack([sig17, sig18A, sig18B])
	mszMDM = methAverage(sigMDM)

	# Find the mean sigma-z for the entire dataset

	sigAll = np.hstack([sig92, sig98, sig99, sig17, sig18A, sig18B])
	mszAll = np.mean(sigAll)

	# Get the mean values ready to be plotted

	jdAll = np.hstack([jd92, jd98, jd99, jd17, jd18A, jd18B])

	meanZ04 = np.repeat(mszZ04, len(jdAll))
	meanMDM = np.repeat(mszMDM, len(jdAll))
	meanAll = np.repeat(mszAll, len(jdAll))

	mZ04 = np.mean([0.040, 0.034, 0.025, 0.042, 0.038, 0.041]) # From Table 3 of Z04
	aZ04 = np.repeat(mZ04, len(jdAll))

	jdAllBut18B = np.hstack([jd92, jd98, jd99, jd17, jd18A])
	sigAllBut18B = np.hstack([sig92, sig98, sig99, sig17, sig18A])

	# Plot sigma-z vs. time for each night

        # 2019-04-13 WIC - hack to provide an identifier for the dataset
        bZ04all = jdAllBut18B < 57000
        
	f1 = plt.figure(1)
	f1.clf()
	ax1 = f1.add_subplot(111)
	ax1.scatter(jdAllBut18B[bZ04all], sigAllBut18B[bZ04all], color='k', alpha=0.5, marker='v', label='Z04', s=25, zorder=3)
        ax1.scatter(jdAllBut18B[~bZ04all], sigAllBut18B[~bZ04all], color='b', alpha=0.5, marker='o', label='MDM', s=25, zorder=3)
	ax1.scatter(jd18B, sig18B, color='r', alpha=0.6, label='MDM (sparse)', s=16, zorder=3)
	ax1.hlines(meanAll, 47000, 60000, color='0.2', lineStyles='-.', label="Mean (All data)", lw=1, zorder=0)

        # show the Z04 value?
        if showCannedMean:
                ax1.hlines(0.037, 47000, 60000, color='0.2', lineStyles='dashed', label='Mean (literature, Z04)', lw=1, zorder=1)

        if not figPrep:
		ax1.plot(jdAll, meanZ04, 'g--', label="Mean (Zurita data)")
		#ax1.plot(jdAll, meanMDM, 'b--', label="Mean (MDM data)")
		ax1.hlines(meanMDM, 47000, 60000, color='blue', lineStyles='--', label="Mean(MDM Data)")
		ax1.plot(jdAll, aZ04, 'k--', label="Mean (Z04 Figure 4)")
	if showOutburst:
		plt.vlines(47667, ymin=0., ymax=0.08, color='g', label="1989 Outburst") # Date obtained from X-Ray Binaries Book
		plt.vlines(x=57194, ymin=0., ymax=0.08, color='b', label="2015 Outburst") # Date obtained from Munoz-Darias et al (2016)
	ax1.legend(loc=9)
	plt.xlabel("JD - 2 400 000 (days)")
	plt.ylabel(r"$\sigma_{z}$ (mag)")
	ax1.set_xlim(47000, 60000)
	ax1.set_ylim(0,0.08)
	if figPrep:
		plt.title("Nightly Flare Activity")
	else:
		if binned:
			plt.title("Nightly Flare Activity (Binned); Degree %i" % degPoly)
		else:
			plt.title("Nightly Flare Activity; Degree %i" % degPoly)
	plt.show()

	# Find the average sigma-z for each year

        # WIC 2019-04-16: updated to include the Z04 table values
        z04_sz_92 = 0.040 
        z04_sz_98 = 0.034 
        z04_sz_99 = 0.025

        z04_sz_01 = 0.042
        z04_sz_02 = 0.038
        z04_sz_03 = 0.041

	mjd92 = methAverage(jd92)
	msz92 = methAverage(sig92)
	print "msz92: %.5f, z04_sz_92: %.5f: msz/z04 = %.5f" % ( msz92, z04_sz_92, msz92/z04_sz_92)

	mjd98 = methAverage(jd98)
	msz98 = methAverage(sig98)
	#print "msz98: ", msz98
	print "msz98: %.5f, z04_sz_98: %.5f: msz/z04 = %.5f" % ( msz98, z04_sz_98, msz98/z04_sz_98)

	mjd99 = methAverage(jd99)
	msz99 = methAverage(sig99)
	#print "msz99: ", msz99
	print "msz92: %.5f, z04_sz_99: %.5f: msz/z04 = %.5f" % ( msz99, z04_sz_99, msz99/z04_sz_99)

        # refactored into an array for calculation:
        v_msz = np.array([msz92, msz98, msz99])
        v_z04 = np.array([z04_sz_92, z04_sz_98, z04_sz_99])
        
        print "Mean (Stddev) of msz/z04: %.5f (%.5f)" % \
                (np.mean(v_msz/v_z04), np.std(v_msz/v_z04))

        # stack into array (use the rough midpoint in MJD for the times)
        sf_z04_to_own = np.mean(v_msz/v_z04)
        v_z04late = np.array([z04_sz_01, z04_sz_02, z04_sz_03])
        mjd_z04late = np.array([52131., 52426., 52839.])
        
	mjd17 = methAverage(jd17)
	msz17 = methAverage(sig17)

	mjd18A = methAverage(jd18A)
	msz18A = methAverage(sig18A)
        
	mjd18B = np.mean(jd18B)
	msz18B = np.mean(sig18B)

	mjdAll = np.hstack([mjd92, mjd98, mjd99, mjd17, mjd18A, mjd18B])
	msAll = np.hstack([msz92, msz98, msz99, msz17, msz18A, msz18B])

	mjdAllBut18B = np.hstack([mjd92, mjd98, mjd99, mjd17, mjd18A])
	msAllBut18B = np.hstack([msz92, msz98, msz99, msz17, msz18A])

        # 2019-04-13 WIC - hack to separate out the plot symbols for Z04 and for our data
        bZ04mean = mjdAllBut18B < 57000

	# Find dashed lines

	m92and98 = np.mean([msz92, msz98, msz99])
	mMDM = np.mean([msz17, msz18A, msz18B])
	mAll = np.mean([msz92, msz98, msz99, msz17, msz18A, msz18B])

	dl92and98 = np.repeat(m92and98, len(mjdAll))
	dlMDM = np.repeat(mMDM, len(mjdAll))
	dlAll = np.repeat(mAll, len(mjdAll))
	dlZ04 = np.repeat(mZ04, len(mjdAll))
        
	# Plot sigma-z vs. time for each year

	f2 = plt.figure(2)
	f2.clf()
	ax2 = f2.add_subplot(111)

        # 2019-04-11: include the individual nightly points underneath the averages?
        if underlayNight:
                ax2.scatter(jdAllBut18B[bZ04all], sigAllBut18B[bZ04all], color='k', alpha=0.25, marker='v', s=9, zorder=2)
                ax2.scatter(jdAllBut18B[~bZ04all], sigAllBut18B[~bZ04all], color='b', alpha=0.25, marker='o', s=9)
                ax2.scatter(jd18B, sig18B, color='r', alpha=0.25, s=9)
        
	ax2.scatter(mjdAllBut18B[bZ04mean], msAllBut18B[bZ04mean], color='k', alpha=1.0, marker='v', label='Z04 (recalculated)', s=49, zorder=3)
        # show the re-scaled late Z04 epochs (2001-2003)
        if showScaledZ04:
                dum4 = ax2.scatter(mjd_z04late, v_z04late * sf_z04_to_own, color='k', alpha=0.8, marker='^', s=40, zorder=3, \
                                   label='Z04 (literature) *%.2f' % (sf_z04_to_own))
                
        ax2.scatter(mjdAllBut18B[~bZ04mean], msAllBut18B[~bZ04mean], color='b', alpha=1.0, marker='o', label='MDM', s=49, zorder=4)
	ax2.scatter(mjd18B, msz18B, color='r', label='MDM (sparse)', s=36, zorder=4)
	#ax2.hlines(dlAll, 47000, 60000, color='0.2', lineStyles='--', label="Mean (All data)", lw=1)

        # 2019-04-18 WIC - made this optional (since if we have the rescaled 200? as well, it's no longer all the data
        if not showScaledZ04:
                ax2.hlines(dlAll, 47000, 60000, color='0.2', lineStyles='-.', label="Mean (All data)", lw=1, zorder=0)

        # show the Z04 value?
        if showCannedMean:
                meanCanned = 0.037
                labelCanned = 'Literature mean (Z04)'
                if showScaledZ04:
                        meanCanned *= sf_z04_to_own
                        labelCanned = '%s * %.2f' % (labelCanned, sf_z04_to_own)
                ax2.hlines(meanCanned, 47000, 60000, color='0.2', lineStyles='dashed', label=labelCanned, lw=1, zorder=1)
        
	if not figPrep:
		ax2.plot(mjdAll, dlZ04, 'k--', label="Mean (Z04 Figure 4)")
		ax2.plot(mjdAll, dl92and98, 'g--', label="Mean (Zurita data)")
		ax2.plot(mjdAll, dlMDM, 'b--', label="Mean (MDM Data)")
		ax2.plot(mjdAll, dlAll, 'r--', label="Mean (All data)")
	if showOutburst:

                label89 = '1989 Outburst'
                label15 = '2015 Outburst'

                # 2019-04-16 WIC - I prefer to annotate the figure rather than add to the legend:
                if annoOutburst:
                        fsz=12
                        anno1 = ax2.annotate('1989 Outburst', (47667.+100, 0.001), xycoords='data', ha='left', va='bottom', color='g', \
                                             fontsize=fsz)
                        anno2 = ax2.annotate('2015 Outburst', (57194.-100, 0.001), xycoords='data', ha='right', va='bottom', color='b', \
                                             fontsize=fsz)
                        label89=''
                        label15=''
                
		plt.vlines(47667, ymin=0., ymax=0.08, color='g', label=label89) # Date obtained from X-Ray Binaries Book
		plt.vlines(x=57194, ymin=0., ymax=0.08, color='b', label=label15) # Date obtained from Munoz-Darias et al (2016)

        leg = ax2.legend(loc=9, frameon=True)
        leg.get_frame().set_facecolor('w')
        
	ax2.set_xlim(47000, 60000)
	ax2.set_ylim(0, 0.08)
	plt.xlabel("JD - 2 400 000 (days)")
	plt.ylabel(r"$\sigma_{z}$  (mag)")
	if figPrep:

                # adjust string to use for title
                sMean = 'Mean'
                if useMedian:
                        sMean = 'Median'
                
		plt.title("%s flaring within each observing run" % (sMean))
	else:
		if binned:
			plt.title ("Yearly Flare Activity (Binned); Degree %i" % degPoly)
		else:
			plt.title("Yearly Flare Activity; Degree %i" % degPoly)
	plt.show()

        # save the figures
        f1.savefig('sigma_perNight.png')
        f2.savefig('sigma_perRun.png')

#
#plotDIA.py
#
#2018-11-09 AMB
#

import readphotometry
from astropy.table import Table, Column
import cPickle as pickle
import numpy as np
import matplotlib.pylab as plt
plt.ion()
plt.style.use('ggplot')

def load17table():
	DDIA17 = pickle.load(open('/Users/Shared/Data/MDM/2017B/proc/DIA/alignedCropped/Output_dir/lightcurves.cPickle', 'r'))
	
	tbl17 = Table()
	tbl17['time'] = np.asarray(DDIA17['times'])
	tbl17['flux'] = np.asarray(DDIA17['fluxAll'])
	tbl17['fluxErr'] = np.asarray(DDIA17['unctAll'])

	flux17 = tbl17['flux']
	fluxErr17 = tbl17['fluxErr']

	tbl17B = Table()

	tbl17B['time'] = np.asarray(DDIA17['times'])

	v17 = flux17[:,15]
	col_v17 = Column(name='fluxV404', data=v17)
	tbl17B.add_column(col_v17)

	c17 = flux17[:,21]
	col_c17 = Column(name='fluxContam', data=c17)
	tbl17B.add_column(col_c17)

	c417 = flux17[:,3]
	col_c417 = Column(name='fluxC4', data=c417)
	tbl17B.add_column(col_c417)

	n17 = flux17[:,17]
	col_n17 = Column(name='fluxNeighbor', data=n17)
	tbl17B.add_column(col_n17)


	v17e = fluxErr17[:,15]
	col_v17e = Column(name='fluxErrV404', data=v17e)
	tbl17B.add_column(col_v17e)

	c17e = fluxErr17[:,21]
	col_c17e = Column(name='fluxErrContam', data=c17e)
	tbl17B.add_column(col_c17e)

	c417e = fluxErr17[:,3]
	col_c417e = Column(name='fluxErrC4', data=c417e)
	tbl17B.add_column(col_c417e)

	n17e = fluxErr17[:,17]
	col_n17e = Column(name='fluxErrNeighbor', data=n17e)
	tbl17B.add_column(col_n17e)

	return tbl17B


#test = load17table()

#print test


def load18table():
	DDIA18 = pickle.load(open('/Users/Shared/Data/MDM/2018A/proc/DIA/allFiveNights/alignedCropped/Output_dir_backup/lightcurves.cPickle', 'r'))

	tbl18 = Table()

	tbl18['time'] = np.asarray(DDIA18['times'])
	tbl18['flux'] = np.asarray(DDIA18['fluxAll'])
	tbl18['fluxErr'] = np.asarray(DDIA18['unctAll'])

	flux18 = tbl18['flux']
	fluxErr18 = tbl18['fluxErr']

	tbl18A = Table()

	tbl18A['time'] = np.asarray(DDIA18['times'])

	v18 = flux18[:,4]
	col_v18 = Column(name='fluxV404', data=v18)
	tbl18A.add_column(col_v18)

	c18 = flux18[:,8]
	col_c18 = Column(name='fluxContam', data=c18)
	tbl18A.add_column(col_c18)

	c418 = flux18[:,0]
	col_c418 = Column(name='fluxC4', data=c418)
	tbl18A.add_column(col_c418)

	n18 = flux18[:,6]
	col_n18 = Column(name='fluxNeighbor', data=n18)
	tbl18A.add_column(col_n18)

	r18 = flux18[:,5]
	col_r18 = Column(name='fluxRandomStar', data=r18)
	tbl18A.add_column(col_r18)


	v18e = fluxErr18[:,4]
	col_v18e = Column(name='fluxErrV404', data=v18e)
	tbl18A.add_column(col_v18e)

	c18e = fluxErr18[:,8]
	col_c18e = Column(name='fluxErrContam', data=c18e)
	tbl18A.add_column(col_c18e)

	c418e = fluxErr18[:,0]
	col_c418e = Column(name='fluxErrC4', data=c418e)
	tbl18A.add_column(col_c418e)

	n18e = fluxErr18[:,6]
	col_n18e = Column(name='fluxErrNeighbor', data=n18e)
	tbl18A.add_column(col_n18e)

	r18e = fluxErr18[:,5]
	col_r18e = Column(name='fluxErrRandomStar', data=r18e)
	tbl18A.add_column(col_r18e)

	return tbl18A

#test2 = load18table()

#print test2


# Note to self: If/when DIA is done on 2018B photometry, call the first table 'tbl2018' and the second 'tbl18B'


def plot(data='2018', errorbars=True, nonEssentials=False, \
         c404=0., cneib=0., ccontam=0.):
	if data == '2017':
		tbl = load17table()
	elif data == '2018':
		tbl = load18table()
	else:
		print "The value for 'data' must be either '2017' or '2018'."

	#print tbl


	# Defining Variables

	time = tbl['time']
	fluxV = tbl['fluxV404']
	fluxC = tbl['fluxContam']
	fluxC4 = tbl['fluxC4']
	fluxN = tbl['fluxNeighbor']
	dyV = tbl['fluxErrV404']
	dyC = tbl['fluxErrContam']
	dyC4 = tbl['fluxErrC4']
	dyN = tbl['fluxErrNeighbor']
	fluxR = tbl['fluxRandomStar']
	dyR = tbl['fluxErrRandomStar']


	# Plot of time vs. flux(V404)

	fig1 = plt.figure(1)
	fig1.clf()
	plt.scatter(time, fluxV, alpha=0.4, s=16, zorder=25)
	if errorbars:
		plt.errorbar(time, fluxV, yerr=dyV, fmt='o', ms=1, ecolor='0.3', alpha=0.2, zorder=10)
	plt.xlabel('Time')
	plt.ylabel('Flux (V404)')

	# Plot of time vs. flux(Contaminant)

	fig2 = plt.figure(2)
	fig2.clf()
	plt.scatter(time, fluxC, alpha=0.4, s=16, zorder=25)
	if errorbars:
		plt.errorbar(time, fluxC, yerr=dyC, fmt='o', ms=1, ecolor='0.3', alpha=0.2, zorder=10)
	plt.xlabel('Time')
	plt.ylabel('Flux (Contaminant)')

	# Plot of time vs. flux(C4)

	fig3 = plt.figure(3)
	fig3.clf()
	plt.scatter(time, fluxC4, alpha=0.4, s=16, zorder=25)
	if errorbars:
		plt.errorbar(time, fluxC4, yerr=dyC4, fmt='o', ms=1, ecolor='0.3', alpha=0.2, zorder=10)
	plt.xlabel('Time')
	plt.ylabel('Flux (C4)')

	# Plot of time vs. flux(Neighbor)

	fig4 = plt.figure(4)
	fig4.clf()
	plt.scatter(time, fluxN, alpha=0.4, s=16, zorder=25)
	if errorbars:
		plt.errorbar(time, fluxN, yerr=dyN, fmt='o', ms=1, ecolor='0.3', alpha=0.2, zorder=10)
	plt.xlabel('Time')
	plt.ylabel('Flux (Neighbor)')


	# Plot of Flux(V404) vs. Flux(Contaminant)

	fig5 = plt.figure(5)
	fig5.clf()
	dum5 = plt.scatter(fluxV+c404, fluxC+ccontam, alpha=0.4, s=16, zorder=25, c=time, cmap='hsv')
	if errorbars:
		plt.errorbar(fluxV+c404, fluxC+ccontam, xerr=dyV, yerr=dyC, fmt='o', ms=1, ecolor='0.3', alpha=0.2, zorder=10)
	plt.xlabel('Flux (V404)')
	plt.ylabel('Flux (Contaminant)')
	plt.colorbar(dum5)

        # let's do a ratio plot
        if ccontam > 0. and c404 > 0.:
                absV = fluxV + c404
                absC = fluxC + ccontam
                
                ratioContam = absC / (absV + absC)
                unctyRatioSquared = (dyV**2 + dyC**2)/(absV + absC)**2 + \
                                    (dyC/absC)**2
                unctyOfRatio = ratioContam * np.sqrt(unctyRatioSquared)

                # now let's plot this
                fig55=plt.figure(55)
                fig55.clf()
                dum55 = plt.scatter(absV, ratioContam, alpha=0.4, s=16, \
                                    zorder=25, c=time, cmap='hsv')

                if errorbars:
                        dum55b = plt.errorbar(absV, ratioContam, xerr=dyV, \
                                              yerr=unctyOfRatio, fmt='o', ms=1, \
                                              ecolor='0.3', alpha=0.2, zorder=10)
                plt.xlabel('Flux (V404)')
                plt.ylabel('(Contaminant / (V404 + contaminant))')
                plt.colorbar(dum55)
                
	# Plot of Flux(V404) vs. Flux(Neighbor)

	fig6 = plt.figure(6)
	fig6.clf()
	dum6 = plt.scatter(fluxV, fluxN, alpha=0.4, s=16, zorder=25, c=time, cmap='hsv')
	if errorbars:
		plt.errorbar(fluxV, fluxN, xerr=dyV, yerr=dyN, fmt='o', ms=1, ecolor='0.3', alpha=0.2, zorder=10)
	plt.xlabel('Flux (V404)')
	plt.ylabel('Flux (Neighbor)')
	plt.colorbar(dum6)

	# Plot of Flux(V404) vs. Flux(C4)

	fig7 = plt.figure(7)
	fig7.clf()
	dum7 = plt.scatter(fluxV, fluxC4, alpha=0.4, s=16, zorder=25, c=time, cmap='hsv')
	if errorbars:
		plt.errorbar(fluxV, fluxC4, xerr=dyV, yerr=dyC4, fmt='o', ms=1, ecolor='0.3', alpha=0.2, zorder=10)
	plt.xlabel('Flux (V404)')
	plt.ylabel('Flux (C4)')
	plt.colorbar(dum7)

	# Plot of Flux(Contaminant) vs. Flux(C4)

	if nonEssentials:
		fig8 = plt.figure(8)
		fig8.clf()
		dum8 = plt.scatter(fluxC, fluxC4, alpha=0.4, s=16, zorder=25, c=time, cmap='hsv')
		if errorbars:
			plt.errorbar(fluxC, fluxC4, xerr=dyC, yerr=dyC4, fmt='o', ms=1, ecolor='0.3', alpha=0.2, zorder=10)
		plt.xlabel('Flux (Contaminant)')
		plt.ylabel('Flux (C4)')
		plt.colorbar(dum8)

		# Plot of Flux(Contaminant) vs. Flux(Neighbor)

		fig9 = plt.figure(9)
		fig9.clf()
		dum9 = plt.scatter(fluxC, fluxN, alpha=0.4, s=16, zorder=25, c=time, cmap='hsv')
		if errorbars:
			plt.errorbar(fluxC, fluxN, xerr=dyC, yerr=dyN, fmt='o', ms=1, ecolor='0.3', alpha=0.2, zorder=10)
		plt.xlabel('Flux (Contaminant)')
		plt.ylabel('Flux (Neighbor)')
		plt.colorbar(dum9)

		# Plot of Flux(Neighbor) vs. Flux(C4)

		fig10 = plt.figure(10)
		fig10.clf()
		dum10 = plt.scatter(fluxN, fluxC4, alpha=0.4, s=16, zorder=25, c=time, cmap='hsv')
		if errorbars:
			plt.errorbar(fluxN, fluxC4, xerr=dyN, yerr=dyC4, fmt='o', ms=1, ecolor='0.3', alpha=0.2, zorder=10)
		plt.xlabel('Flux (Neighbor)')
		plt.ylabel('Flux (C4)')
		plt.colorbar(dum10)

		# Plot of Flux(Random Star) vs. Flux(C4)

		fig11 = plt.figure(11)
		fig11.clf()
		dum11 = plt.scatter(fluxR, fluxC4, alpha=0.4, s=16, zorder=25, c=time, cmap='hsv')
		if errorbars:
			plt.errorbar(fluxR, fluxC4, xerr=dyR, yerr=dyC4, fmt='o', ms=1, ecolor='0.3', alpha=0.2, zorder=10)
		plt.xlabel('Flux (Random Star (#6 on Map))')
		plt.ylabel('Flux (C4)')
		plt.colorbar(dum11)


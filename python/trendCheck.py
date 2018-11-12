#
#trendCheck.py
#
#2018-11-12 AMB
#
#Quick script: No arguments or definitions; plots on import/reload!
#

import matplotlib.pylab as plt
plt.ion()
plt.style.use('ggplot')
from astropy.table import Table, Column
import numpy as np

aDum = np.genfromtxt("/Users/amblevin/Desktop/RvC4.txt", unpack=False)

tPho = Table()

nCols = np.shape(aDum)[-1]
time = aDum[:,0]
fR = aDum[:,1]
feR = aDum[:,2]
fC = aDum[:,3]
feC = aDum[:,4]

tPho['JD - 2 400 000'] = Column(time)
tPho['Relative Flux(Star #6/C5)'] = Column(fR)
tPho['Relative Flux Error(Star #6/C5)'] = Column(feR)
tPho['Relative Flux(C4/C5)'] = Column(fC)
tPho['Relative Flux Error(C4/C5)'] = Column(feC)

#print tPho

# Defining Variables

jd = tPho['JD - 2 400 000']
fluxR = tPho['Relative Flux(Star #6/C5)']
dyR = tPho['Relative Flux Error(Star #6/C5)']
fluxC = tPho['Relative Flux(C4/C5)']
dyC = tPho['Relative Flux Error(C4/C5)']

# Plot of flux(R/C5) vs. Flux(C4/C5)

fig1 = plt.figure(1)
fig1.clf()
dum = plt.scatter(fluxR, fluxC, alpha=0.4, s=16, zorder=25, c=jd, cmap='hsv')
plt.errorbar(fluxR, fluxC, xerr=dyR, yerr=dyC, fmt='o', ms=1, ecolor='0.3', alpha=0.2, zorder=10)
plt.xlabel('Relative Flux (Star #6 / C5 )')
plt.ylabel('Relative Flux (C4 / C5)')
plt.colorbar(dum)
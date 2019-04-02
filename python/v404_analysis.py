#
#v404_analysis.py
#
#2019-04-02 AMB
#
#Script for analyzing sigma-z values from written data from lcPlot.
#
#No arguments or definitions; plots on import/reload!
#

# Import everything needed and set plot settings

from astropy.table import Table
import numpy as np
import matplotlib.pylab as plt
plt.ion()
plt.style.use('ggplot')

# Load the data in format [JD, sigma-z]

try:
	t92 = Table.read('/Users/amblevin/Desktop/sigma92.txt', format='ascii') # Still treating this file as suspect
	t98 = Table.read('/Users/amblevin/Desktop/sigma98.txt', format='ascii')
	t17 = Table.read('/Users/amblevin/Desktop/sigma17.txt', format='ascii')
	t18A = Table.read('/Users/amblevin/Desktop/sigma18A.txt', format='ascii')
	t18B = Table.read('/Users/amblevin/Desktop/sigma18B.txt', format='ascii')
except UnboundLocalError:
	print "One or more files not found. Check the directories indicated in the source code."

# Extract all columns into arrays

jd92 = np.asarray(t92['JD'])
sig92 = np.asarray(t92['std'])

jd98 = np.asarray(t98['JD'])
sig98 = np.asarray(t98['std'])

jd17 = np.asarray(t17['JD'])
sig17 = np.asarray(t17['std'])

jd18A = np.asarray(t18A['JD'])
sig18A = np.asarray(t18A['std'])

jd18B = np.asarray(t18B['JD'])
sig18B = np.asarray(t18B['std'])

# Find the mean sigma-z for Z04 data

sigZ04 = np.hstack([sig92, sig98])
mszZ04 = np.mean(sigZ04)

# Find the mean sigma-z for MDM data

sigMDM = np.hstack([sig17, sig18A, sig18B])
mszMDM = np.mean(sigMDM)

# Find the mean sigma-z for the entire dataset

sigAll = np.hstack([sig92, sig98, sig17, sig18A, sig18B])
mszAll = np.mean(sigAll)

# Get the mean values ready to be plotted

jdAll = np.hstack([jd92, jd98, jd17, jd18A, jd18B])

meanZ04 = np.repeat(mszZ04, len(jdAll))
meanMDM = np.repeat(mszMDM, len(jdAll))
meanAll = np.repeat(mszAll, len(jdAll))

mZ04 = np.mean([0.040, 0.034, 0.025, 0.042, 0.038, 0.041]) # From Table 3 of Z04
aZ04 = np.repeat(mZ04, len(jdAll))

# Plot sigma-z vs. time for each night

f1 = plt.figure(1)
f1.clf()
ax1 = f1.add_subplot(111)
ax1.scatter(jdAll, sigAll, color='k')
ax1.plot(jdAll, meanAll, 'r--', label="Mean (All data)")
ax1.plot(jdAll, meanZ04, 'g--', label="Mean (Zurita data)")
ax1.plot(jdAll, meanMDM, 'b--', label="Mean (MDM data)")
ax1.plot(jdAll, aZ04, 'k--', label="Mean (Z04 Figure 4)")
ax1.legend()
plt.xlabel("JD - 2 400 000")
plt.ylabel(r"$\sigma_{z}$")
plt.title("Nightly Flare Activity")
plt.show()

# Find the average sigma-z for each year

mjd92 = np.mean(jd92)
msz92 = np.mean(sig92)

mjd98 = np.mean(jd98)
msz98 = np.mean(sig98)

mjd17 = np.mean(jd17)
msz17 = np.mean(sig17)

mjd18A = np.mean(jd18A)
msz18A = np.mean(sig18A)

mjd18B = np.mean(jd18B)
msz18B = np.mean(sig18B)

mjdAll = np.hstack([mjd92, mjd98, mjd17, mjd18A, mjd18B])
msAll = np.hstack([msz92, msz98, msz17, msz18A, msz18B])

# Find dashed lines

m92and98 = np.mean([msz92, msz98])
mMDM = np.mean([msz17, msz18A, msz18B])
mAll = np.mean([msz92, msz98, msz17, msz18A, msz18B])

dl92and98 = np.repeat(m92and98, len(mjdAll))
dlMDM = np.repeat(mMDM, len(mjdAll))
dlAll = np.repeat(mAll, len(mjdAll))
dlZ04 = np.repeat(mZ04, len(mjdAll))

# Plot sigma-z vs. time for each year

f2 = plt.figure(2)
f2.clf()
ax2 = f2.add_subplot(111)
ax2.scatter(mjdAll, msAll, color='k')
ax2.plot(mjdAll, dlZ04, 'k--', label="Mean (Z04 Figure 4)")
ax2.plot(mjdAll, dl92and98, 'g--', label="Mean (Zurita data)")
ax2.plot(mjdAll, dlMDM, 'b--', label="Mean (MDM Data)")
ax2.plot(mjdAll, dlAll, 'r--', label="Mean (All data)")
ax2.legend()
ax2.set_ylim(0, 0.07)
plt.xlabel("JD - 2 400 000")
plt.ylabel(r"$\sigma_{z}$")
plt.title("Yearly Flare Activity")
plt.show()
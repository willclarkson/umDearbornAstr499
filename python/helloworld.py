#cleaning up angle data


from astropy.table import Table
import num.py as np
import matplotlib.pylab as plt
import seaborn as sns
from astropy.stats import sigma_clip
tDum= Table.read('Gaiacatalog0.ASC', format='ascii.sextractor')
xCut=5000.
bGood=tDum['FLUX_ISO']
dum= ax.scatter(tDum['FLUX_ISO'][bGood],tDum['THETA_IMAGE'][bGood],color='c', alpha=0.7,marker='v')
filtered_data=sigma_clip(tDum['THETA_IMAGE'][bGood],sigma=4,iters=5)


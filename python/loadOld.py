#
# loadOld.py
#

# 2018-05-07

from astropy.table import Table, Column
import numpy as np

def loadPhot(fileIn='199208elip_R.dat'):

    """Quick utility to read in the zurita photometry, using the ad hoc
rules for reading in the casares collection of zurita04 data

    """

    # what might the comments strings be...
    lCommens = ['%', '!', '#']

    aData = np.array([])
    for sCommen in lCommens:
        try:
            aDum = np.genfromtxt(fileIn, comments=sCommen, unpack=False)
        except:
            badRead = True

    # Now we have to decide what this all means...
    tPho = Table()

    # are we dealing with times or phases here?
    nCols = np.shape(aDum)[-1]
    timeOrPhase = aDum[:,0]
    magn = aDum[:,1]

    sCol = 'time'
    if np.max(timeOrPhase) < 1.0:
        sCol = 'phase'

    tPho[sCol] = Column(timeOrPhase)

    # Flux or magnitude?
    if nCols > 2: # For some reason, 2003 data was displaying incorrectly when using the values themselves..
        tPho['mag'] = Column(magn)
    else:
        tPho['flux'] = Column(magn)
        
    # now for the other columns
    if nCols > 2:
        tPho['magErr'] = Column(aDum[:,2])

    if nCols > 3:
        tPho['magC'] = Column(aDum[:,3])
        tPho['magCerr'] = Column(aDum[:,4])

    # attach metadata
    tPho.meta['file'] = fileIn[:]
        
    return tPho
        

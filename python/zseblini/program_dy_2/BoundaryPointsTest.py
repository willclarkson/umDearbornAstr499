#Boundary points testing
from astropy.table import Table, Column
import matplotlib.pylab as plt
from astropy.io import ascii
from astropy.io import fits
#import seaborn as sns
from astropy.stats import sigma_clip
import numpy as np
import astropy.units as u
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from matplotlib.pylab import quiver
from numpy import multiply

# For world coordinate system
from astropy.wcs import WCS

i=0
plt.ion()

 
def go(fCat='GaiaCatalog0.ASC', \
       fHeader='V404_Cyg_adOFF-012_R_120sTest_MAPPED.fit', \
       colBlobLength='A_IMAGE', blobLenDefault=5.):
    
    """Plotting the Co-ordinates in a Click-Animated Sequence"""
    
   
    #Objects with Flux > 5000
    # tDUM = Table.read('GaiaCatalog0.ASC', format='ascii.sextractor')
    tDUM = Table.read(fCat, format='ascii.sextractor')
   
    
    # Let's set a conditional on whether our blob length column is present
    if colBlobLength in tDUM.colnames:
        vBlobLength = tDUM[colBlobLength]
    else:
        print("programday_1 INFO - column %s not found. Generating a dummy length" %(colBlobLength))
    # let's add our dummy length column to the table.
        colBlobLength = "%s_GEN" % (colBlobLength)
        tDUM[colBlobLength] = np.repeat(blobLenDefault, len(tDUM))

    # return

    #finding time from the fits header
    #hdul = fits.open('V404_Cyg_adOFF-012_R_120sTest_MAPPED.fit')
    ##hdul.info()
    #myHeader = fits.getheader('V404_Cyg_adOFF-012_R_120sTest_MAPPED.fit')
    myHeader = fits.getheader(fHeader)
    time=myHeader['DATE-OBS']
    #print(time)

    # Let's promote the world coordinate system parsing out to here,
    # since we're going to need it whatever we do
    wcs = WCS(myHeader)


    #Define the good data
    xCut = 5000
    bGood = tDUM['FLUX_ISO'] > xCut
    # let's get the image dimensions from the header for our frame
    # boundary rectangle:
    nX = myHeader['NAXIS1']
    nY = myHeader['NAXIS2']

    boundsX = np.asarray([0., nX, nX, 0.], 'float')
    boundsY = np.asarray([0., 0., nY, nY], 'float')

    boundsX = np.hstack((boundsX, boundsX[0]))
    boundsY = np.hstack((boundsY, boundsY[0]))

    print(boundsX)
    print(boundsY)
        # Now that we've read in the header for the image, use it to
        # convert pixel coords to sky coords if they didn't come in with
        # the catalog. Reference:
    
     # https://docs.astropy.org/en/stable/wcs/
 
    def fig6(fig):
            fig.suptitle('Object Locations in Azimuth and Altitude')
            fig.clf()
            yLabel6=plt.ylabel('Altitude')
            xLabel6=plt.xlabel('Azimuth')
            ax6=fig.add_subplot(111)
            #dum6 = ax6.scatter(myAz1,myAlt1,marker='.')
            #dum66 = ax6.scatter(myAz2,myAlt2,marker='.')

            #Finding the border in ALTAZ
    
            FnY = myHeader['NAXIS1']
            FnX = myHeader['NAXIS2']
            boundsX0 = np.asarray([0.1, FnX, FnX, 0.1], 'float')
            boundsY0 = np.asarray([0.1, 0.1, FnY, FnY], 'float')
            SnX,SnY =wcs.all_pix2world(boundsX0, boundsY0, 0)
            #tDUM['TnX'] = Column(SnX, unit=u.deg)
            #tDUM['TnY'] = Column(SnY, unit=u.deg)
    
    
            #boundsX1 = np.asarray([0., SnX, SnX, 0.], 'float')
            #boundsY1 = np.asarray([0., 0., SnY, SnY], 'float')

            boundsX2 = np.hstack((SnX, SnX[0]))
            boundsY2 = np.hstack((SnY, SnY[0]))
            dum666= ax6.plot(boundsX2,boundsY2)
            print(FnX)
            print(FnY)
            print(SnX)
            print(SnY)
            print(boundsX2)
            print(boundsY2)
    #SWITCH FIGS
    switch_figs =  {
        0: fig6
    }


    def onclick1(fig):
        global i
        print(i)
        fig.clear()
        i += 1
        i %= 1
        switch_figs[i](fig)
        plt.draw()

    
    fig = plt.figure()
    switch_figs[0](fig)
    #fig.canvas.mpl_connect('button_press_event', lambda event: onclick1(fig))

    plt.show()


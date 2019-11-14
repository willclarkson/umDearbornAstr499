from astropy.table import Table
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

plt.ion()

def go(fCat='GaiaCatalog0.ASC', fHeader='V404_Cyg_adOFF-012_R_120sTest_MAPPED.fit', \
       colBlobLength='A_IMAGE', blobLenDefault=5.):
    
    """Try quiver plot"""

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

#        print(tDUM.colnames)

#        print tDUM[0:3]

# return

    fig1=plt.figure(1)
    fig1.suptitle('Objects with Flux Greater than 5000')
    fig1.clf()
    yLabel1=plt.ylabel('Theta')
    xLabel1=plt.xlabel('Flux')
    ax=fig1.add_subplot(111)
    xCut = 5000.
    bGood = tDUM['FLUX_ISO'] > xCut
    dum = ax.scatter(tDUM['FLUX_ISO'][bGood],tDUM['THETA_IMAGE'][bGood])


    #Histogram of Good Data
    fig2 = plt.figure(2)
    fig2.suptitle('Angle of Trail of Objects With Flux Greater than 5000')
    fig2.clf()
    ax2 = fig2.add_subplot(211)
    dum2 = ax2.hist(tDUM['THETA_IMAGE'][bGood], bins=150)
    label2=plt.ylabel('Selected Data')

    ax22 = fig2.add_subplot(212)
    filtered_data = sigma_clip(tDUM['THETA_IMAGE'][bGood],sigma=4, iters=5)
    dum22=ax22.hist(filtered_data, bins=150)
    xLabel2=plt.xlabel('Theta-')
    label22=plt.ylabel('Filtered Data')


    #Statistical Parameters
    print np.mean(filtered_data)
    print np.std(filtered_data)

    #finding time from the fits header
    #hdul = fits.open('V404_Cyg_adOFF-012_R_120sTest_MAPPED.fit')
    ##hdul.info()
    #myHeader = fits.getheader('V404_Cyg_adOFF-012_R_120sTest_MAPPED.fit')
    myHeader = fits.getheader(fHeader)
    time=myHeader['DATE-OBS']
    print(time)
    UMD_Observatory= EarthLocation(lat=41.32*u.deg, lon=-83.24*u.deg )
    objPos = SkyCoord(ra=tDUM['ALPHA_J2000'][bGood], dec=tDUM['DELTA_J2000'][bGood], frame='fk5')

    MappedAltAz=objPos.transform_to(AltAz(obstime=time,location=UMD_Observatory))
    myAz = np.asarray(MappedAltAz.az)
    myAlt = np.asarray(MappedAltAz.alt)
    print(myAlt)
    print(myAz)

    fig4=plt.figure(4)
    fig4.suptitle('Object Locations in Azimuth and Altitude')
    fig4.clf()
    yLabel4=plt.ylabel('Altitude')
    xLabel4=plt.xlabel('Azimuth')
    ax4=fig4.add_subplot(111)
    dum4 = ax4.scatter(myAz,myAlt,marker='.')


    # let's try a quiver plot. We will assemble the x- and y- components of the vectors
    myOrientation = tDUM['THETA_IMAGE'][bGood]
    myLen = tDUM[colBlobLength][bGood]
    uX = myLen * np.cos(np.radians(myOrientation))
    vY = myLen * np.sin(np.radians(myOrientation))

    dumArrows = ax4.quiver(myAz, myAlt, uX, vY)

#Vector Plot Test
#x=MappedAltAz.Az
#y=MappedAltAz.Alt
#u=tDUM['THETA_IMAGE'][bGood]
#fig3 = plt.figure(3)
#fig3.suptitle('Vector plot')
#xLabel3=plt.xlabel('Theta-')
#ax3 = fig3.add_subplot(111)

#x=np.array([9 ,10])
#y=np.array([11, 12])
#u=np.array([3.14])
#v=np.array([2])
#vectors=ax3.quiver([x,y],u,v)
#plt.show()
#fits.close()
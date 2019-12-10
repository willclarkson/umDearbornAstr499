
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

    #print(boundsX)
    #print(boundsY)
        # Now that we've read in the header for the image, use it to
        # convert pixel coords to sky coords if they didn't come in with
        # the catalog. Reference:
    
     # https://docs.astropy.org/en/stable/wcs/
 
            
    fig=plt.figure()
    fig.suptitle('Object Locations in Azimuth and Altitude')
    fig.clf()
    yLabel6=plt.ylabel('Altitude')
    xLabel6=plt.xlabel('Azimuth')
    ax6=fig.add_subplot(111)
    UMD_Observatory= EarthLocation(lat=41.32*u.deg, lon=-83.24*u.deg )
       #GEOMETRY OF THE END POINTS OF EACH ELIPSE
    EL=tDUM['A_IMAGE']*[bGood]
    ELVector=np.array(EL)
    arL=len(EL)
    Divisor=[2]*arL
    Len2=np.divide(ELVector,Divisor)
    ET=tDUM['THETA_IMAGE']*[bGood]
    Thet=np.array(ET)
    cosT=np.cos(Thet)
    sinT=np.sin(Thet)
    AS=np.multiply(Len2,sinT)
    AC=np.multiply(Len2,cosT)
    ECX=tDUM['X_IMAGE']*[bGood]
    ECY=tDUM['Y_IMAGE']*[bGood]
    CX=np.array(ECX)
    CY=np.array(ECY)
    
    #PLOTTING ENDPOINTS
    EX1=(CX+AC)
    EY1=(CY+AS)
    EX2=(CX-AC)
    EY2=(CY-AS)
    #ENDPOINT1 in ALTAZ
    RA1, DEC1 = wcs.all_pix2world(EX1, EY1, 0)
    tDUM['END_POINT1_RA'] = Column(RA1, unit=u.deg)
    tDUM['END_POINT1_DEC']= Column(DEC1, unit=u.deg)
    objPos1 = SkyCoord(ra= tDUM['END_POINT1_RA'], dec=tDUM['END_POINT1_DEC'], frame='fk5')    
    MappedAltAz1=objPos1.transform_to(AltAz(obstime=time,location=UMD_Observatory))
    myAz1 = np.asarray(MappedAltAz1.az)
    myAlt1 = np.asarray(MappedAltAz1.alt)
    #ENDPOINT2 in ALTAZ
    RA2, DEC2 = wcs.all_pix2world(EX2, EY2, 0)
    tDUM['END_POINT2_RA'] = Column(RA2, unit=u.deg)
    tDUM['END_POINT2_DEC']= Column(DEC2, unit=u.deg)
    objPos2 = SkyCoord(ra= tDUM['END_POINT2_RA'], dec=tDUM['END_POINT2_DEC'], frame='fk5')
    MappedAltAz2=objPos2.transform_to(AltAz(obstime=time,location=UMD_Observatory))
    myAz2 = np.asarray(MappedAltAz2.az)
    myAlt2 = np.asarray(MappedAltAz2.alt)

    dum6 = ax6.scatter(myAz1,myAlt1,marker='.')
    dum66 = ax6.scatter(myAz2,myAlt2,marker='.')

            #Testing Connectors

    AZ12=np.hstack((myAz1,myAz2))
    Con=[] 
    for i in range(0, len(AZ12), 1):
        if i % 2 ==0:
            Con.append(myAz1[i]) 
        else:
            Con.append(myAz2[i])  

    print(Con)
    ALT12=np.hstack((myAlt1,myAlt2))
    Con1=[] 
    for i in range(0, len(ALT12),1):
        if i % 2 ==0:
            Con1.append(myAlt1[i]) 
        else:
            Con1.append(myAlt2[i])  
    print(Con1)
    def connectpoints(x,y,p1,p2):
        x1, x2 = x[p1], x[p2]
        y1, y2 = y[p1], y[p2]
        ax6.plot([x1,x2],[y1,y2],'k-')

    for i in np.arange(len(Con)):
        connectpoints(Con,Con1,i,i+1)
    plt.show()

        #plt.show()
    #https://matplotlib.org/3.1.1/api/backend_bases_api.html#matplotlib.backend_bases.KeyEvent
    #
    #x, y =
    #for i in range(0, len(x), 2):
     #   plt.plot(x[i:i+2], y[i:i+2])
    #plt.show()

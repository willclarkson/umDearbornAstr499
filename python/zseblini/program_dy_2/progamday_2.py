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

    
    # Now that we've read in the header for the image, use it to
    # convert pixel coords to sky coords if they didn't come in with
    # the catalog. Reference:
    
    # https://docs.astropy.org/en/stable/wcs/

    # read in the header anyway
    if not 'ALPHA_J2000' in tDUM.colnames:
        #wcs = WCS(myHeader)
        xPix = tDUM['X_IMAGE']
        yPix = tDUM['Y_IMAGE']
        RA, DEC = wcs.all_pix2world(xPix, yPix, 0)
        # now that we've calculated RA, DEC, let's add them to the
        # table in-memory so that the rest of the routine can carry on
        # as if they'd come in with the data:
        tDUM['ALPHA_J2000'] = Column(RA, unit=u.deg)
        tDUM['DELTA_J2000'] = Column(DEC, unit=u.deg)
        
        #print tDUM['ALPHA_J2000'][0:5]
        
    UMD_Observatory= EarthLocation(lat=41.32*u.deg, lon=-83.24*u.deg )
    objPos = SkyCoord(ra=tDUM['ALPHA_J2000'][bGood], dec=tDUM['DELTA_J2000'][bGood], frame='fk5')

    MappedAltAz=objPos.transform_to(AltAz(obstime=time,location=UMD_Observatory))
    myAz = np.asarray(MappedAltAz.az)
    myAlt = np.asarray(MappedAltAz.alt)
    #print(myAlt)
    #print(myAz)

    #Scatter plot of detected objects
    def fig4(fig):
        fig.suptitle('Object Locations in Azimuth and Altitude')
        fig.clf()
        yLabel4=plt.ylabel('Altitude')
        xLabel4=plt.xlabel('Azimuth')
        ax4 = fig.add_subplot(111)
        dum4 = ax4.scatter(myAz,myAlt,marker='.')

         # let's try a quiver plot. We will assemble the x- and y- components of the vectors
        myOrientation = tDUM['THETA_IMAGE'][bGood]
        myLen = tDUM[colBlobLength][bGood]
        uX = myLen * np.cos(np.radians(myOrientation))
        vY = myLen * np.sin(np.radians(myOrientation))
        dumArrows = ax4.quiver(myAz, myAlt, uX, vY)

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
    def fig5(fig):
        fig.suptitle('Object End points')
        fig.clf()
        ax5=fig.add_subplot(111)
        dum5 = ax5.scatter(EX1,EY1,marker='.')
        dum55=ax5.scatter(EX2,EY2,marker='.')
        dum555= ax5.plot(boundsX,boundsY)

    #From IMAGE to RA and DEC to ALT AND AZ
    
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

    #Plotting the points in ALTAZ

    def fig6(fig):
        fig.suptitle('Object Locations in Azimuth and Altitude')
        fig.clf()
        yLabel6=plt.ylabel('Altitude')
        xLabel6=plt.xlabel('Azimuth')
        ax6=fig.add_subplot(111)
        dum6 = ax6.scatter(myAz1,myAlt1,marker='.')
        dum66 = ax6.scatter(myAz2,myAlt2,marker='.')

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
    #print(FnX)
    #print(FnY)
    #print(SnX)
    #print(SnY)
    #print(boundsX2)
    #print(boundsY2)


    #SWITCH FIGS
    switch_figs =  {
        0: fig4,
        1: fig5,
        2: fig6
    }


    def onclick1(fig):
        global i
        print(i)
        fig.clear()
        i += 1
        i %= 3
        switch_figs[i](fig)
        plt.draw()

#x = np.linspace(0, 2*np.pi, 1000)
    fig = plt.figure()
    switch_figs[0](fig)
    fig.canvas.mpl_connect('button_press_event', lambda event: onclick1(fig))

    plt.show()




    
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

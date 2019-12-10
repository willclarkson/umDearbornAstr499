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
    Con=[]
    Con1=[] 
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
    Dx=(EX1-EX2)
    Dy=(EY1-EY2)
    def fig5(fig):
        fig.suptitle('Object End points')
        fig.clf()
        ax5=fig.add_subplot(111)

        
        dum5 = ax5.scatter(EX1,EY1,marker='.')
        dum55=ax5.scatter(EX2,EY2,marker='.')
        dum555=ax5.arrow(EX1,EY1,Dx,Dy)
        dum5555= ax5.plot(boundsX,boundsY)

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
    DAZ=np.array(myAz1-myAz2)
    DALT=np.array(myAlt1-myAlt2)
    #Plotting the points in ALTAZ
    
    def fig6(fig):
        fig.suptitle('Object Locations in Azimuth and Altitude')
        fig.clf()
        yLabel6=plt.ylabel('Altitude')
        xLabel6=plt.xlabel('Azimuth')
        ax6=fig.add_subplot(111)
        #dum6 = ax6.scatter(myAz1,myAlt1,marker='.')
        #dum66 = ax6.scatter(myAz2,myAlt2,marker='.')
        #dum6666=ax6.arrow(myAz1,myAlt1,DAZ,DALT)
        #dum6666.all()

        #Finding the border in ALTAZ
        UMD_Observatory= EarthLocation(lat=41.32*u.deg, lon=-83.24*u.deg )
        FnY = myHeader['NAXIS1']
        FnX = myHeader['NAXIS2']
        SnXRA,SnYDEC =wcs.all_pix2world(FnX, FnY, 0)
        objPos = SkyCoord(ra=SnXRA, dec=SnYDEC, unit='deg', frame='fk5')
        MappedAltAz3=objPos.transform_to(AltAz(obstime=time,location=UMD_Observatory))
        myAz3 = np.asarray(MappedAltAz3.az)
        myAlt3 = np.asarray(MappedAltAz3.alt)
        Zeros1=2
        Zeros=2
        ZerosRA,ZerosDEC=wcs.all_pix2world(Zeros1,Zeros,0)
        objPosZeros = SkyCoord(ra=ZerosRA, dec=ZerosDEC, unit='deg', frame='fk5')
        MappedAltAz0=objPosZeros.transform_to(AltAz(obstime=time,location=UMD_Observatory))
        myAz0 = np.asarray(MappedAltAz0.az)
        myAlt0 = np.asarray(MappedAltAz0.alt)
        print(myAz0)
        print(myAlt0)
        boundsX0 = np.array([myAz0, myAz3, myAz3, myAz0], 'float')
        boundsY0 = np.array([myAlt0,myAlt0, myAlt3, myAlt3], 'float')
        boundsX2 = np.hstack((boundsX0, boundsX0[0]))
        boundsY2 = np.hstack((boundsY0, boundsY0[0]))
        print(boundsX0)
        print(boundsY0)
        dum666= ax6.plot(boundsX2,boundsY2)Pr
        
        #CONNECTIONS
        AZ12=np.hstack((myAz1,myAz2))
        
        for i in range(0, len(AZ12), 1):
            if i % 2 ==0:
                Con.append(myAz1[i]) 
            else:
                Con.append(myAz2[i])  
        ALT12=np.hstack((myAlt1,myAlt2))
        
        for i in range(0, len(ALT12),1):
            if i % 2 ==0:
                Con1.append(myAlt1[i]) 
            else:
                Con1.append(myAlt2[i])  
        def connectpoints(x,y,p1,p2):
            x1, x2 = x[p1], x[p2]
            y1, y2 = y[p1], y[p2]
            ax6.plot([x1,x2],[y1,y2],'k-')
        for i in np.arange(len(Con)):
            connectpoints(Con,Con1,i,i+1)

    def fig7(fig):
        fig.suptitle('Object End points')
        fig.clf()
        ax7=fig.add_subplot(111)
        UMD_Observatory= EarthLocation(lat=41.32*u.deg, lon=-83.24*u.deg )
        FnY = myHeader['NAXIS1']
        FnX = myHeader['NAXIS2']
        SnXRA,SnYDEC =wcs.all_pix2world(FnX, FnY, 0)
        objPos = SkyCoord(ra=SnXRA, dec=SnYDEC, unit='deg', frame='fk5')
        MappedAltAz3=objPos.transform_to(AltAz(obstime=time,location=UMD_Observatory))
        myAz3 = np.asarray(MappedAltAz3.az)
        myAlt3 = np.asarray(MappedAltAz3.alt)
        Zeros1=2
        Zeros=2
        ZerosRA,ZerosDEC=wcs.all_pix2world(Zeros1,Zeros,0)
        objPosZeros = SkyCoord(ra=ZerosRA, dec=ZerosDEC, unit='deg', frame='fk5')
        MappedAltAz0=objPosZeros.transform_to(AltAz(obstime=time,location=UMD_Observatory))
        myAz0 = np.asarray(MappedAltAz0.az)
        myAlt0 = np.asarray(MappedAltAz.alt)
        boundsX0r = np.asarray([ZerosRA, SnXRA, SnXRA, ZerosRA], 'float')
        boundsY0r = np.asarray([ZerosDEC, ZerosDEC, SnYDEC, SnYDEC], 'float')
        boundsX2r = np.hstack((boundsX0r, boundsX0r[0]))
        boundsY2r = np.hstack((boundsY0r, boundsY0r[0]))
        print(boundsX2r)
        print(boundsY2r)
        
        
        dum7 = ax7.scatter(RA1, DEC1,marker='.')
        dum77=ax7.scatter(RA2, DEC2,marker='.')
        dum777= ax7.plot(boundsX2r,boundsY2r)
        #dum777=ax7.arrow(EX1,EY1,Dx,Dy)
        #dum7777= ax7.plot(ZerosRA,ZerosDEC)
    
        
    #print(FnX)
    #print(FnY)
    #print(SnX)
    #print(SnY)
    #print(boundsX2)
    #print(boundsY2)


    #SWITCH FIGS
#    switch_figs =  {
 #       0: fig6,
  #      1: fig5,
   #     2: fig6
    #}


    #def onclick1(fig):
     #   global i
      #  print(i)
       # fig.clear()
    #    i += 1
     #   i %= 3
      #  switch_figs[i](fig)
       # plt.draw()

    
    fig = plt.figure()
    plt.draw()
    fig6(fig)
    fig7(fig)
    #switch_figs[0](fig)
    #fig.canvas.mpl_connect('button_press_event', lambda event: onclick1(fig))


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

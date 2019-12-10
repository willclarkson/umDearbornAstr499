#Boundary points testing
from astropy.table import Table, Column
import matplotlib.pylab as plt
from astropy.io import ascii
from astropy.io import fits
#import seaborn as sns
from astropy.stats import sigma_clip
import numpy as np
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
 
    def fig6(fig):
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

            
            
            #plt.show()



            #Finding the border in ALTAZ
    
            FnY = myHeader['NAXIS1']
            FnX = myHeader['NAXIS2']
            SnXRA,SnYDEC =wcs.all_pix2world(FnX, FnY, 0)
            Zeros1=1
            Zeros=1
            ZerosRA,ZerosDEC=wcs.all_pix2world(Zeros1,Zeros,0)
            #boundsX1 = np.asarray([0., SnX, SnX, 0.], 'float')
            #boundsY1 = np.asarray([0., 0., SnY, SnY], 'float')
            #boundsX2 = np.hstack((SnX,SnX[0]))
            #boundsY2 = np.hstack((SnY,SnY[0]))
            #boundsX2 = np.hstack((SnX, SnX[0]))
            #boundsY2 = np.hstack((SnY, SnY[0]))
            #objPos = SkyCoord(ra=boundsX2, dec=boundsX2, frame='fk5')
            #tDUM['TnX'] = Column(SnX, unit=u.deg)
            #tDUM['TnY'] = Column(SnY, unit=u.deg)
            objPos = SkyCoord(ra=SnXRA, dec=SnYDEC, unit='deg', frame='fk5')
            objPosZeros = SkyCoord(ra=ZerosRA, dec=ZerosDEC, unit='deg', frame='fk5')
            MappedAltAz4=objPosZeros.transform_to(AltAz(obstime=time,location=UMD_Observatory))
            UMD_Observatory= EarthLocation(lat=41.32*u.deg, lon=-83.24*u.deg )
            MappedAltAz3=objPos.transform_to(AltAz(obstime=time,location=UMD_Observatory))
            myAz3 = np.asarray(MappedAltAz3.az)
            myAlt3 = np.asarray(MappedAltAz3.alt)
            myAzzeros = np.asarray(MappedAltAz4.az)
            myAltzeros = np.asarray(MappedAltAz4.alt)
            #print(myAzzeros)
            #print(myAltzeros)
            print(FnX)
            print(FnY)
            print(SnXRA)
            print(SnYDEC)
            print(myAz3)
            print(myAlt3)
            #boundsX0 = np.asarray([0.1, MappedAltAz3.az, MappedAltAz3.az, 0.1], 'float')
            #boundsY0 = np.asarray([0.1, 0.1, MappedAltAz3.alt, MappedAltAz3.alt], 'float')
            boundsX0 = np.asarray([myAzzeros, myAz3, myAz3,myAzzeros], 'float')
            boundsY0 = np.asarray([myAltzeros,myAltzeros, myAlt3, myAlt3], 'float')
            print(boundsX0)
            print(boundsY0)
            #boundsX2 = np.hstack((myAz3,myAz3[0]))
            #boundsY2 = np.hstack((myAlt3,myAlt3[0]))
            boundsX2 = np.hstack((boundsX0, boundsX0[0]))
            boundsY2 = np.hstack((boundsY0, boundsY0[0]))
            dum666= ax6.plot(boundsX2,boundsY2)

           
            #dum666= ax6.plot(boundsX2,boundsY2)
            #print(boundsX2)
            #print(boundsY2)
    def fig7(fig):
            fig.suptitle('Object Locations in Image Co-Ords')
            fig.clf()
            yLabel7=plt.ylabel('Y-Image')
            xLabel7=plt.xlabel('X-Image')
            ax7=fig.add_subplot(111)
            FnY = myHeader['NAXIS1']
            FnX = myHeader['NAXIS2']
            boundsX0 = np.asarray([0.1, FnX, FnX, 0.1], 'float')
            boundsY0 = np.asarray([0.1, 0.1, FnY, FnY], 'float')
            boundsX3 = np.hstack((boundsX0,boundsX0[0]))
            boundsY3 = np.hstack((boundsY0,boundsY0[0]))
            dum7= ax7.plot(boundsX3,boundsY3)


            print(FnX)
            print(FnY)


    def fig8(fig):
            fig.suptitle('Object Locations in Image Co-Ords')
            fig.clf()
            yLabel8=plt.ylabel('Y-Image')
            xLabel8=plt.xlabel('X-Image')
            ax8=fig.add_subplot(111)
            FnY = myHeader['NAXIS1']
            FnX = myHeader['NAXIS2']
            TnX,TnY =wcs.all_pix2world(FnX, FnY, 0)
            boundsX4 = np.asarray([0, TnX, TnX, 0], 'float')
            boundsY4 = np.asarray([0, 0, TnY, TnY], 'float')
            boundsX5 = np.hstack((boundsX4,boundsX4[0]))
            boundsY5 = np.hstack((boundsY4,boundsY4[0]))
            dum8= ax8.plot(boundsX5,boundsY5)
            print(TnX)
            print(TnY)
            print(boundsX5)
            print(boundsY5)


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
    fig.canvas.mpl_connect('button_press_event', lambda event: onclick1(fig))

    plt.show()


   # def on_key(event):
    #    global i
     #   print(i)
      #  fig.clear()
       # i += 1
       # i %= 3
      #  switch_figs[i](fig)
      #  plt.draw()

    
    #fig = plt.figure()
    #switch_figs[0](fig)
    #fig.canvas.mpl_connect('key_press_event', on_key)

    #plt.show()
from astropy.table import Table, Column
import matplotlib.pylab as plt
from astropy.io import ascii
from astropy.io import fits
#import seaborn as sns
from astropy.stats import sigma_clip
from astropy.stats import signal_to_noise_oir_ccd
import numpy as np
import astropy.units as u
#import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, ITRS, AltAz
from matplotlib.pylab import quiver
from numpy import multiply
import matplotlib.animation as animation
from matplotlib.animation import FFMpegWriter
from astropy.utils import iers
from astropy.stats import histogram



iers.Conf.iers_auto_url.set('ftp://cddis.gsfc.nasa.gov/pub/products/iers/finals2000A.all')

# for timing execution
import time as systemTime

# For world coordinate system
from astropy.wcs import WCS
Tzero = systemTime.time()

i=0
plt.ioff()
UMD_Observatory= EarthLocation(lat=41.32*u.deg, lon=-83.24*u.deg )


def go(fCat='GaiaCatalog0.ASC', \
       fHeader='V404_Cyg_adOFF-012_R_120sTest_MAPPED.fit', \
       colFilteredData='FilteredData',colBlobLength='A_IMAGE', blobLenDefault=5., blobSF=7., \
       frameName=''):
    
    """Plotting the Co-ordinates in a Click-Animated Sequence"""

    # for ''safety'' for the figures, close everything plotted here
    plt.close('all')
    
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
    SciImg = fits.open(fHeader)
    Img=SciImg[0].data
    #myHeader = fits.getheader('V404_Cyg_adOFF-012_R_120sTest_MAPPED.fit')
    myHeader = SciImg[0].header #fits.getheader(fHeader)
    time=myHeader['DATE-OBS']


    # Look at the Fits Header to find the Alt Az of where the telescope thinks
    # its pointing
    

    # Let's promote the world coordinate system parsing out to here,
    # since we're going to need it whatever we do
    wcs = WCS(myHeader)


    #Define the good data
    xCut = 4000
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
                
    #objPos = SkyCoord(ra=tDUM['ALPHA_J2000'][bGood], dec=tDUM['DELTA_J2000'][bGood], frame='fk5')
    objPosSTARS = SkyCoord(ra=tDUM['ALPHA_J2000'][bGood], dec=tDUM['DELTA_J2000'][bGood], frame='fk5')

    MappedAltAz=objPosSTARS.transform_to(AltAz(obstime=time,location=UMD_Observatory))
    myAz = np.asarray(MappedAltAz.az)
    myAlt = np.asarray(MappedAltAz.alt)
    
    # Replicate zahra's computation of the head of the arrow, 
    # but vectorially and only for the [bgood] items
    deltaX = tDUM['A_IMAGE']*np.cos(tDUM['THETA_IMAGE'])*blobSF
    deltaY = tDUM['A_IMAGE']*np.sin(tDUM['THETA_IMAGE'])*blobSF

    xTail = tDUM['X_IMAGE']
    yTail = tDUM['Y_IMAGE']


    xHead = xTail + deltaX
    yHead = yTail + deltaY

    # now we have our positions - in IMAGE space - for the head of the arrow. 
    # Use our new-found expertise in transformation to map these onto the sky...
    tailRA,tailDEC =wcs.all_pix2world(xTail, yTail, 0)
    objTailSTARS =  SkyCoord(ra=tailRA, dec=tailDEC, frame='fk5', unit='deg')
    tailAltAz=objTailSTARS.transform_to(AltAz(obstime=time,location=UMD_Observatory))

    # do the same for the heads of the vectors
    headRA,headDEC =wcs.all_pix2world(xHead, yHead, 0)
    objHeadSTARS =  SkyCoord(ra=headRA, dec=headDEC, frame='fk5', unit='deg')
    headAltAz=objHeadSTARS.transform_to(AltAz(obstime=time,location=UMD_Observatory))


    # let's extract the variables we want
    myAzHead = np.asarray(headAltAz.az)
    myAltHead = np.asarray(headAltAz.alt)

    myAzTail = np.asarray(tailAltAz.az)
    myAltTail = np.asarray(tailAltAz.alt)
    
    # now let's populate our data table with them
    dAz = (myAzHead - myAzTail)/blobSF
    dAl = (myAltHead - myAltTail)/blobSF
    tDUM['A_ALTAZ'] = np.sqrt(dAz**2 + dAl**2)*u.deg
    tDUM['THETA_ALTAZ'] = np.degrees(np.arctan(dAl/dAz))*u.deg

    
    thisFits=fHeader
    Name=thisFits.split(".")[0]

    print("DEBUG -- input catalog name:",Name)
    
    # Looking at the Objects that we want to Trace
   
    Brightest= np.amax(tDUM['FLUX_ISO'] )
    bBright=np.equal(tDUM['FLUX_ISO'],  Brightest)

    #p
    #return
    Tracking_ObjectX=xTail[bBright]
    Tracking_ObjectY=yTail[bBright]
    #print(Tracking_ObjectX)
    #print(Tracking_ObjectY)
    trackRA,trackDEC = wcs.all_pix2world(Tracking_ObjectX, Tracking_ObjectY, 0)
    objTrackingSTARS =  SkyCoord(ra=trackRA, dec=trackDEC, frame='fk5', unit='deg')
    c_ITRS = objTrackingSTARS.transform_to(ITRS(obstime=time))
    # Calculate local apparent Hour Angle (HA), wrap at 0/24h
    local_ha = UMD_Observatory.lon - c_ITRS.spherical.lon
    local_ha.wrap_at(24*u.hourangle, inplace=True)
    # Calculate local apparent Declination
    local_dec = c_ITRS.spherical.lat
    
    #print("Local apparent HA, Dec={} {}".format(local_ha.to_string(unit=u.hourangle, sep=':'), local_dec.to_string(unit=u.deg, sep=':', alwayssign=True) ))
    TrackingStarsAltAz=objTrackingSTARS.transform_to(AltAz(obstime=time,location=UMD_Observatory))
    trackAz=TrackingStarsAltAz.az
    trackAlt=TrackingStarsAltAz.alt
    print("Tracking Star DEBUG:", trackAlt.degree, trackAz.degree)
    #print(trackAlt)
    #print(trackAz)

    #return

    
    # Statistics Gathering
    # SNR of the image
    ImgMean=np.mean(Img)
    ImgStd=np.std(Img)
    #print(ImgMean)
    #print(ImgStd)
    SNR=np.log10(ImgMean/ImgStd)
    #print(SNR)
    # Using AstroPyStats to find the SNR a diff way
    #ExposureT=myHeader['EXPTIME']
    
    
    #SNRAstr=signal_to_noise_oir_ccd(t=ExposureT,source_eps, sky_eps, dark_eps, rd, npix=9)


    
    #return
    

    # initialize an empty table
    TableOfStats = Table()


    
    # let's specify a column name:
    sCol = 'FLUX_ISO'

    for sCol in ['FLUX_ISO','A_IMAGE','THETA_IMAGE', 'A_ALTAZ','THETA_ALTAZ']:
    
        FilteredData1= sigma_clip(tDUM[sCol][bGood],sigma=4, iters=5)
        filtered_data=np.array(FilteredData1)
        ###print(filtered_data)
        fdmean=np.mean(filtered_data)
        fdStDev=np.std(filtered_data)
        #statsVector=np.array([fdmean,fdStDev])
        #print("FilteredData %3 " % (filtered_data))
        ###print("FilteredData Mean %.3e " % (fdmean))
        ###print("FilteredData Standard Deviation %.3e " % ( fdStDev))
        # StatsColNames=[sCol + '_mu',sCol + '_std']
        #TableOfStats=Table(statsVector*tDUM['FLUX_ISO'].unit,names=StatsColNames)

        # now let's just stick in each column to the existing table:
        TableOfStats[sCol+'_mu'] = [fdmean]*tDUM[sCol].unit
        TableOfStats[sCol+'_std'] = [fdStDev]*tDUM[sCol].unit
    
   # Few things to Add to the Table after the processed Statistics
    TableOfStats['Target_Object_Alt']= trackAlt.degree 
    TableOfStats['Target_Object_Az']= trackAz.degree
    TableOfStats['SNR']=(SNR)*u.dB    

    # let's call fig6 here.
    #SaveUnder=Name + "figure" +".jpg"
    #fig6(SaveUnder, figNum=6, closeAllFigsFirst=True)
   
        
    #print tDUM[0:3]
    def fig5(fig):
        yLabel5=plt.ylabel('Filtered Data')
        #xLabel5=plt.xlabel('Azimut')
        ax5=fig.add_subplot(111)
        dum5=histogram(filtered_data,bins='blocks')
        dum55=histogram(filtered_data,bins='knuth')
        dum555=histogram(filtered_data,bins='scott')
        dum5555=histogram(filtered_data,bins='freedman')
        #dum55=ax5.plot(dum5)
        #print(dum5)
        #print(dum55)
        #print(dum555)
        #print(dum5555)
        
    def fig6(imgName='blah.png', figNum=6, closeAllFigsFirst=True):

        """Updated figure 6 to generate the figure within this method"""

        if closeAllFigsFirst:
            plt.close('all')
            
        fig = plt.figure(figNum)
        ax6=fig.add_subplot(111)
        ax6.set_xlabel('Azimuth')
        ax6.set_ylabel('Altitude')

        # 2019-12-18 - commented out the hardcoding of the limits for now
        #ax6.set_xlim(111., 113.)
        #ax6.set_ylim(29., 30.)
        
        #fig.suptitle('Object Locations in Azimuth and Altitude')
        #yLabel6=plt.ylabel('Altitude')
        #xLabel6=plt.xlabel('Azimuth')
        #ax6=fig.add_subplot(111)
        plt.axis([111,125,27,45])
        

        
        FnY = np.float(myHeader['NAXIS1'])
        FnX = np.float(myHeader['NAXIS2'])
        SnXRA,SnYDEC =wcs.all_pix2world(FnX, FnY, 0)
        objPos = SkyCoord(ra=SnXRA, dec=SnYDEC, unit='deg', frame='fk5')

        boundsSkyra,boundsSkydec=wcs.all_pix2world(boundsX, boundsY, 0)
        boundsSky=SkyCoord(ra=boundsSkyra, dec=boundsSkydec, unit='deg', frame='fk5')
        #print("fig6 DEBUG:", objPosSTARS.ra.degree, boundsSky.ra.degree)
        MappedAltAzSky=boundsSky.transform_to(AltAz(obstime=time,location=UMD_Observatory))
        myAzSky = np.asarray(MappedAltAzSky.az)
        myAltSky = np.asarray(MappedAltAzSky.alt)
        MappedAltAz3=objPos.transform_to(AltAz(obstime=time,location=UMD_Observatory))
        myAz3 = np.asarray(MappedAltAz3.az)
        myAlt3 = np.asarray(MappedAltAz3.alt)
        dum66666= ax6.plot(myAzSky,myAltSky)
        Zeros1=2
        Zeros=2
        ZerosRA,ZerosDEC=wcs.all_pix2world(Zeros1,Zeros,0)
        objPosZeros = SkyCoord(ra=ZerosRA, dec=ZerosDEC, unit='deg', frame='fk5')
        MappedAltAz0=objPosZeros.transform_to(AltAz(obstime=time,location=UMD_Observatory))
        myAz0 = np.asarray(MappedAltAz0.az)
        myAlt0 = np.asarray(MappedAltAz0.alt)
        boundsX0 = np.array([myAz0, myAz3, myAz3, myAz0], 'float')
        boundsY0 = np.array([myAlt0, myAlt3,myAlt0, myAlt3], 'float')
        boundsX2 = np.hstack((boundsX0, boundsX0[0]))
        boundsY2 = np.hstack((boundsY0, boundsY0[0]))
        trackingstars=ax6.plot(trackAz,trackAlt, markersize=5,color='r')

        #dum666= ax6.plot(boundsX2,boundsY2)

        myAzHeadTrim = myAzHead[bGood]
        myAzTailTrim = myAzTail[bGood]
        myAltHeadTrim = myAltHead[bGood]
        myAltTailTrim = myAltTail[bGood]
        
        
        # Let's try building a quiver plot out of the start and end points
        dAz = myAzHeadTrim - myAzTailTrim
        dAlt = myAltHeadTrim - myAltTailTrim

        # let's try adding a quiver plot in data units now:
        #quiv6 = ax6.quiver(myAzTail, myAltTail, dAlt, dAz, units='xy', angles='uv', color='r')

        # get out of the function before tripping the part that currently breaks
        #return

        # let's try a less clever method to plot up the [1,2] pairs:

        #print("DEBUG - azTail, azHead", np.shape(myAzTail), np.shape(myAzHead))
        for iRow in range(np.size(myAzTailTrim)):
            blah = ax6.plot([myAzTailTrim[iRow], myAzHeadTrim[iRow]], \
                            [myAltTailTrim[iRow],myAltHeadTrim[iRow]]\
                , color='c')
            blah = ax6.plot(myAzTailTrim[iRow], myAltTailTrim[iRow]\
                , color='c', marker='o', markersize=2)
        #lala=ax6.plot(myHeader['AZ_OBJ'],myHeader['ALT_OBJ'], marker='o', markersize=5)
                

        ax6.set_title(Name)

        fig.show()
        fig.savefig(imgName)
        #fig = plt.figure()
        #plt.draw()
        #fig6(fig)
        #plt.show()

        #CONNECTIONS

   # print("INFO - time to execute RA->AltAz: %.3e seconds" % (systemTime.time()-tZero))


    # 2019-12-18 updated this to accept an input argument
    if len(frameName) < 3:
           SaveUnder=Name + "figure" +".jpg"
    else:
           SaveUnder = frameName[:] # copy rather than reference
    fig6(SaveUnder, figNum=6, closeAllFigsFirst=True)
    #plt.close('all')

    ##fig = plt.figure(6)
    ##fig.clf()
    
    #plt.draw()
    ##fig6(fig, imgName=SaveUnder)
    #plt.show()
    #fig.savefig(SaveUnder)
    #fig5(fig)
    #plt.show()
    return TableOfStats

    # File: slides.py

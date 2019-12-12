.. _supported_observatories:

Supported Observatories
~~~~~~~~~~~~~~~~~~~~~~~

The pipeline is currently set up to work on data from the following
observatories/instruments:

+--------------------------+--------------------+----------------+
| Telescope                | Instrument         | PP Keyword     |
+==========================+====================+================+
| Apache Point ARC 3.5m    | AGILE              | ARC35AGILE     |
+--------------------------+--------------------+----------------+
| Apache Point ARC 3.5m    | ARCTIC             | ARC35ARCTIC    |
+--------------------------+--------------------+----------------+
| Apache Point ARC 3.5m    | SPICAM             | ARC35SPICAM    |
+--------------------------+--------------------+----------------+
| Calar Alto 1.23m         | DLR-MkIII          | CA123DLRMKIII  |
+--------------------------+--------------------+----------------+
| CTIO 0.9m                | CFCCD              | CTIO09         |
+--------------------------+--------------------+----------------+
| CTIO 1.0m                | Y4KCam             | CTIO10         |
+--------------------------+--------------------+----------------+
| CTIO 1.3m                | ANDICAM (CCD)      | CTIO13CCD      |
+--------------------------+--------------------+----------------+
| Danish 1.54m             | DFOSC (imaging)    | DFOSC          |
+--------------------------+--------------------+----------------+
| Discovery Channel        | Large Monolithic   | DCTLMI         |
| Telescope                | Imager             |                |
+--------------------------+--------------------+----------------+
| Gemini North 8.2m        | GMOS               | GMOSN          |
+--------------------------+--------------------+----------------+
| Gran Canaria 1.5m        | MUSCAT2            | TCS15MUSCAT2   |
+--------------------------+--------------------+----------------+
| INT                      | WFC                | INTWFC         |
+--------------------------+--------------------+----------------+
| IRSF 1.4m                | SIRIUS             | IRSFSIRIUS     |
+--------------------------+--------------------+----------------+
| KMTnet SAAO              | --- (*)            | KMTNETS        |
+--------------------------+--------------------+----------------+
| KPNO 4m                  | Mosaic 1           | KPNO4MOS1      |
+--------------------------+--------------------+----------------+
| KPNO 4m                  | Mosaic 3           | KPNO4MOS3      |
+--------------------------+--------------------+----------------+
| KPNO 4m                  | NEWFIRM            | KPNO4NEWF      |
+--------------------------+--------------------+----------------+
| LCOGT 1m (LSC)           | SBIG (kb78)        | LCOSBIGKB78    |
+--------------------------+--------------------+----------------+
| LCOGT 1m (LSC)           | Sinistro (fl03)    | LCOSINFL03     |
+--------------------------+--------------------+----------------+
| LCOGT 1m (CPT)           | Sinistro (fl06)    | LCOSINFL06     |
+--------------------------+--------------------+----------------+
| LCOGT 1m (LSC)           | Sinistro (fa03)    | LCOSINFA03     |
+--------------------------+--------------------+----------------+
| LCOGT XXXm (XXX)         | Sinistro (fa15)    | LCOSINFA15     |
+--------------------------+--------------------+----------------+
| LCOGT 1m (COJ)           | Sinistro (fl11)    | LCOSINFL11     |
+--------------------------+--------------------+----------------+
| LCOGT 1m (CPT)           | Sinistro (fl16)    | LCOSINFL16     |
+--------------------------+--------------------+----------------+
| LCOGT 2m (COJ)           | Spectral (fs01)    | LCOSPECFS01    |
+--------------------------+--------------------+----------------+
| Lowell 31"               | NASACAM            | LOWELL31       |
+--------------------------+--------------------+----------------+
| Lowell 42"               | NASA42             | LOWELL42       |
+--------------------------+--------------------+----------------+
| Lowell 42"               | SITE               | LOWELL42SITE   |
+--------------------------+--------------------+----------------+
| Lowell 24"               | LONEOS             | LONEOS         |
+--------------------------+--------------------+----------------+
| Magellan                 | IMACS long camera  | MAGIMACSL      |
+--------------------------+--------------------+----------------+
| Magellan                 | IMACS short camera | MAGIMACSS      |
+--------------------------+--------------------+----------------+
| McDonald Struve 2.1m     | CQUEAN             | STRUVECQUEAN   |
+--------------------------+--------------------+----------------+
| MMT                      | MMTCAM             | MMTCAM         |
+--------------------------+--------------------+----------------+
| NEXT 0.6m                | FLI                | NEXT           |
+--------------------------+--------------------+----------------+
| Nordic Optical           | ALFOSC             | NOTALFOSC      |
| Telescope                |                    |                |
+--------------------------+--------------------+----------------+
| Observatoire Haute-      | CCD                | OHP120         |
| Provence 1.2m            |                    |                |
+--------------------------+--------------------+----------------+
| Palomar 60-inch          | optical facility   | P60OPT         |
|                          | camera             |                |
+--------------------------+--------------------+----------------+
| Palomar 60-inch          | SED Machine        | P60SEDM        |
+--------------------------+--------------------+----------------+
| San Pedro Martir 84cm    | Mexman (E2V CCD)   | MEXMAN         |
+--------------------------+--------------------+----------------+
| Steward Observatory 90"  | Spacewatch Camera  | STEWARD90SCC   |
| Bok Telescope            |                    |                |
+--------------------------+--------------------+----------------+
| SOAR 4.1m                | Goodman (**)       | SOARGOODMAN    |
+--------------------------+--------------------+----------------+
| Sutherland 40-inch       | SHA                | SL40IN         |
+--------------------------+--------------------+----------------+
| Telescopio Nazionale     | DOLORES            | TNGDOLORES     |
| Galileo                  |                    |                |
+--------------------------+--------------------+----------------+
| UKIRT                    | WFCAM              | UKIRTWFCAM     |
+--------------------------+--------------------+----------------+
| Vatican Advanced         | VATT4k             | VATT4K         |
| Technology Telescope     |                    |                |
+--------------------------+--------------------+----------------+
| Very Large Telescope     | FORS2              | VLTFORS2       |
+--------------------------+--------------------+----------------+
| WIYN 0.9m                | Half Degree Imager | WIYN09HDI      |
+--------------------------+--------------------+----------------+
| ZTF                      | Mosaic (***)       | ZTFMOSAIC      |
+--------------------------+--------------------+----------------+
| Generic Telescope        | any                | GENERIC        |
+--------------------------+--------------------+----------------+

(*): This camera is a multi-detector instrument; it is recommended to
split multi-extension FITS frames from this instrument into individual
single-extension FITS images and to run the pipeline on these
individual FITS images.

(**): SOAR Goodman image header have a number of non-standard keywords
that force astropy to crash; in order to use PP, please remove header
keywords ``PARAM0``, ``PARAM61``, ``PARAM62``, and ``PARAM63`` prior
to running PP, e.g., using the ``delhead`` command provided by
`WCSTools`_.

(***): wcs provided with telescope data is usually excellent; in this case
I suggest skipping registration and to use ``pp_run`` with the
``-keep_wcs`` option.

It is recommended to stitch images from cameras with multiple amplifiers
together using the correct orientations. Furthermore, it is recommended to
flatten FITS files with multiple headers by simply combining the individual
headers under a single extension.
 
If you would like to use the pipeline for other observatories, please
contact me.

.. _supported_catalogs:

Supported Reference Catalogs 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PP is currently able to use the following catalogs for astrometric and
photometric calibration:

+------------------------+--------------------------+---------------+--------------------------+------------------------------------------------------------+
| Catalog Name           | Catalog Type             | Registration? | Filter Bands             | Comments                                                   |
+========================+==========================+===============+==========================+============================================================+
| Gaia DR1 (`Gaia`_)     | astrometric              | yes           | G                        | all-sky catalog, excellent astrometry, G < 20 mag          |
+------------------------+--------------------------+---------------+--------------------------+------------------------------------------------------------+
| Gaia/TGAS (`TGAS`_)    | astrometric              | yes           | G                        | all-sky catalog, excellent astrometry, G < 12 mag          |
+------------------------+--------------------------+---------------+--------------------------+------------------------------------------------------------+
| USNO-B1.0              | astrometric              | yes           | R2                       | good astrometry, V < 21 mag                                |
+------------------------+--------------------------+---------------+--------------------------+------------------------------------------------------------+
| 2MASS (`2MASS`_)       | astrometric/photometric  | yes           | J, H, Ks, K* (Vega)      | all-sky NIR catalog, good astrometry                       |
+------------------------+--------------------------+---------------+--------------------------+------------------------------------------------------------+
| URAT-1 (`URAT-1`_)     | astrometric/photometric  | yes (SCAMP    | g, r, i (SDSS AB);       | good coverage over the Northern hemisphere, photometry from|
|                        |                          | >= trunk.r345)| B, V, R*, I* (Vega)      | APASS (see below)                                          |
+------------------------+--------------------------+---------------+--------------------------+------------------------------------------------------------+
| Sloan Digital Sky      | astrometric/photometric  | yes           | u, g, r, i, z (SDSS AB); | excellent photometry, Northern hemisphere, patchy coverage |
| Survey Release 9       |                          |               | U*, B*, V*, R*, I* (Vega)|                                                            | 
| (`SDSS-R9`_)           |                          |               |                          |                                                            |
+------------------------+--------------------------+---------------+--------------------------+------------------------------------------------------------+
| AAVSO Photometric All  | photometric              | no            | g, r, i (SDSS AB);       | good coverage, good photometry for stars with V<17         | 
| Sky Survey Release 9   |                          |               | B, V, R*, I* (Vega)      |                                                            |
| (`APASS9`_)            |                          |               |                          |                                                            |
+------------------------+--------------------------+---------------+--------------------------+------------------------------------------------------------+
| Pan-STARRS DR1         | photometric              | no            | g, r, i, z, y (SDSS AB); | good coverage, good photometry for stars with V<20;        | 
| (`PANSTARRS`_)         |                          |               | B*, V*, R*, I* (Vega)    | currently, only cone searches with radius < 0.5 deg        |
|                        |                          |               |                          | supported                                                  |
+------------------------+--------------------------+---------------+--------------------------+------------------------------------------------------------+
| SkyMapper DR1          | photometric              | no            | u, v, g, r, i, z         | good coverage in the Southern hemisphere                   | 
| (`SkyMapper`_)         |                          |               | (SDSS AB)                |                                                            |
|                        |                          |               |                          |                                                            |
+------------------------+--------------------------+---------------+--------------------------+------------------------------------------------------------+

The catalog name in brackets is the identifier used by PP; e.g., if
you want to use URAT-1 for the registration of your images, use option
`-cat URAT-1` with the ``pp_register`` command. The "Registration?"
column identifies if this catalog is supported by SCAMP in order to
use it for image registration. All filter bands marked with an
asterisk (*) are obtained through transformations from other bands;
the respective photometric system is shown in brackets.


If you are interested in using catalogs other than those listed,
please let me know.



.. _supported filters:

Supported Catalog Transformations
---------------------------------

PP supports the following catalog transformations:

* ugriz -> BVRI: `Chonis & Gaskell 2008`_
* JHKs (2MASS) -> JHK (UKIRT): `Hodgkin et al. 2009`_
* Pan-STARRS grizy -> SDSS griz + BVRI: `Tonry et al. 2012`_
  
Independent checks indicate that these transformations are reliable and accurate. More quantitative results coming soon...


.. _Chonis & Gaskell 2008: http://adsabs.harvard.edu/abs/2008AJ....135..264C
.. _Hodgkin et al. 2009: http://adsabs.harvard.edu/abs/2009MNRAS.394..675H
.. _Tonry et al. 2012: http://adsabs.harvard.edu/abs/2012ApJ...750...99T

.. _Gaia: http://sci.esa.int/gaia/
.. _TGAS: https://gaia.esac.esa.int/documentation/GDR1/Data_processing/chap_cu3tyc/
.. _2MASS: http://www.ipac.caltech.edu/2mass/
.. _URAT-1: http://cdsads.u-strasbg.fr/cgi-bin/nph-bib_query?2015AJ....150..101Z&db_key=AST&nosetcookie=1
.. _SDSS-R9: http://www.sdss3.org/dr9/
.. _APASS9: http://www.aavso.org/apass
.. _PANSTARRS: http://panstarrs.stsci.edu/
.. _SkyMapper: http://skymapper.anu.edu.au


.. _WCSTools: http://tdc-www.harvard.edu/wcstools/

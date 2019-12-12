Changelog
=========

Major changes to the pipeline since 2016-10-01 (see `Mommert 2017`_) are
documented here.

* 2018-12-02: major overhaul of diagnostic output

* 2018-11-23: implementation of ``pp_setup.py``, which will eventually
  replace ``_pp_conf.py``; photometric catalog data used in the photometric
  calibration can now by output as ascii table and/or into final .db file;
  ``pillow`` has been replaced by ``skimage``; ``pandas`` is now a required
  Python module

* 2017-11-24: implementation of `-solar` option in ``pp_run`` and
  ``pp_calibrate`` to obtain photometric calibration only from stars
  with Sun-like colors

* 2017-10-20: implementation of ``pp_stackedphotometry``, providing
  automated image stacking and subsequent photometry

* 2017-10-04: implementation of ``pp_combine``, which enables
  automated image combination

* 2017-08-29: implementation of Gaia/TGAS as an alternative
  astrometric catalog for shallow widefield observations

* 2017-06-01: extraction of serendipitously observed targets
  (asteroids and variable stars) implemented in ``pp_distill``

* 2017-03-19: if no filtername is provided (``None``) or the
  ``-instrumental`` option of ``pp_calibrate`` is used, ``pp_run``
  will complete all pipeline tasks using these instrumental magnitudes

* 2017-03-16: implementation of `Pan-STARRS DR1`_ for the photometric
  calibration; currently, PP uses a home-built access of MAST at
  STScI, which is limited to catalog queries with a maximum cone
  radius of 0.5 deg; please note that this kind of query is rather
  slow compared to Vizier queries of SDSS or APASS

* 2017-02-24: ``pillow`` is now a required python module; the pipeline
  now supports default distortion parameters for wide-field cameras

* 2017-02-04: catalog.data is now an astropy.table, catalog downloads
  using astroquery.vizier (no effect to the user), pp now supports
  Gaia DR1 (CMC, USNOB1, and PPMX have been removed)


  
.. _Mommert 2017: http://adsabs.harvard.edu/abs/2017A%26C....18...47M
.. _Pan-STARRS DR1: http://panstarrs.stsci.edu/




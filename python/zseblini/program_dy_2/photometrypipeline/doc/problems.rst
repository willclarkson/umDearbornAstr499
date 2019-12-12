Problems?
=========

Frequently Asked Questions
--------------------------

General Problems
~~~~~~~~~~~~~~~~

**The pipeline does not run and gives me an error similar to `pp_run: Command not found`. What is wrong?**
   There are two things to check: (1) did you properly setup the
   `PHOTPIPEDIR` environment variable (check this by typing ``echo
   $PHOTPIPEDIR`` in a terminal), or (2) the command ``pp_run`` uses a
   symbolic link to ``pp_run.py``; if the former does not work, try
   the latter.
   
**Running the pipeline creates warnings (e.g., FutureWarning); what
  should I do?**
   Warnings are usually caused by Python module version issues and
   cause no harm whatsoever. In theory, there should be no warnings,
   but they happen occasionally. If you are haunted by some warnings,
   let me know and I will try to resolve the issue.

**I get an error message like `Intel MKL FATAL ERROR: Cannot load
 libmkl_avx2.so or libmkl_def.so.`. What do I have to do?**
   This seems to be a problems with your numpy or scipy
   installation. Try to install the latest versions of both
   packages. If you are using anaconda, try `conda install -f numpy`
   and `conda install -f scipy`.
   
  
pp_calibrate (Photometric Calibration)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**I keep getting output like `zeropoint for image.fits: Warning: 0 reference stars after source matching for frame image.ldac`. What does it mean?**
   It means that none of the reference stars with measured magnitudes
   in your field of view could be matched with a source in your
   image. As a result, the magnitudes in the photometry output files
   are simply instrumental magnitudes, not calibrated ones. Try using
   a different photometric catalog in :func:`pp_calibrate`. If your
   field of view is small (<2 arcmin), there might just be no stars
   with known magnitudes in the field, in which case there is not a
   lot that can be done...

**The pipeline fails to derive useful magnitude zeropoints from 2MASS data**
   This problem might be solved by installing the latest version of astropy
   (currently 2.0.2).
   
**The pipeline crashes in ``pp_calibrate`` with the following error
message: `IndexError: too many indices for array`.
** Well, this is embarassing... I am familiar with this problem, but I
haven't found a way to solve it, yet. The problem is that Source
Extractor runs in ``pp_photometry`` to produce an array of aperture
photometry results for the curve-of-growth analysis. After that,
Source Extractor is supposed to run then once again using the optimum
aperture radius. Sometimes, this second Source Extractor runs seems to
fail, causing ``pp_calibrate`` to fail. A manual workaround is to run
``pp_photometry`` again using the ``-aprad`` option and the optimum
aperture derived from the last run. Running ``pp_calibrate`` again
will then succeed. - **Update: this problem should now be fixed. If
you still encounter problems, please let me know!**

   
pp_distill (Target Photometry Extraction)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Why does the target photometry vary wildly, although the magnitude zeropoints are consistent with each other?**

   Check the offset of the photometry position to the expected target
   position. If the offsets are not consistent, it is likely that PP
   picks up noise or nearby sources (also check the target thumbnail
   images). This can happen if the target is very faint (in this case,
   lower the `snr` limit in :func:`pp_photometry`), or if the `snr`
   limit and `minarea` paramters are picked too small and too much
   noise is picked up. In either case, redo the :func:`pp_photometry`
   step and play with the `snr` and `minarea` parameters.
   

... This Does Not Solve My Problem...
-------------------------------------

If you encounter problems, e.g., PP stops unexpectedly or crashes,
please don't hesitate to contact me (michael . mommert (at) nau . edu)
and attach the following things to your email:

* the **LOG file** of your PP process (can be found in the
  ``.diagnostics/`` directory), and 

* the **error message** that was printed on the screen.

I will try to get back to you as soon as possible.

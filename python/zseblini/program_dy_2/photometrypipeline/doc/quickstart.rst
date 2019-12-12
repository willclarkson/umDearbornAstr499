.. _quickstart:

Quickstart
==========


Prerequisits
------------

Image data should be properly reduced before using the pipeline for
best results, including cropping the data section. Bias subtraction
and flat fielding improves photometry results but is not absolutely
necessary. PP's ability to provide astrometric and photometric
calibration puts some constraints on the way data is stored: data from
separate fields, as well as data using different instrument settings
(e.g., different binning modes) should be stored in individual
directories, which in turn should be separated by filters::

  -+- all_data -+
                |
                +- field_1 -+- filter_1
                |           +- filter_2
                |           +- filter_3
                |
                +- field_2 -+- filter_1
                |           +- filter_2
               ...

Separate fields are defined as having gaps between individual frames
that are comparable to, or larger than, the field of view. Series of
frames that were tracked on a moving target can be put in the same
directory if the total track is smaller than 3-5 times the size of a
single-frame field of view.

Moving objects are currently only identified based on the images'
``OBJECT`` keywords. The object name should be as simple as possible,
consisting either of the bodies official number or designation; please
use either a blank or an underscore to separate the designation's year
from the identifier.


Running PP
----------

PP can be run in a fully automated or semi-automated mode, providing
different levels of user interaction.


Fully Automated Mode
~~~~~~~~~~~~~~~~~~~~

In the directory tree example above, PP can be run in different
places, treating the data differently. If you want to run PP only on
data for one field and filter, you can change in that directory and
**run PP locally on all fits files in that directory**::

  cd all_data/field_1/filter_1
  pp_run *fits

If your data are organized as shown in the example above, you can **run
PP from any higher level directory to analyze all underlying
directories in a consecutive way**::
  
  cd all_data
  pp_run all

Passing ``all`` signalizes PP to walk through underlying directories,
starting from the current one. What happens is that PP creates a PP
subprocess for each data directory. In case you want PP to only run on
a subset of fits files starting with a certain prefix, you can use
option `-prefix`, e.g., ::

  pp_run -prefix reduced all 

is the equivalent of using only files that are included in
``reduced*.fits``.


Semi-Automated Mode Walkthrough
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This section describes the individual steps PP takes to analyze the
data. While ``pp_run`` performs these steps automatically, each of the
following functions can be called manually, which allows to tweak the
analysis process. If you intend to perform the analysis fully
manually, please note the individual functions have to be called in
the following order:

* :func:`pp_prepare`: prepare the input images and implant rough
  WCS information into the image header

* :func:`pp_register`: use ``SCAMP`` to register all input images
  based on the implanted rough WCS information; different catalogs are
  tried until all images have been registered ; this function calls
  :func:`pp_extract` automatically

* :func:`pp_photometry`: derive instrumental magnitudes using a
  curve-of-growth analysis, or a manually provided aperture radius

* :func:`pp_calibrate`: photometrically calibrate instrumental
  magnitudes and create a `SQLite` database file for each image; note
  that this function has to be called even if you plan on using
  instrumental magnitudes only (use the `-instrumental` option)

* :func:`pp_distill`: extract target information from the photometry
  databases created by the previous task; see the function reference
  for the different options of target identification

.. _manual target identification:

Manual Target Identification
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In case the target has no identifier, or positions/ephemerides cannot
be obtained automatically (e.g., space debris, newly discovered
asteroids, etc.), or you want to verify the calibration accuracy using
a manually selected control star in the field, the target has to be
identified manually from the image data using :func:`pp_manident`.

Image data are at minimum required to have passed :func:`pp_prepare`,
:func:`pp_photometry`, and :func:`pp_calibrate`; :func:`pp_manident`
may also be called after a full :func:`pp_run` call. In order to
identify the target in all images, :func:`pp_manident` allows you to
browse through all images and click on the target. The trajectory of
the target is fit using a spline function. Quitting
:func:`pp_manident` creates a ``positions.dat`` file, which can be
used as input for :func:`pp_distill` using the `-positions` option.

The manual target identification also allows the user to extract
photometry from images with highly trailed background stars. In that
case, the resulting photometry will consist of instrumental
magnitudes. Hence, :func:`pp_register` does not have to be called and
:func:`pp_calibrate` should be called using the `-instrumental`
option. Positions used in the target identification and listed in the
final photometry file are based on the rough WCS information implanted
by :func:`pp_prepare` and should not be trusted!


PP Diagnostics
--------------

PP generates by default significant amounts of diagnostic information
on each run. These information can be accessed in the individual
directories where the data resides with any web browser, e.g., ::

  firefox all_data/field_2/filter_3/diagnostics.html

If you ran PP with the `all` argument (see above), a file
``summary.html`` will be generated in the root directory (``all_data``),
which provides links to the individual ``index.html`` files.


More information on the diagnostic output is available here:
:ref:`diagnostics`.


Results
-------

PP derives the calibrated photometry for the target that it finds in
the ``OBJECT`` header keyword, as well as one rather bright 'control
star' that is used to check the consistency of the photometric
calibration. Results are written to files
``photometry_<objectname>.dat`` in the respective filter directory.


Although PP is designed to run mostly automatically, some common sense
is required to make sure the results are reliable. 

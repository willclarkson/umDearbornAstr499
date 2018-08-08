
import sys
import os
import fileinput
import glob

def helloWorld(gain = 2.):

    """Quick script to print to terminal"""

    print "Hello, world -- %.2f" % (gain)
# EDIT 6/14/18 -- Gain was found to be ~4.0
def launchDIA(gain=4.0, forcePhotom=True):

    """Our launcher script for DIA"""

    try:
        from pyDIA import calibration_functions as cal
    except:
        print "launchDIA WARN - problem importing pyDIA"

        return

    print "startDIA launching... using gain: ", gain

    ####################

    # 3/8/18
    # set package install location
    #
    sys.path.append(os.getcwd())

    #  import high-level pipeline routines
    use_GPU = False
    if use_GPU:
        from pyDIA import DIA_GPU as DIA
    else:
        from pyDIA import DIA_CPU as DIA

    from pyDIA import calibration_functions as cal
    from pyDIA.data_structures import Parameters


#params = DIA.parameters() # program was giving error, "name 'params' is not defined"; so maybe defining it before will solve this? -- It did not.

    # loading the default parameters UPDATE 2018-03-19 "Parameters" seems to be in pyDIA.data_structures...
    params = Parameters()
#params = DIA.Parameters()
    params.use_GPU = use_GPU


    # override parameter defaults if necessary

    params.gain = gain # MDM website does not give a default value for this parameter
    params.readnoise = 5.0 # 2.5 MHz according to MDM website for Andor camera; default value was 5.0
    params.pixel_min = 10
    params.pixel_max = 140000 # consider bringing lower if segfault again
    params.name_pattern = '*.fits'
    #params.datekey = 'PHJDMID' # default given by pyDIA documentation; some inexplicable, non-program-breaking errors arise when using this
    params.datekey = 'HJD'
    params.pdeg = 1
    params.sdeg = 1
    params.bdeg = 1
    #params.use_stamps = False
    params.use_stamps = True
    params.nstamps = 9  # if we get segfault again, try bringing this down to 10 (from 100)
    params.loc_output = 'Output_dir'

    print "INFO - about to start DIA processing"

    # Perform difference imaging, photometry, and calibration

#    if not(os.path.exists(params.loc_output) or forcePhotom):
    if forcePhotom:
            print "HERE"
            DIA.imsub_all_fits(params)
            print "INFO - tried imsub_all_fits"
    if not(os.path.exists(params.loc_output+'/calibration.png')):
            cal.calibrate(params.loc_output)

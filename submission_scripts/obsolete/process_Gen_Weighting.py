#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Before submission you need to (or before submitting dag which takes the current environment):
# source /cvmfs/icecube.opensciencegrid.org/users/BeyondStandardModel/bsm-py2-v3.1.1-1/setup.sh
# /cvmfs/icecube.opensciencegrid.org/users/lfischer/oscnext_meta_V01-00-05/build/env-shell.sh

### Script to add the I3MCWeightDict to the generation level files  ###

#
# Author: Leander Fischer
#

import time
t0 = time.time()

from optparse import OptionParser

usage = "usage: %prog [options]"
parser = OptionParser(usage)
parser.add_option("-i", "--infile",
                 dest="INFILE", help="Read input from INFILE (.i3{.gz/.zst} format)")
parser.add_option("-o", "--outfile",
                 dest="OUTFILE", help="Write output to OUTFILE (.i3{.gz/.zst} format)")
parser.add_option("-g", "--gcdfile", default="/cvmfs/icecube.opensciencegrid.org/data/GCD/GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE_withScaledNoise.i3.gz",
                 dest="GCDFILE", help="Read in GCD file")
parser.add_option("-d","--datatype",type="string", default="LeptonInjector",
                  dest="datatype", help="Generator data type (used for weighting).")
parser.add_option("-v", "--verbose", default=False, action="store_true",
                 dest="VERBOSE", help="print info level output")
(options,args) = parser.parse_args()
if not options.INFILE:
    parser.error('Infile not specified.')

datatype = options.datatype

import os
from I3Tray import *
from icecube import dataio, dataclasses
import icecube.icetray.i3logging as logging

if options.VERBOSE:
    logging.set_level('INFO')
else:
    logging.set_level('WARN')

logging.log_info('Options: {}'.format(options))

from icecube.oscNext.frame_objects.weighting import WEIGHT_DICT_KEY
from icecube.oscNext.frame_objects.weighting import weighting

OUTPATH = os.path.split(options.OUTFILE)[0]
if not os.path.exists(OUTPATH):
    logging.log_info('Creating output filepath: {}'.format(OUTPATH))
    os.makedirs(OUTPATH)

##### Run the weighting module #####

tray = I3Tray()
    
tray.context['I3FileStager'] = dataio.get_stagers()

tray.AddModule('I3Reader',
               'reader',
               FilenameList=[options.GCDFILE, options.INFILE],
              )

# Store the weighting info
tray.Add(weighting,
    "weighting",
    data_type=datatype,
    )

# tray.Add( do_something, 'test_segment', streams = [icetray.I3Frame.DAQ, icetray.I3Frame.Physics] )

tray.Add("I3Writer",
    Streams=[icetray.I3Frame.DAQ, icetray.I3Frame.Physics, icetray.I3Frame.TrayInfo, icetray.I3Frame.Simulation],
    filename=options.OUTFILE)


tray.Execute()
tray.Finish()

##### Done #####

t1 = time.time()
logging.log_info('Time it took: {:.3f} s'.format(t1-t0))

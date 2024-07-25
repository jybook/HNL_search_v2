#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/icetray-start
#METAPROJECT /cvmfs/icecube.opensciencegrid.org/users/lfischer/oscnext_hnl_meta/build

#!/usr/bin/env python
# -*- coding: utf-8 -*-

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
parser.add_option("-v", "--verbose", default=False, action="store_true",
                 dest="VERBOSE", help="print info level output")
(options,args) = parser.parse_args()
if not options.INFILE:
    parser.error('Infile not specified.')
if not options.GCDFILE:
    parser.error('GCDfile not specified.')


import os
import numpy as np
from icecube import icetray, dataclasses, dataio
from I3Tray import *
from icecube.dataclasses import I3Double, I3Particle, I3Constants
from icecube.icetray import I3Bool

from icecube import photonics_service

from icecube.millipede import TaupedeFit

import icecube.icetray.i3logging as logging


if options.VERBOSE:
    logging.set_level('INFO')
else:
    logging.set_level('WARN')

logging.log_info('Options: {}'.format(options))

OUTPATH = os.path.split(options.OUTFILE)[0]
if not os.path.exists(OUTPATH):
    logging.log_info('Creating output filepath: {}'.format(OUTPATH))
    os.mkdir(OUTPATH)


"""
Perform reconstruction and pre-processing
"""

# use uncleaned pulses since millipede accounts for noise
pulses='SplitInIcePulses'

tray = I3Tray()
    
tray.context['I3FileStager'] = dataio.get_stagers()

tray.AddModule('I3Reader',
               'reader',
               FilenameList=[options.GCDFILE, options.INFILE],
              )

# add table base/cascade service (newest photonics tables including effective distande parametrization by Marcel Usner)
# needs this photonics-service project http://code.icecube.wisc.edu/svn/projects/photonics-service/branches/ice-properties
table_base = os.path.expandvars('$I3_DATA/photon-tables/splines/cascade_single_spice_lea_flat_z20_a10.{}.fits')
cascade_service = photonics_service.I3PhotoSplineService(amplitudetable = table_base.format('abs'),
                                                         timingtable = table_base.format('prob'),
                                                         effectivedistancetable = os.path.expandvars('$I3_DATA/photon-tables/splines/cascade_effectivedistance_spice_lea_z20.eff.fits'),
                                                         timingSigma = 0,
                                                         maxRadius=400.
                                                        )

# set TaupedeFit parameters
common_params = dict(
    Pulses=pulses,
    CascadePhotonicsService=cascade_service,
    # minimizer='SIMPLEX',
    PhotonsPerBin = 1,
    )

# logging.log_info('common_params: {}'.format(common_params))

name = 'MuMillipedeFit'
Seed = 'TaupedeFit'


# run the fit (with the TaupedeFit as Seed)
tray.AddModule("MuMillipede",
        name,
        Output=name,
        ShowerSpacing=5,
        MuonSpacing=0,
        SeedTrack = Seed,
        **common_params)


tray.Add("I3Writer",
    # Streams=map(icetray.I3Frame.Stream, 'QP'),
    Streams=[icetray.I3Frame.DAQ, icetray.I3Frame.Physics, icetray.I3Frame.TrayInfo, icetray.I3Frame.Simulation],
    DropOrphanStreams=[icetray.I3Frame.DAQ],
    filename=options.OUTFILE)

tray.Execute()
tray.Finish()

t1 = time.time()
logging.log_info('Time it took: {:.3f} s'.format(t1-t0))
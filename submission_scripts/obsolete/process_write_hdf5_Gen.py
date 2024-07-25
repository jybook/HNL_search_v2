#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Before submission you need to (or before submitting dag which takes the current environment):
# source /cvmfs/icecube.opensciencegrid.org/users/BeyondStandardModel/bsm-py2-v3.1.1-1/setup.sh
# /cvmfs/icecube.opensciencegrid.org/users/lfischer/oscnext_meta_V01-00-05/build/env-shell.sh

### Script to store i3 files as hdf5 ###

#
# Author: Leander Fischer
#

import time

t0 = time.time()

from optparse import OptionParser

usage = "usage: %prog [options]"
parser = OptionParser(usage)
parser.add_option(
    "-i",
    "--infile",
    dest="INFILE",
    help="Read input from INFILE (.i3{.gz/.zst} format)",
)
parser.add_option(
    "-o",
    "--outfile",
    dest="OUTFILE",
    help="Write output to OUTFILE (.i3{.gz/.zst} format)",
)
parser.add_option(
    "-g",
    "--gcdfile",
    default="/cvmfs/icecube.opensciencegrid.org/data/GCD/GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE_withScaledNoise.i3.gz",
    dest="GCDFILE",
    help="Read in GCD file",
)
parser.add_option(
    "-v",
    "--verbose",
    default=False,
    action="store_true",
    dest="VERBOSE",
    help="print info level output",
)
(options, args) = parser.parse_args()
if not options.INFILE:
    parser.error("Infile not specified.")

import os
import numpy as np
from I3Tray import *
from icecube import (
    dataio,
    dataclasses,
    icetray,
    phys_services,
    LeptonInjector,
)
import icecube.icetray.i3logging as logging

if options.VERBOSE:
    logging.set_level("INFO")
else:
    logging.set_level("WARN")

logging.log_info("Options: {}".format(options))


from icecube.hdfwriter import I3HDFTableService, I3SimHDFWriter
from icecube.tableio import I3TableWriter

# write to the same directory as the infiles are in
options.OUTFILE = options.INFILE

HDF5_OUTFILE = options.OUTFILE.split(".i3")[0] + ".h5"

hdf5 = I3HDFTableService(HDF5_OUTFILE)

OUTPATH = os.path.split(options.OUTFILE)[0]
if not os.path.exists(OUTPATH):
    logging.log_info("Creating output filepath: {}".format(OUTPATH))
    os.makedirs(OUTPATH)

HDF5_KEYS = [
    # EventProperties
    "mHNL",
    "mHNL_min",
    "mHNL_max",
    "distance",
    "distanceMin",
    "distanceMax",
    "lifetime",
    "totalEnergy",
    "finalStateX",
    "finalStateY",
    # I3MCTree
    # true HNL variables
    "HNL_true_x",
    "HNL_true_y",
    "HNL_true_z",
    "HNL_true_energy",
    "HNL_true_zenith",
    "HNL_true_azimuth",
    "HNL_true_time",
    # true primary variables
    "true_x",
    "true_y",
    "true_z",
    "true_energy",
    "true_zenith",
    "true_azimuth",
    "true_time",
    # true first (DIS) cascade variables
    "casc0_true_x",
    "casc0_true_y",
    "casc0_true_z",
    "casc0_true_energy",
    "casc0_true_zenith",
    "casc0_true_azimuth",
    "casc0_true_time",
    # true second (HNL decay) cascade variables
    "casc1_true_x",
    "casc1_true_y",
    "casc1_true_z",
    "casc1_true_energy",
    "casc1_true_zenith",
    "casc1_true_azimuth",
    "casc1_true_time",
    "nan_decay_energy",
]

# check hdf5 keys
for k in HDF5_KEYS:
    assert isinstance(k, basestring), "`hdf5_keys` must be strings, found %s %s" % (
        k,
        type(k),
    )

# remove duplicates from the HDF5 keys
HDF5_KEYS = list(
    set(HDF5_KEYS)
)  # note that this does NOT preserve order, but doesn't matter


# function to get HNL true energy
def store_mc_true_variables(frame):
    if not frame.Has("I3MCTree"):
        return False
    mctree = frame["I3MCTree"]
    p_true = mctree.primaries[0]

    p_daughters = mctree.get_daughters(p_true)
    assert len(p_daughters) == 2

    for p_daughter in p_daughters:
        if p_daughter.type == dataclasses.I3Particle.Hadrons:
            casc_0_true = p_daughter
        else:
            hnl_true = p_daughter

    hnl_daughters = mctree.get_daughters(hnl_true)
    assert len(hnl_daughters) > 0

    nan_energy = False

    for count_hnl_daughters, hnl_daughter in enumerate(hnl_daughters):
        if not count_hnl_daughters:
            casc_1_true = hnl_daughter
        else:
            assert casc_1_true.pos == hnl_daughter.pos
            if np.isnan(hnl_daughter.energy):
                nan_energy = True
            casc_1_true.energy = casc_1_true.energy + hnl_daughter.energy

    frame["true_x"] = dataclasses.I3Double(p_true.pos.x)
    frame["true_y"] = dataclasses.I3Double(p_true.pos.y)
    frame["true_z"] = dataclasses.I3Double(p_true.pos.z)
    frame["true_energy"] = dataclasses.I3Double(p_true.energy)
    frame["true_zenith"] = dataclasses.I3Double(p_true.dir.zenith)
    frame["true_azimuth"] = dataclasses.I3Double(p_true.dir.azimuth)
    frame["true_time"] = dataclasses.I3Double(p_true.time)

    frame["HNL_true_x"] = dataclasses.I3Double(hnl_true.pos.x)
    frame["HNL_true_y"] = dataclasses.I3Double(hnl_true.pos.y)
    frame["HNL_true_z"] = dataclasses.I3Double(hnl_true.pos.z)
    frame["HNL_true_energy"] = dataclasses.I3Double(hnl_true.energy)
    frame["HNL_true_zenith"] = dataclasses.I3Double(hnl_true.dir.zenith)
    frame["HNL_true_azimuth"] = dataclasses.I3Double(hnl_true.dir.azimuth)
    frame["HNL_true_time"] = dataclasses.I3Double(hnl_true.time)

    frame["casc0_true_x"] = dataclasses.I3Double(casc_0_true.pos.x)
    frame["casc0_true_y"] = dataclasses.I3Double(casc_0_true.pos.y)
    frame["casc0_true_z"] = dataclasses.I3Double(casc_0_true.pos.z)
    frame["casc0_true_energy"] = dataclasses.I3Double(casc_0_true.energy)
    frame["casc0_true_zenith"] = dataclasses.I3Double(casc_0_true.dir.zenith)
    frame["casc0_true_azimuth"] = dataclasses.I3Double(casc_0_true.dir.azimuth)
    frame["casc0_true_time"] = dataclasses.I3Double(casc_0_true.time)

    frame["casc1_true_x"] = dataclasses.I3Double(casc_1_true.pos.x)
    frame["casc1_true_y"] = dataclasses.I3Double(casc_1_true.pos.y)
    frame["casc1_true_z"] = dataclasses.I3Double(casc_1_true.pos.z)
    frame["casc1_true_energy"] = dataclasses.I3Double(casc_1_true.energy)
    frame["casc1_true_zenith"] = dataclasses.I3Double(casc_1_true.dir.zenith)
    frame["casc1_true_azimuth"] = dataclasses.I3Double(casc_1_true.dir.azimuth)
    frame["casc1_true_time"] = dataclasses.I3Double(casc_1_true.time)

    frame["nan_decay_energy"] = icetray.I3Bool(nan_energy)

    return True


# function to get EventProperties (LI)
def store_LI_event_properties(frame):
    if not frame.Has("EventProperties"):
        return False
    event_properties = frame["EventProperties"]

    frame["mHNL"] = dataclasses.I3Double(event_properties.mHNL)
    frame["mHNL_min"] = dataclasses.I3Double(event_properties.mHNL_min)
    frame["mHNL_max"] = dataclasses.I3Double(event_properties.mHNL_max)
    frame["distance"] = dataclasses.I3Double(event_properties.distance)
    frame["distanceMin"] = dataclasses.I3Double(event_properties.distanceMin)
    frame["distanceMax"] = dataclasses.I3Double(event_properties.distanceMax)
    frame["lifetime"] = dataclasses.I3Double(event_properties.lifetime)
    frame["totalEnergy"] = dataclasses.I3Double(event_properties.totalEnergy)
    frame["finalStateX"] = dataclasses.I3Double(event_properties.finalStateX)
    frame["finalStateY"] = dataclasses.I3Double(event_properties.finalStateY)
    return True


##### Run the module #####

tray = I3Tray()

tray.context["I3FileStager"] = dataio.get_stagers()

tray.AddModule(
    "I3Reader",
    "reader",
    FilenameList=[options.GCDFILE, options.INFILE],
)

# get hnl true energy
tray.AddModule(
    store_mc_true_variables,
    "store_mc_true_variables",
    Streams=[icetray.I3Frame.DAQ],
)

# get event properties (LeptonInjector)
tray.AddModule(
    store_LI_event_properties,
    "store_LI_event_properties",
    Streams=[icetray.I3Frame.DAQ],
)

tray.AddModule(
    "I3NullSplitter",
    "p_frame_writer",
    SubEventStreamName="fullevent",
)

# write to HDF5 file
tray.AddModule(
    I3TableWriter,
    "hdf5_file_writer",
    TableService=hdf5,
    keys=HDF5_KEYS,
    SubEventStreams=["fullevent"],
)

tray.Execute()
tray.Finish()

##### Done #####

t1 = time.time()
logging.log_info("Time it took: {:.3f} s".format(t1 - t0))

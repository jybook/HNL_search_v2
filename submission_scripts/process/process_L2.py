#!/bin/sh /cvmfs/icecube.opensciencegrid.org/users/BeyondStandardModel/bsm-py2-v3.1.1-1/icetray-start-standard
#METAPROJECT /data/user/lfischer/software/oscnext_hnl/build/


### The Level 2 simulation processing script ###

import os
import subprocess, logging

from I3Tray import I3Tray
from icecube import icetray, dataio, dataclasses

import numpy as np

from icecube.filterscripts.offlineL2.level2_all_filters import OfflineFilter
from icecube.filterscripts.offlineL2 import SpecialWriter

from icecube import LeptonInjector  # needed for hdf5 writing
from icecube.hdfwriter import I3HDFWriter


########## Functions needed to store as hdf5 ##########

HDF5_KEYS = [
    # EventProperties
    "decay_channel",
    "distance",
    "distanceMax",
    "distanceMin",
    "finalStateX",
    "finalStateY",
    "final_state_particle0",
    "final_state_particle1",
    "primary_type",
    "lifetime",
    "mHNL",
    "outgoing_neutrino_energy",
    "totalEnergy",
    "physical",
    "total_column_depth",
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
    # weights
    "LeptonInjectorWeight",
    "LifetimeWeight_1e-03",
    "OneWeight",
    "ReferenceWeight_1e-03",
]

# remove duplicates from the HDF5 keys
HDF5_KEYS = list(
    set(HDF5_KEYS)  # note that this does NOT preserve order, but doesn't matter
)

# function to get HNL true variables
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
    frame["distance"] = dataclasses.I3Double(event_properties.distance)
    frame["distanceMin"] = dataclasses.I3Double(event_properties.distanceMin)
    frame["distanceMax"] = dataclasses.I3Double(event_properties.distanceMax)
    frame["lifetime"] = dataclasses.I3Double(event_properties.lifetime)
    frame["totalEnergy"] = dataclasses.I3Double(event_properties.totalEnergy)
    frame["finalStateX"] = dataclasses.I3Double(event_properties.finalStateX)
    frame["finalStateY"] = dataclasses.I3Double(event_properties.finalStateY)
    frame['primary_type'] =  dataclasses.I3Double(event_properties.initialType)
    frame['final_state_particle0'] =  dataclasses.I3Double(event_properties.finalType1)
    frame['final_state_particle1'] =  dataclasses.I3Double(event_properties.finalType2)
    frame['total_column_depth'] =  dataclasses.I3Double(event_properties.totalColumnDepth)
    frame["decay_channel"] = dataclasses.I3Double(event_properties.decay_channel)
    frame["outgoing_neutrino_energy"] = dataclasses.I3Double(event_properties.outgoingNeutrinoEnergy)
    frame["physical"] = dataclasses.I3Double(event_properties.physical)

    return True


# function to get weights
def store_weights(frame):
    if not frame.Has("I3MCWeightDict"):
        return False
    weight_dict = frame["I3MCWeightDict"]

    frame["LeptonInjectorWeight"] = dataclasses.I3Double(weight_dict['LeptonInjectorWeight'])
    frame["LifetimeWeight_1e-03"] = dataclasses.I3Double(weight_dict['LifetimeWeight_1e-03'])
    frame["OneWeight"] = dataclasses.I3Double(weight_dict['OneWeight'])
    frame["ReferenceWeight_1e-03"] = dataclasses.I3Double(weight_dict['ReferenceWeight_1e-03'])

    return True


def make_parser():
    """Make the argument parser"""
    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option(
        "-s",
        "--simulation",
        action="store_true",
        default=True,
        dest="mc",
        help="Mark as simulation (MC) [default: %default]",
    )
    parser.add_option(
        "-i",
        "--input",
        action="store",
        type="string",
        default="test_L1.i3.zst",
        dest="infile",
        help="Input i3 file(s)  (use comma separated list for multiple files) [default: %default]",
    )
    parser.add_option(
        "-g",
        "--gcd",
        action="store",
        type="string",
        default="/cvmfs/icecube.opensciencegrid.org/data/GCD/GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE_withScaledNoise.i3.gz",
        dest="gcdfile",
        help="GCD file for input i3 file [default: %default]",
    )
    parser.add_option(
        "-o",
        "--output",
        action="store",
        type="string",
        default="test_L2.i3.zst",
        dest="outfile",
        help="Output i3 file [default: %default]",
    )
    parser.add_option(
        "--osg",
        action="store",
        type="string",
        default="False",
        dest="osg",
        help="Do you want to run on the OSG?? [default: %default]",
    )
    parser.add_option(
        "-n",
        "--num",
        action="store",
        type="int",
        dest="num",
        help="Number of frames to process [default: %default]",
    )
    parser.add_option(
        "--dstfile",
        action="store",
        type="string",
        default=None,
        dest="dstfile",
        help="DST root file (should be .root) [default: %default]",
    )
    parser.add_option(
        "--gapsfile",
        action="store",
        type="string",
        default=None,
        dest="gapsfile",
        help="gaps text file (should be .txt) [default: %default]",
    )
    parser.add_option(
        "--icetopoutput",
        action="store",
        type="string",
        default=None,
        dest="icetopoutput",
        help="Output IceTop file [default: %default]",
    )
    parser.add_option(
        "--eheoutput",
        action="store",
        type="string",
        default=None,
        dest="eheoutput",
        help="Output EHE i3 file [default: %default]",
    )
    parser.add_option(
        "--slopoutput",
        action="store",
        type="string",
        default=None,
        dest="slopoutput",
        help="Output SLOP file [default: %default]",
    )
    parser.add_option(
        "--rootoutput",
        action="store",
        type="string",
        default=None,
        dest="rootoutput",
        help="Output root file [default: %default]",
    )
    parser.add_option(
        "--photonicsdir",
        action="store",
        type="string",
        default=os.path.expandvars("$I3_DATA/photon-tables"),
        dest="photonicsdir",
        help="Directory with photonics tables [default: %default]",
    )
    parser.add_option(
        "--log-level",
        default="WARN",
        dest="LOG_LEVEL",
        help="Sets the logging level (ERROR, WARN, INFO, DEBUG, TRACE) [default: %default]",
    )
    parser.add_option(
        "--identifier_out",
        type="string",
        default="test",
        dest="identifier_out",
        help="Set name (outfiles). [default: %default]",
    )
    parser.add_option(
        "--hdf5",
        action="store_false",
        dest="write_hdf5",
        default=True,
        help="Write hdf5 file [default: %default]"
    )
    return parser


def main(options, stats={}):
    """The main L2 processing script"""
    tray = I3Tray()

    log_levels = {
        "error": icetray.I3LogLevel.LOG_ERROR,
        "warn": icetray.I3LogLevel.LOG_WARN,
        "info": icetray.I3LogLevel.LOG_INFO,
        "debug": icetray.I3LogLevel.LOG_DEBUG,
        "trace": icetray.I3LogLevel.LOG_TRACE,
    }

    if options["LOG_LEVEL"].lower() in log_levels.keys():
        icetray.set_log_level(log_levels[options["LOG_LEVEL"].lower()])
    else:
        logging.warning("log level option %s not recognized.")
        logging.warning("Options are ERROR, WARN, INFO, DEBUG, and TRACE.")
        logging.warning("Sticking with default of WARN.")
        icetray.set_log_level(icetray.I3LogLevel.LOG_WARN)

    # make list of input files from GCD and infile
    osg = "".join(options["osg"])
    infile = "".join(options["infile"])
    outfile = "".join(options["outfile"])
    gcdfile = "".join(options["gcdfile"])
    if osg == "True":

        def copy_to_OSG(NPX_file):
            subprocess.check_call(
                [
                    "globus-url-copy",
                    "-nodcau",
                    "-rst",
                    NPX_file,
                    "file:" + os.getcwd() + "/" + str(NPX_file.split("/")[-1]),
                ]
            )
            print("copy worked")

        def copy_to_NPX(NPX_file):
            subprocess.check_call(
                [
                    "globus-url-copy",
                    "-nodcau",
                    "-rst",
                    "file:" + os.getcwd() + "/" + str(NPX_file.split("/")[-1]),
                    NPX_file,
                ]
            )

        gcdfile_NPX = str("gsiftp://gridftp.icecube.wisc.edu" + gcdfile)
        gcdfile = str(os.getcwd() + "/" + gcdfile.split("/")[-1])
        infile_NPX = str("gsiftp://gridftp.icecube.wisc.edu" + infile)
        infile = str(os.getcwd() + "/" + infile.split("/")[-1])
        copy_to_OSG(gcdfile_NPX)
        copy_to_OSG(infile_NPX)
        outfile = str("gsiftp://gridftp.icecube.wisc.edu" + outfile)
        outfile_temp = str(os.getcwd() + "/" + outfile.split("/")[-1])
    else:
        temp_dir = os.path.join(
            "/data/ana/BSM/HNL/MC/scripts/temp/", options["identifier_out"]
        )
        if not os.path.exists(temp_dir):
            os.makedirs(temp_dir)
        outfile_temp = os.path.join(temp_dir, outfile.split("/")[-1])

    infiles = [gcdfile, infile]
    logging.info("infiles: ", infiles)

    # test access to input and output files
    for f in infiles:
        if not os.access(f, os.R_OK):
            raise Exception("Cannot read from %s" % f)

    def test_write(f):
        if f:
            try:
                open(f, "w")
            except IOError:
                raise Exception("Cannot write to %s" % f)
            finally:
                os.remove(f)

    test_write(options["outfile"])
    test_write(options["dstfile"])
    test_write(options["gapsfile"])
    test_write(options["icetopoutput"])
    test_write(options["eheoutput"])
    test_write(options["slopoutput"])
    test_write(options["rootoutput"])

    # read input files
    tray.Add(dataio.I3Reader, "Reader", Filenamelist=infiles)

    tray.AddSegment(
        OfflineFilter,
        "OfflineFilter",
        dstfile=options["dstfile"],
        mc=options["mc"],
        doNotQify=options["mc"],
        photonicsdir=options["photonicsdir"],
    )

    ###################################################################
    ########### WRITE STUFF                  ##########################
    ###################################################################

    # Write the physics and DAQ frames
    tray.AddModule(
        "I3Writer",
        "EventWriter",
        filename=outfile_temp,
        Streams=[
            icetray.I3Frame.DAQ,
            icetray.I3Frame.Physics,
            icetray.I3Frame.TrayInfo,
            icetray.I3Frame.Simulation,
        ],
        DropOrphanStreams=[icetray.I3Frame.DAQ],
    )

    # special outputs
    if options["eheoutput"]:
        tray.AddSegment(
            SpecialWriter.EHEWriter, "write_ehe", Filename=options["eheoutput"]
        )
    if options["slopoutput"]:
        tray.AddSegment(
            SpecialWriter.SLOPWriter, "write_slop", Filename=options["slopoutput"]
        )
    if options["rootoutput"]:
        tray.AddSegment(
            SpecialWriter.RootWriter, "write_root", Filename=options["rootoutput"]
        )
    if options["gapsfile"]:
        tray.AddSegment(
            SpecialWriter.GapsWriter,
            "write_gaps",
            Filename=options["gapsfile"],
            MinGapTime=1,
        )
    if options["icetopoutput"]:
        # this needs to be the last output
        tray.AddModule(
            "Delete",
            "deleteicetopextrakeys",
            keys=["InIceRawData", "CleanInIceRawData"],
        )
        tray.AddSegment(
            SpecialWriter.IceTopWriter, "write_icetop", Filename=options["icetopoutput"]
        )

    if options['write_hdf5']:

        outfile_hdf5_temp = outfile_temp.replace('.i3.zst',".hdf5")

        # get hnl true energy
        tray.AddModule(
            store_mc_true_variables,
            "store_mc_true_variables",
            If=lambda frame: frame["I3EventHeader"].sub_event_stream == "InIceSplit",
        )

        # get event properties (LeptonInjector)
        tray.AddModule(
            store_LI_event_properties,
            "store_LI_event_properties",
            If=lambda frame: frame["I3EventHeader"].sub_event_stream == "InIceSplit",
        )

        # get weights
        tray.AddModule(
            store_weights,
            "store_weights",
            If=lambda frame: frame["I3EventHeader"].sub_event_stream == "InIceSplit",
        )

        tray.AddSegment(
            I3HDFWriter,
            output = outfile_hdf5_temp,
            keys = HDF5_KEYS,
            SubEventStreams=["InIceSplit"],
        )

    # make it go
    if options["num"]:
        tray.Execute(options["num"])
    else:
        tray.Execute()

    # print more CPU usage info. than speicifed by default
    tray.PrintUsage(fraction=1.0)
    for entry in tray.Usage():
        stats[entry.key()] = entry.data().usertime

    # clean up forcefully in case we're running this in a loop
    if osg == "True":
        copy_to_NPX(outfile)
        copy_to_NPX(outfile.replace('.i3.zst',".hdf5"))
    else:
        os.system(str("mv {} {}".format(outfile_temp, outfile)))
        os.system(str("mv {} {}".format(outfile_hdf5_temp, outfile.replace('.i3.zst',".hdf5"))))

    del tray


if __name__ == "__main__":
    # run as script from the command line
    # get parsed args
    parser = make_parser()
    (options, args) = parser.parse_args()
    opts = {}
    # convert to dictionary
    for name in parser.defaults:
        value = getattr(options, name)
        if name == "infile" and "," in value:
            value = value.split(",")  # split into multiple inputs
        opts[name] = value

    # call main function
    main(opts)

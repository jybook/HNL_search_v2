#!/bin/sh /cvmfs/icecube.opensciencegrid.org/users/BeyondStandardModel/bsm-py2-v3.1.1-1/icetray-start-standard
#METAPROJECT /n/holylfs05/LABS/arguelles_delgado_lab/Everyone/jbook/i3/build


### The Level 1 simulation processing script ###

import os
import subprocess, logging

from I3Tray import I3Tray
from icecube import icetray, dataio, filter_tools, dataclasses
from icecube import phys_services
from icecube.filterscripts import filter_globals
from icecube.filterscripts.all_filters import OnlineFilter
from icecube.phys_services.which_split import which_split

import numpy as np

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
#     # weights
#     "LeptonInjectorWeight",
#     "LifetimeWeight_1e-03",
#     "OneWeight",
#     "ReferenceWeight_1e-03",
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
        "-i",
        "--input",
        action="store",
        type="string",
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
        default="test_L1.i3.zst",
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
        "--qify",
        action="store_true",
        default=False,
        dest="qify",
        help="Apply QConverter, use if file is P frame only [default: %default]",
    )
    parser.add_option(
        "--MinBiasPrescale",
        action="store",
        type="int",
        default=None,
        dest="MinBiasPrescale",
        help="Set the Min Bias prescale to something other than default [default: %default]",
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
        "--enable-gfu",
        action="store_true",
        default=False,
        dest="GFU",
        help="Do not run GFU filter [default: %default]",
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
    """The main L1 processing script"""
    tray = I3Tray()

    if not options["infile"]:
        raise Exception(
            "You need to specify an input file with the -i or --input parameter."
        )

    if not options["gcdfile"]:
        raise Exception(
            "You need to specify a GCD file with the -g or --gcd parameter."
        )

    if not options["outfile"]:
        raise Exception(
            "You need to specify an output file with the -o or --output parameter."
        )

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
            "/n/holylfs05/LABS/arguelles_delgado_lab/Lab/HNL_MC/temp/", options["identifier_out"]
        )
        if not os.path.exists(temp_dir):
            os.makedirs(temp_dir)
        outfile_temp = os.path.join(temp_dir, outfile.split("/")[-1])

    infiles = [gcdfile, infile]
    # logging.info("infiles: ", infiles)

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

    ########################################

    tray = I3Tray()

    tray.Add(dataio.I3Reader, "reader", filenamelist=infiles)

    # run online filters
    online_kwargs = {}
    if options["photonicsdir"]:
        online_kwargs.update(
            {
                "SplineRecoAmplitudeTable": os.path.join(
                    options["photonicsdir"], "splines", "InfBareMu_mie_abs_z20a10.fits"
                ),
                "SplineRecoTimingTable": os.path.join(
                    options["photonicsdir"], "splines", "InfBareMu_mie_prob_z20a10.fits"
                ),
                #            'alert_followup_base_GCD_filename': options['gcdfile'],
            }
        )
    if options["GFU"] is not None:
        online_kwargs["gfu_enabled"] = options["GFU"]
    tray.AddSegment(
        OnlineFilter,
        "OnlineFilter",
        decode=False,
        simulation=True,
        vemcal_enabled=False,
        alert_followup=False,
        **online_kwargs
    )

    # make random service
    seed = os.getpid()
    filter_mask_randoms = phys_services.I3GSLRandomService(seed)

    # override MinBias Prescale
    filterconfigs = filter_globals.filter_pairs + filter_globals.sdst_pairs
    if options["MinBiasPrescale"]:
        for i, filtertuple in enumerate(filterconfigs):
            if filtertuple[0] == filter_globals.FilterMinBias:
                del filterconfigs[i]
                filterconfigs.append((filtertuple[0], options["MinBiasPrescale"]))
                break
    icetray.logging.log_info("%s" % str(filterconfigs))

    # Generate filter Masks for all P frames
    tray.AddModule(
        filter_tools.FilterMaskMaker,
        "MakeFilterMasks",
        OutputMaskName=filter_globals.filter_mask,
        FilterConfigs=filterconfigs,
        RandomService=filter_mask_randoms,
    )

    # Merge the FilterMasks
    tray.AddModule(
        "OrPframeFilterMasks",
        "make_q_filtermask",
        InputName=filter_globals.filter_mask,
        OutputName=filter_globals.qfilter_mask,
    )

    # Q+P frame specific keep module needs to go first, as KeepFromSubstram
    # will rename things, let's rename post keep.
    def is_Q(frame):
        return frame.Stop == frame.DAQ

    simulation_keeps = [
        "BackgroundI3MCTree",
        "BackgroundI3MCTreePEcounts",
        "BackgroundI3MCPESeriesMap",
        "BackgroundI3MCTree_preMuonProp",
        "BackgroundMMCTrackList",
        "BeaconLaunches",
        "CorsikaInteractionHeight",
        "CorsikaWeightMap",
        "EventProperties",
        "GenerationSpec",
        "I3LinearizedMCTree",
        "I3MCTree",
        "I3MCTreePEcounts",
        "I3MCTree_preMuonProp",
        "I3MCPESeriesMap",
        "I3MCPulseSeriesMap",
        "I3MCPulseSeriesMapParticleIDMap",
        "I3MCWeightDict",
        "LeptonInjectorProperties",
        "MCHitSeriesMap",
        "MCPrimary",
        "MCPrimaryInfo",
        "MMCTrackList",
        "PolyplopiaInfo",
        "PolyplopiaPrimary",
        "RNGState",
        "SignalI3MCPEs",
        "SimTrimmer",  # for SimTrimmer flag
        "TimeShift",  # the time shift amount
        "WIMP_params",  # Wimp-sim
        "noise_weight",  # weights for noise-only vuvuzela simulations
        "I3GENIEResultDict",  # weight informaition for GENIE simulations
    ]

    keep_before_merge = (
        filter_globals.q_frame_keeps
        + [
            "InIceDSTPulses",  # keep DST pulse masks
            "IceTopDSTPulses",
            "CalibratedWaveformRange",  # keep calibration info
            "UncleanedInIcePulsesTimeRange",
            "SplitUncleanedInIcePulses",
            "SplitUncleanedInIcePulsesTimeRange",
            "SplitUncleanedInIceDSTPulsesTimeRange",
            "CalibrationErrata",
            "SaturationWindows",
            "InIceRawData",  # keep raw data for now
            "IceTopRawData",
        ]
        + simulation_keeps
    )

    tray.AddModule("Keep", "keep_before_merge", keys=keep_before_merge, If=is_Q)

    ## second set of prekeeps, conditional on filter content, based on newly created Qfiltermask
    # Determine if we should apply harsh keep for events that failed to pass any filter
    ##  Note: excluding the sdst_streams entries

    tray.AddModule(
        "I3IcePickModule<FilterMaskFilter>",
        "filterMaskCheckAll",
        FilterNameList=filter_globals.filter_streams,
        FilterResultName=filter_globals.qfilter_mask,
        DecisionName="PassedAnyFilter",
        DiscardEvents=False,
        Streams=[icetray.I3Frame.DAQ],
    )

    def do_save_just_superdst(frame):
        if frame.Has("PassedAnyFilter"):
            if not frame["PassedAnyFilter"].value:
                return True  #  <- Event failed to pass any filter.
            else:
                return False  # <- Event passed some filter

        else:
            icetray.logging.log_error("Failed to find key frame Bool!!")
            return False

    keep_only_superdsts = (
        filter_globals.keep_nofilterpass
        + [
            "PassedAnyFilter",
            "InIceDSTPulses",
            "IceTopDSTPulses",
            "SplitUncleanedInIcePulses",
            "SplitUncleanedInIcePulsesTimeRange",
            "SplitUncleanedInIceDSTPulsesTimeRange",
            "RNGState",
        ]
        + simulation_keeps
    )
    tray.AddModule(
        "Keep", "KeepOnlySuperDSTs", keys=keep_only_superdsts, If=do_save_just_superdst
    )

    ## Now clean up the events that not even the SuperDST filters passed on.
    tray.AddModule(
        "I3IcePickModule<FilterMaskFilter>",
        "filterMaskCheckSDST",
        FilterNameList=filter_globals.sdst_streams,
        FilterResultName=filter_globals.qfilter_mask,
        DecisionName="PassedKeepSuperDSTOnly",
        DiscardEvents=False,
        Streams=[icetray.I3Frame.DAQ],
    )

    def dont_save_superdst(frame):
        if frame.Has("PassedKeepSuperDSTOnly") and frame.Has("PassedAnyFilter"):
            if frame["PassedAnyFilter"].value:
                return False  #  <- these passed a regular filter, keeper
            elif not frame["PassedKeepSuperDSTOnly"].value:
                return True  #  <- Event failed to pass SDST filter.
            else:
                return False  # <- Event passed some  SDST filter
        else:
            icetray.logging.log_error("Failed to find key frame Bool!!")
            return False

    tray.AddModule(
        "Keep",
        "KeepOnlyDSTs",
        keys=filter_globals.keep_dst_only
        + ["PassedAnyFilter", "PassedKeepSuperDSTOnly", filter_globals.eventheader],
        If=dont_save_superdst,
    )

    ## Frames should now contain only what is needed.  now flatten, write/send to server
    ## Squish P frames back to single Q frame, one for each split:
    tray.AddModule(
        "KeepFromSubstream",
        "null_stream",
        StreamName=filter_globals.NullSplitter,
        KeepKeys=filter_globals.null_split_keeps,
    )

    in_ice_keeps = (
        filter_globals.inice_split_keeps + filter_globals.onlinel2filter_keeps
    )
    in_ice_keeps = in_ice_keeps + [
        "I3EventHeader",
        "SplitUncleanedInIcePulses",
        "SplitUncleanedInIcePulsesTimeRange",
        "TriggerSplitterLaunchWindow",
        "I3TriggerHierarchy",
        "GCFilter_GCFilterMJD",
    ]
    tray.AddModule(
        "Keep",
        "inice_keeps",
        keys=in_ice_keeps,
        If=which_split(split_name=filter_globals.InIceSplitter),
    )

    tray.AddModule(
        "KeepFromSubstream",
        "icetop_split_stream",
        StreamName=filter_globals.IceTopSplitter,
        KeepKeys=filter_globals.icetop_split_keeps,
    )

    # Apply small keep list (SuperDST/SmallTrig/DST/FilterMask for non-filter passers
    # Remove I3DAQData object for events not passing one of the 'filters_keeping_allraw'
    tray.AddModule(
        "I3IcePickModule<FilterMaskFilter>",
        "filterMaskCheck",
        FilterNameList=filter_globals.filters_keeping_allraw,
        FilterResultName=filter_globals.qfilter_mask,
        DecisionName="PassedConventional",
        DiscardEvents=False,
        Streams=[icetray.I3Frame.DAQ],
    )

    ## Clean out the Raw Data when not passing conventional filter
    def I3RawDataCleaner(frame):
        if not (
            (
                "PassedConventional" in frame
                and frame["PassedConventional"].value == True
            )
            or ("SimTrimmer" in frame and frame["SimTrimmer"].value == True)
        ):
            frame.Delete("InIceRawData")
            frame.Delete("IceTopRawData")

    tray.AddModule(
        I3RawDataCleaner, "CleanErrataForConventional", Streams=[icetray.I3Frame.DAQ]
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

    if not options.infile:
        parser.error("No input file specified")

    # convert to dictionary
    opts = vars(options)
    opts["infile"] = opts["infile"].split(",")

    # call main function
    main(opts)

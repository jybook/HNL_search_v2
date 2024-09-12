#!/bin/sh /cvmfs/icecube.opensciencegrid.org/users/BeyondStandardModel/bsm-py2-v3.1.1-1/icetray-start-standard
#METAPROJECT /n/holylfs05/LABS/arguelles_delgado_lab/Everyone/jbook/i3/build


### Takes an I3 file at any level and extracts all information to HDF5 ###

import os
import subprocess, logging

from I3Tray import I3Tray
from icecube import icetray, dataio, dataclasses, simclasses
from icecube.simclasses import I3CLSimEventStatistics

import numpy as np

from icecube.filterscripts.offlineL2.level2_all_filters import OfflineFilter
from icecube.filterscripts.offlineL2 import SpecialWriter

from icecube import LeptonInjector  # needed for hdf5 writing
from icecube.hdfwriter import I3HDFWriter
from icecube.hdfwriter import I3SimHDFWriter
########## Functions needed to store as hdf5 ##########

HDF5_KEYS_GEN = [
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

HDF5_KEYS_PHOT = [
    "casc0_photons",
    "casc0_photons_at_doms",
    "casc0_photon_weights",
    "casc1_photons",
    "casc1_photons_at_doms",
    "casc1_photon_weights",
]
HDF5_KEYS_DET = [
    "MCPESeriesMap",
    "MCPESeriesMap_withNoise",
#     "triggers_fired",
#     "trigger_sources",
#     "trigger_types",
]
HDF5_KEYS_L1 = ["SplitUncleanedInIcePulses"]
HDF5_KEYS_L2 = []

HDF_KEYS = {"gen": HDF5_KEYS_GEN, "phot": HDF5_KEYS_PHOT, "det": HDF5_KEYS_DET, "L1": HDF5_KEYS_L1, "L2": HDF5_KEYS_L2}
      
# function to get HNL true variables
def store_mc_true_variables(frame):
    print("gen reached")
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

# Add photon level information
def store_photon_information(frame):
    if not frame.Has("clsim_stats"):
        return False
           
    photon_stats=frame["clsim_stats"]     
    mctree = frame["I3MCTree"]
    p_true = mctree.primaries[0]
    p_daughters = mctree.get_daughters(p_true)
    
    assert len(p_daughters) == 2

    for p_daughter in p_daughters:
        if p_daughter.type == dataclasses.I3Particle.Hadrons:
            frame["casc0_photons"] = I3CLSimEventStatistics.GetNumberOfPhotonsGeneratedForParticle(photon_stats, p_daughter)#photon_stats.GetNumberOfPhotonsGeneratedForParticle(p_daughter.type)
            frame["casc0_photons_at_doms"] = photon_stats.GetNumberOfPhotonsAtDOMsForParticle(p_daughter.type)
            frame["casc0_totalweights"]=photon_stats.GetSumOfWeightsPhotonsGeneratedForParticle(p_daughter.type)          
        else:
            hnl_true = p_daughter

    hnl_daughters = mctree.get_daughters(hnl_true)
    assert len(hnl_daughters) > 0
    
    casc1_photons=0
    casc1_photons_at_doms=0
    casc1_photon_weights=0
    for hnl_daughter in hnl_daughters:
        casc1_photons += frame["I3MCTree_clsim"].GetNumberOfPhotonsGeneratedForParticle(hnl_daughter.type)
        casc1_photons_at_doms += frame["I3MCTree_clsim"].GetNumberOfPhotonsAtDOMsForParticle(hnl_daughter.type)
        casc1_photon_weights +=frame["I3MCTree_clsim"].GetSumOfWeightsPhotonsGeneratedForParticle(hnl_daughter.type)
        
    frame["casc1_photons"]= casc1_photons   
    frame["casc1_photons_at_doms"]= casc1_photons_at_doms
    frame["casc1_photon_weights"]= casc1_photons_total_weights

    return True

def store_det_information(frame):
    if not frame.Has("MCPESeriesMap"):
        return False
    
#     frame["MCPESeriesMap"]= simclasses.I3MCPESeriesMap(frame["MCPESeriesMap"])
#     frame["MCPESeriesMap_withNoise"]= simclasses.I3MCPESeriesMap(frame["MCPESeriesMap_withNoise"])
#     print("loading maps")
#     The above aren't necesary because the frame contains them already, so including the keys in the list of HDF5 keys to write should be enough

# this doesn't work because of issues writing vectors to I3FrameObjects
#     triggers = frame["I3Triggers"]
#     trigger_types = []
#     trigger_sources = []
#     for trigger in triggers:
#         if trigger.fired:
#             key = trigger.key
#             trigger_types.append(str(dataclasses.TriggerKey.get_type_string(key, key.type)))
#             trigger_sources.append(str(dataclasses.TriggerKey.get_source_string(key, key.source)))
            
#     frame["trigger_sources"] = trigger_sources#dataclasses.I3VectorString(trigger_sources)
#     frame["trigger_types"] = trigger_types#dataclasses.I3VectorString(trigger_types)

#     triggers = frame["I3Triggers"]
#     triggers_fired=0
#     for trigger in triggers:
#         if trigger.fired:
#             triggers_fired += 1
#     frame["triggers_fired"] = dataclasses.I3Double(triggers_fired)
    
    return True

def store_L1_information(frame):
    if not frame.Has("SplitUncleanedInIcePulses"):
        return False
    
    frame["SplitUncleanedInIcePulses"]= simclasses.I3RecoPulseSeriesMapMask(frame["SplitUncleanedInIcePulses"])


def make_parser():
    """Make the argument parser"""
    from optparse import OptionParser
    
    usage = "usage: %prog [options]"
    parser = OptionParser(usage)
   
    parser.add_option(
        "-l",
        "--processing-level",
        action="store",
        type="string",
        default="gen",
        dest="level",
        help="Processing level (gen, phot, det, L1, L2). [default: %default]",
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
        "-o",
        "--output",
        action="store",
        type="string",
        default="test_L2.hdf5",
        dest="outfile",
        help="Output HDF5 file [default: %default]",
    )
    
    parser.add_option(
        "-n",
        "--num",
        action="store",
        type="int",
        dest="num",
        help="Number of frames to process, if not all [default: %default]",
    )
   
    parser.add_option(
        "--identifier_out",
        type="string",
        default="test",
        dest="identifier_out",
        help="Set name (outfiles). [default: %default]",
    )
  
    return parser


def main(options, stats={}):
    """Write information from higher processing levels to hdf5"""
    tray = I3Tray()

    # make list of input files from GCD and infile
    infile = "".join(options["infile"])
    outfile = "".join(options["outfile"])
    
    temp_dir = os.path.join(
        "/n/holylfs05/LABS/arguelles_delgado_lab/Lab/HNL_MC/temp/", options["identifier_out"]
    )
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)
    outfile_temp = os.path.join(temp_dir, outfile.split("/")[-1])

    infiles = [infile]
    
    # Get appropriate HDF keys for desired level

    level = options["level"]
    HDF5_KEYS = HDF_KEYS[level]


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

    # read input files
    tray.Add(dataio.I3Reader, "Reader", Filenamelist=infiles)

    print(HDF5_KEYS)
    ###################################################################
    ########### WRITE STUFF                  ##########################
    ###################################################################

    if (level == "gen"):
        # get hnl true energy
        tray.AddModule(
            store_mc_true_variables,
            "store_mc_true_variables",
            Streams=[icetray.I3Frame.DAQ]
            )

        # get event properties (LeptonInjector)
        tray.AddModule(
            store_LI_event_properties,
            "store_LI_event_properties",
            Streams=[icetray.I3Frame.DAQ]
            )
    if (level == "phot"):
        tray.AddModule(
            store_photon_information,
            "store_photon_information",
            Streams=[icetray.I3Frame.DAQ]
            )

    if (level == "det"):
        print("det level file")
        tray.AddModule(
            store_det_information,
            "store_det_information",
            Streams=[icetray.I3Frame.DAQ, icetray.I3Frame.Physics]
    #         If=lambda frame: frame["I3EventHeader"].sub_event_stream == "InIceSplit",
            )
        print("module added")

    if (level == "L1"):
        tray.AddModule(
            store_L1_information,
            "store_L1_information",           
            Streams=[icetray.I3Frame.DAQ]
            )


#         # get weights
#         tray.AddModule(
#             store_weights,
#             "store_weights",
#             If=lambda frame: frame["I3EventHeader"].sub_event_stream == "InIceSplit",
#         )

    tray.AddSegment(
        I3SimHDFWriter,
        output = outfile_temp,
        keys = HDF5_KEYS,
#         SubEventStreams=["InIceSplit"],
    )

    # make it go
    if options["num"]:
        tray.Execute(options["num"])
    else:
        tray.Execute()
        print("tray executed")
        
    tray.Finish()
    
    os.system(str("mv {} {}".format(outfile_temp, outfile)))

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

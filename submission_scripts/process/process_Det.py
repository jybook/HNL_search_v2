#!/bin/sh /cvmfs/icecube.opensciencegrid.org/users/BeyondStandardModel/bsm-py2-v3.1.1-1/icetray-start-standard
#METAPROJECT /data/user/lfischer/software/oscnext_hnl/build/


### The Detector Level simulation processing script ###

import time

import os
import numpy as np

from optparse import OptionParser

import subprocess

from I3Tray import I3Tray

from icecube import (
    icetray,
    dataclasses,
    dataio,
    simclasses,
    phys_services,
    sim_services,
    DOMLauncher,
    DomTools,
    genie_icetray,
    clsim,
    trigger_sim,
)
from icecube import LeptonInjector  # needed for hdf5 writing
from icecube.icetray import I3Units

from icecube.hdfwriter import I3SimHDFWriter

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


class RemoveLatePhotons(icetray.I3Module):
    """Removes very late photons from the photon series generated by clsim"""

    def __init__(self, context):
        icetray.I3Module.__init__(self, context)
        self.AddOutBox("OutBox")

        self.AddParameter(
            "InputPhotonSeries", "Name of the photon map", "UnweightedPhotons"
        )
        self.AddParameter(
            "TimeLimit", "Time after which the photons should be discarded", 1e5
        )  # In nanoseconds

    def Configure(self):
        self.pmap_in = self.GetParameter("InputPhotonSeries")
        self.time_limit = self.GetParameter("TimeLimit")
        self.photons_dropped = 0

    def DAQ(self, frame):
        latePhotons = False

        # Check if there are any late photons, raise flag if there are
        for domkey, pseries in frame[self.pmap_in]:
            for photon in pseries:
                if photon.GetTime() > self.time_limit:
                    latePhotons = True
                    break

        # If there are NO late photons, just return the frame without touching it
        if not latePhotons:
            self.PushFrame(frame)
            return True

        # If there were late photons, remove them
        else:
            photonTuples = []
            # print("frame type before: {}".format(type(frame[self.pmap_in])))
            for domkey, pseries in frame[self.pmap_in]:
                newPhotonList = []
                for photon in pseries:
                    if photon.GetTime() <= self.time_limit:
                        newPhotonList.append(photon)
                    else:
                        self.photons_dropped += (
                            1  # Keeping track how many photons I discard
                        )
                # Finished checking the photons in the DOM, pass them to the new list
                if newPhotonList != []:
                    photonTuples.append((domkey, newPhotonList))
            frame.Delete(self.pmap_in)


def BasicHitFilter(frame):
    hits = 0
    if frame.Has("MCPESeriesMap"):
        hits = len(frame.Get("MCPESeriesMap"))
    if hits > 0:
        #        print("has photons")
        return True
    else:
        #       print("does NOT have photons")
        return False


def BasicDOMFilter(frame):
    if frame.Has("InIceRawData"):
        if len(frame["InIceRawData"]) > 0:
            return True
        else:
            return False
    else:
        return False


def main():

    ### START ###

    usage = "usage: %prog [options]"
    parser = OptionParser(usage)
    parser.add_option(
        "-i", "--infile", dest="INFILE", help="Read input from INFILE (.i3{.gz} format)"
    )
    parser.add_option(
        "-o",
        "--outfile",
        default="/data/ana/BSM/HNL/MC/test_files/test_output_Det.i3.zst",
        dest="OUTFILE",
        help="Write output to OUTFILE (.i3{.gz} format) [default: %default]",
    )
    parser.add_option(
        "-r",
        "--runnumber",
        type="string",
        default="1",
        dest="RUNNUMBER",
        help="The run/dataset number for this simulation, is used as seed for random generator [default: %default]",
    )
    parser.add_option(
        "-f",
        "--filenr",
        type="string",
        default="1",
        dest="FILENR",
        help="File number, stream of I3SPRNGRandomService [default: %default]",
    )
    parser.add_option(
        "-g",
        "--gcdfile",
        default="/cvmfs/icecube.opensciencegrid.org/data/GCD/GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE_withScaledNoise.i3.gz",
        dest="GCDFILE",
        help="Read in GCD file [default: %default]",
    )
    parser.add_option(
        "-e",
        "--efficiency",
        type="float",
        default=1.0,
        dest="EFFICIENCY",
        help="DOM Efficiency ... the same as UnshadowedFraction [default: %default]",
    )
    parser.add_option(
        "-m",
        "--icemodel",
        default="spice_3.2.1",
        dest="ICEMODEL",
        type="str",
        help="Should be the same ice model as used for photon propagation (step2) [default: %default]",
    )
    parser.add_option(
        "-n",
        "--noise",
        default="vuvuzela",
        dest="NOISE",
        help="Noise model (vuvuzela/poisson/none) [default: %default]",
    )
    parser.add_option(
        "-l",
        "--holeice",
        default="angsens/as.flasher_p1_0.30_p2_-1",
        dest="HOLEICE",
        help="Pick the hole ice parameterization, corresponds to a file name path relative to $I3_SRC/ice-models/resources/models/ [default: %default]",
    )
    parser.add_option(
        "--osg", type="string", help="Is this job running on the OSG?", default="False"
    )
    parser.add_option(
        "--identifier_out",
        type="string",
        default="test",
        dest="identifier_out",
        help="Set name (outfiles). [default: %default]",
    )

    parser.add_option(
        "--hdf5", action="store_false", dest="write_hdf5", default=True, help="Write hdf5 file [default: %default]"
    )

    (options, args) = parser.parse_args()
    if len(args) != 0:
        crap = "Got undefined options:"
        for a in args:
            crap += a
            crap += " "
        parser.error(crap)

    t0 = time.clock()
    print('Starting ...')

    osg = options.osg
    options.FILENR = int(options.FILENR)
    options.RUNNUMBER = int(options.RUNNUMBER)
    write_hdf5  = options.write_hdf5


    ### START ###

    if osg == "True":

        gcdfile_NPX = str("gsiftp://gridftp.icecube.wisc.edu" + options.GCDFILE)
        gcdfile = str(os.getcwd() + "/" + options.GCDFILE.split("/")[-1])

        infile_NPX = str("gsiftp://gridftp.icecube.wisc.edu" + options.INFILE)
        infile = str(os.getcwd() + "/" + options.INFILE.split("/")[-1])

        copy_to_OSG(gcdfile_NPX)
        copy_to_OSG(infile_NPX)

        infiles = [gcdfile, infile]
        outfile = str("gsiftp://gridftp.icecube.wisc.edu" + options.OUTFILE)
        outfile_temp = str(os.getcwd() + "/" + outfile.split("/")[-1])
    else:
        temp_dir = os.path.join(
            "/data/ana/BSM/HNL/MC/scripts/temp/", options.identifier_out
        )
        if not os.path.exists(temp_dir):
            os.makedirs(temp_dir)
        outfile_temp = os.path.join(temp_dir, str(options.OUTFILE.split("/")[-1]))
        gcdfile = options.GCDFILE
        infiles = [gcdfile, options.INFILE]


    print("Using RUNNR: {}".format(options.RUNNUMBER))
    print("DOM efficiency: {}".format(options.EFFICIENCY))
    print("Using hole ice: {}".format(options.HOLEICE))
    print(
        "Looking for ice model in {}".format(
            os.path.expandvars("$I3_SRC/ice-models/resources/models/")
        )
    )
    if options.ICEMODEL == "" or options.ICEMODEL == None:
        print("No ice model provided. The baseline efficiency can be found in cfg.txt")
        print("of the ice model used for photon propagation. ")
        print("Be very careful! ")
        icemodel_path = None
    else:
        icemodel_path = os.path.expandvars(
            "$I3_SRC/ice-models/resources/models/{}".format(options.ICEMODEL)
        )
        if os.path.isdir(icemodel_path):
            print("Folder with ice model found: ", icemodel_path)
        else:
            print("Error! No ice model with such name found in:")
            print(os.path.expandvars("$I3_SRC/ice-models/resources/models/"))
            exit()
    print("Ice model path: {}".format(icemodel_path))

    # now start the actual processing
    tray = I3Tray()
    tray.AddModule("I3Reader", "reader", FilenameList=infiles)

    # Random service
    # from globals import max_num_files_per_dataset
    max_num_files_per_dataset = 100000
    streamnum = options.FILENR

    tray.AddService("I3SPRNGRandomServiceFactory", "sprngrandom")(
        ("Seed", options.RUNNUMBER),
        ("StreamNum", streamnum),
        ("NStreams", max_num_files_per_dataset),
        ("instatefile", ""),
        ("outstatefile", ""),
    )

    # Now fire up the random number generator with that seed
    randomService = phys_services.I3SPRNGRandomService(
        seed=options.RUNNUMBER, nstreams=max_num_files_per_dataset, streamnum=streamnum
    )


    ####
    ## Remove photons from neutron decay and other processes that take too long (unimportant)
    ####

    tray.AddModule(
        RemoveLatePhotons,
        "RemovePhotons",
        InputPhotonSeries="I3Photons",
        # InputPhotonSeries = "PhotonSeriesMap",
        TimeLimit=1e5,
    )  # nanoseconds


    ####
    ## Make hits from photons (change efficiency here already!)
    ####

    tray.AddModule("I3GeometryDecomposer", "I3ModuleGeoMap")

    gcd_file = dataio.I3File(gcdfile)


    tray.AddSegment(
        clsim.I3CLSimMakeHitsFromPhotons,
        "makeHitsFromPhotons",
        #                MCTreeName="I3MCTree_clsim",
        #                PhotonSeriesName="UnweightedPhotons2",
        PhotonSeriesName="I3Photons",
        MCPESeriesName="MCPESeriesMap",
        RandomService=randomService,
        DOMOversizeFactor=1.0,
        UnshadowedFraction=options.EFFICIENCY,
        IceModelLocation=icemodel_path,
        #               UseHoleIceParameterization=holeice
        HoleIceParameterization=os.path.expandvars(
            "$I3_SRC/ice-models/resources/models/%s" % options.HOLEICE
        ),
        GCDFile=gcd_file,
    )


    # from icecube.BadDomList import bad_dom_list_static
    txtfile = (
        os.path.expandvars("$I3_SRC")
        + "/BadDomList/resources/scripts/bad_data_producing_doms_list.txt"
    )
    tray.AddModule(
        BasicHitFilter,
        "FilterNullMCPE",
        Streams=[icetray.I3Frame.DAQ, icetray.I3Frame.Physics],
    )

    mcpe_to_pmt = "MCPESeriesMap"
    if options.NOISE == "poisson":
        print("Error! Poisson noise is not supported anymore. Exiting")
        exit()
    elif options.NOISE == "vuvuzela":
        from icecube import vuvuzela

        tray.AddModule(
            "Vuvuzela",
            "vuvuzela_noise",
            InputHitSeriesMapName=mcpe_to_pmt,
            OutputHitSeriesMapName=mcpe_to_pmt + "_withNoise",
            StartWindow=-11 * I3Units.microsecond,
            EndWindow=11 * I3Units.microsecond,
            IceTop=False,
            InIce=True,
            ScaleFactor=1.0,
            DeepCoreScaleFactor=1,
            DOMsToExclude=[],  # This will be cleaned later by DOM launch cleaner
            RandomService="I3RandomService",
            SimulateNewDOMs=True,
            DisableLowDTCutoff=True,
            UseIndividual=True,
        )

        mcpeout = mcpe_to_pmt + "_withNoise"
    elif options.NOISE == "none":
        print("\n*******WARNING: Noiseless simulation!!********\n")
        # exit()
        mcpeout = mcpe_to_pmt
    else:
        print("Pick a valid noise model!")
        exit()

    tray.AddModule(
        "PMTResponseSimulator",
        "rosencrantz",
        Input=mcpeout,
        Output=mcpeout + "_weighted",
        MergeHits=True,
    )

    tray.AddModule(
        "DOMLauncher",
        "guildenstern",
        Input=mcpeout + "_weighted",
        Output="InIceRawData_unclean",
        UseTabulatedPT=True,
    )

    tray.AddModule("I3DOMLaunchCleaning", "launchcleaning")(
        ("InIceInput", "InIceRawData_unclean"),
        ("InIceOutput", "InIceRawData"),
        ("FirstLaunchCleaning", False),
    )

    # Dropping frames without InIceRawData
    tray.AddModule(
        BasicDOMFilter,
        "FilterNullInIce",
        Streams=[icetray.I3Frame.DAQ, icetray.I3Frame.Physics],
    )
    ###### triggering
    tray.AddModule(
        "Delete",
        "delete_triggerHierarchy",
        Keys=["I3TriggerHierarchy", "TimeShift", "CleanIceTopRawData"],
    )

    tray.AddSegment(
        trigger_sim.TriggerSim,
        "trig",
        gcd_file=gcd_file,
        time_shift_args={"SkipKeys": ["BundleGen"]},  #
        # added in run_id
        run_id=1,
    )

    SkipKeys = [
        # "MCPMTResponseMap",
        # "MCTimeIncEventID",
        "I3MCTree_clsim",
        "I3MCTree_preMuonProp",
        # "I3Photons",
        # "clsim_stats",
        # "InIceRawData_unclean",
    ]

    tray.AddModule(
        "I3Writer",
        "writer",
        SkipKeys=SkipKeys,
        Filename=outfile_temp,
        Streams=[
            icetray.I3Frame.DAQ,
            icetray.I3Frame.Physics,
            icetray.I3Frame.TrayInfo,
            icetray.I3Frame.Simulation,
        ],
    )

    if write_hdf5:

        outfile_hdf5_temp = outfile_temp.replace('.i3.zst',".hdf5")

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

        # get weights
        tray.AddModule(
            store_weights,
            "store_weights",
            Streams=[icetray.I3Frame.DAQ],
        )

        tray.AddSegment(
            I3SimHDFWriter,
            output = outfile_hdf5_temp,
            existing_header=True,
            keys = HDF5_KEYS,
        )

    tray.AddModule("TrashCan", "adios")

    tray.Execute()
    tray.Finish()

    if osg == "True":
        copy_to_NPX(outfile)
        copy_to_NPX(outfile.replace('.i3.zst',".hdf5"))
    else:
        os.system(str("mv {} {}".format(outfile_temp, options.OUTFILE)))
        os.system(str("mv {} {}".format(outfile_hdf5_temp, options.OUTFILE.replace('.i3.zst',".hdf5"))))

    t1 = time.clock()

    print("Time it took: {}s".format(t1-t0))
    print('done ...')


if __name__ == '__main__':
    main()

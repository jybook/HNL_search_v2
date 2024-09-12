#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.2.1/icetray-start
#METAPROJECT /data/user/jbook/i3/build/


### The Photon Level simulation processing script ###
import time

import os
import numpy as np

from icecube.simprod import segments

from optparse import OptionParser

import subprocess

from I3Tray import I3Tray
from icecube import icetray, dataclasses, LeptonInjector, dataio, phys_services, clsim  # LI is needed for hdf5 writing

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


def BasicHitFilter(frame, photon_series = "I3Photons"):
    hits = 0
    if frame.Has(photon_series):
        hits = len(frame.Get(photon_series))
    if hits > 0:
        return True
    else:
        return False


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
        default="/data/ana/BSM/HNL/MC/test_files/test_output_Phot.i3.zst",
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
        "-l",
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
        default=1.2,  # Using efficiency > 1 as default so we can support systematics sets
        dest="EFFICIENCY",
        help="DOM Efficiency ... the same as UnshadowedFraction [default: %default]",
    )
    parser.add_option(
        "-m",
        "--icemodel",
        default="spice_3.2.1",
        dest="ICEMODEL",
        help="Ice model (spice_3.2.1, spice_mie, spice_lea, etc) [default: %default]",
    )
    parser.add_option(
        "-a",
        "--holeice",
        default="ANGSENS/angsens/as.flasher_p1_0.30_p2_-1",
        dest="HOLEICE",
        help="Pick the hole ice parameterization, corresponds to a file name path relative to $I3_SRC/ice-models/resources/models/ [default: %default]",
    )
    parser.add_option(
        "-c",
        "--crossenergy",
        type="float",
        default=30.0,
        dest="CROSSENERGY",
        help="The cross energy where the hybrid clsim approach will be used [default: %default]",
    )
    parser.add_option(
        "-t", action="store_true", dest="GPU", default=False, help="Run on GPUs or CPUs [default: %default]"
    )
    parser.add_option(
        "--osg", type="string", default="False", help="Is this job running on the OSG? [default: %default]"
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

    t0 = time.perf_counter() #time.clock()
    print('Starting ...')

    osg = options.osg
    options.FILENR = int(options.FILENR)
    options.RUNNUMBER = int(options.RUNNUMBER)
    write_hdf5  = options.write_hdf5

    if options.GPU:
        CPU = False
    else:
        CPU = True

    print("Using RUNNUMBER: {}".format(options.RUNNUMBER))

    # Now fire up the random number generator with that seed
    # from globals import max_num_files_per_dataset
    max_num_files_per_dataset = 100000
    randomService = phys_services.I3SPRNGRandomService(
        seed=options.RUNNUMBER,
        nstreams=max_num_files_per_dataset,
        streamnum=options.FILENR,
    )

    if osg == "True":

        # gcdfile_NPX = str("gsiftp://gridftp.icecube.wisc.edu" + options.GCDFILE)
        # gcdfile = str(os.getcwd() + "/" + options.GCDFILE.split("/")[-1])

        gcdfile = options.GCDFILE  # if this is in cvmfs then all is good

        infile_NPX = str("gsiftp://gridftp.icecube.wisc.edu" + options.INFILE)
        infile = str(os.getcwd() + "/" + options.INFILE.split("/")[-1])

        # copy_to_OSG(gcdfile_NPX)
        copy_to_OSG(infile_NPX)

        infiles = [gcdfile, infile]
        outfile = str("gsiftp://gridftp.icecube.wisc.edu" + options.OUTFILE)
        outfile_temp = str(os.getcwd() + "/" + outfile.split("/")[-1])

        outfile_hdf5 = outfile.replace('.i3.zst',".hdf5")
        outfile_hdf5_temp = outfile_temp.replace('.i3.zst',".hdf5")
    else:
        gcdfile = options.GCDFILE
        infiles = [gcdfile, options.INFILE]

    # # Set DOM efficiency to x2 the default (to procude upgrade style simulation)
    # options.EFFICIENCY = 2.4

    icemodel_path = os.path.join(os.path.expandvars("$I3_SRC/ice-models/resources/models/ICEMODEL/"), options.ICEMODEL)
    print("Ice model {}".format(icemodel_path))
    print("DOM efficiency: {}".format(options.EFFICIENCY))
    print("Setting cross energy: {}".format(float(options.CROSSENERGY), "GeV"))
    print("Using CPUs {}".format(CPU))
    print("Using GPUs {}".format(options.GPU))

    # now start the actual processing
    tray = I3Tray()
    tray.AddModule("I3Reader", "reader", FilenameList=infiles)
    tray.AddModule("I3GeometryDecomposer", "I3ModuleGeoMap")

#     gcd_file = dataio.I3File(gcdfile)

    tray.AddModule("Rename", keys=["I3MCTree", "I3MCTree_preMuonProp"])
    tray.AddSegment(segments.PropagateMuons, "PropagateMuons", RandomService=randomService)

    print("OpenCL devices: ")
    for d in clsim.I3CLSimOpenCLDevice.GetAllDevices():
        icetray.logging.log_warn(str(d))
    print("---------------------")

    photon_series = "I3Photons"
    tray.AddSegment(
        clsim.I3CLSimMakePhotons,
        "goCLSIM",
        UseCPUs=CPU,
        UseGPUs=options.GPU,
        MCTreeName="I3MCTree",
        OutputMCTreeName="I3MCTree_clsim",
        FlasherInfoVectName=None,
        #MMCTrackListName="MMCTrackList",
        PhotonSeriesName=photon_series,
#         ParallelEvents=1000,
        RandomService=randomService,
        IceModelLocation=icemodel_path,
        UseI3PropagatorService=False,
        UseGeant4=True,
        CrossoverEnergyEM=0.1,
        CrossoverEnergyHadron=float(options.CROSSENERGY),
        StopDetectedPhotons=True,
        HoleIceParameterization=os.path.expandvars(
            "$I3_SRC/ice-models/resources/models/%s" % options.HOLEICE
        ),
        DoNotParallelize=False,
        DOMOversizeFactor=1.0,
        UnshadowedFraction=options.EFFICIENCY,
        GCDFile=gcdfile,
#       ExtraArgumentsToI3CLSimModule={this is the default value "DoublePrecision": False,  # will impact performance if true
#         "StatisticsName": "clsim_stats",
#             "IgnoreDOMIDs": [],
#         },
        ExtraArgumentsToI3PhotonPropagationClientModule={
            "StatisticsName": "clsim_stats",
        },

    )

    # Tested that all frames go through CLSIM. Removing the ones without any hits to save space.
    tray.AddModule(
        BasicHitFilter,
        "FilterNullPhotons",
        photon_series=photon_series,
        Streams=[icetray.I3Frame.DAQ, icetray.I3Frame.Physics],
    )

    SkipKeys = [
        "I3MCTree_bak",
        "MCPESeriesMap",
        "MCPESeriesMapParticleIDMap",
        ]

    if osg == "True":
        this_outfile = outfile_temp
        this_outfile_hdf5 = outfile_hdf5_temp
    else:
        this_outfile = options.OUTFILE
        this_outfile_hdf5 = this_outfile.replace('.i3.zst',".hdf5")

    tray.AddModule(
        "I3Writer",
        "writer",
        SkipKeys=SkipKeys,
        Filename=this_outfile,
        Streams=[
            icetray.I3Frame.DAQ,
            icetray.I3Frame.Physics,
            icetray.I3Frame.TrayInfo,
            icetray.I3Frame.Simulation,
        ],
    )

    if write_hdf5:

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
            output = this_outfile_hdf5,
            keys = HDF5_KEYS,
        )

    tray.AddModule("TrashCan", "adios")

    tray.Execute()
    tray.Finish()

    if osg == "True":
        copy_to_NPX(outfile)
        copy_to_NPX(outfile_hdf5)

    t1 = time.perf_counter() #time.clock()

    print("Time it took: {}s".format(t1-t0))
    print('done ...')


if __name__ == '__main__':
    main()

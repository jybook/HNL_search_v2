#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.2.1/icetray-start
#METAPROJECT /data/user/jbook/i3/build


### The Generation Level simulation processing script for the HNL signal set ###
import time
import os
import numpy as np

from icecube.icetray import I3Tray
from icecube import icetray, dataclasses, LeptonInjector, dataio

from icecube.icetray import I3Units
from icecube.LeptonInjector import weight_hnl_generation, load_generation_weighting_files, weight_hnl_lifetime_framewise
from icecube.LeptonInjector.oscnext_frame_objects_weighting import create_weight_dict
from icecube.hdfwriter import I3SimHDFWriter

from optparse import OptionParser

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


########## End


def main():
    '''
    Main function: Run OptionParser, load needed inputs (for weighting), set RNG and then call LeptonInjector. Skip empty frames and then add weighting and write hdf5.
    '''
    
    from optparse import OptionParser

    usage = "usage: %prog [options]"
    parser = OptionParser(usage)
    parser.add_option("-o", "--outfile", dest="outfile", type="string", help="Outfile name")
    parser.add_option("-i", "--infile", dest="infile", type="string", help="Infile name")
    parser.add_option("-m","--mass", dest="HNL_mass", type="string", help="Mass of the HNL", default="0.3")
    parser.add_option("--hdf5", action="store_false", dest="write_hdf5", default=True, help="Write hdf5 file [default: %default]")

    options,args = parser.parse_args()
    infile 		= options.infile
    outfile 	= options.outfile
    HNL_mass    = float(options.HNL_mass)
    write_hdf5  = options.write_hdf5

    t0 = time.perf_counter() # time.clock()

    print('Starting ...')

    print('Outfile: {}'.format(outfile))
    
    if(write_hdf5):print('Also writing hdf5 file')

    # The .lic file contains the simulation data used for weighting.
    generation_data_file_path = outfile.replace('.i3.zst', '.lic')

    xs_location = os.path.expandvars('$I3_SRC/LeptonInjector/resources/new_cross_sections/M_{:04.0f}MeV/splines'.format(HNL_mass*1e3))
    print('Using xs from {}'.format(xs_location))

    # load flux and cross sections for the weighting
    nusquids_flux, leptonweighter_xsec = load_generation_weighting_files(
        flux_file='{}/LeptonInjector/resources/flux_files/atmospheric_0_0.000000_0.000000_0.000000_0.000000_0.000000_0.000000_2006_2014.hdf5'.format(os.environ['I3_SRC']),
        xs_location=xs_location,
    )

    # now start the actual tray processing
    tray = I3Tray()
    tray.AddModule("I3Reader", "reader", FilenameList=infile)

    ##### add weighting functions #####

    # create weight dict
    tray.AddModule(
        create_weight_dict,
        streams=[icetray.I3Frame.DAQ]
    )

    # weight generation (OneWeight and OneWeight+flux+osc)
    tray.AddModule(
        weight_hnl_generation,
        nusquids_flux=nusquids_flux,
        leptonweighter_xsec=leptonweighter_xsec,
        leptoninjector_config=str(generation_data_file_path),
        streams=[icetray.I3Frame.DAQ],
    )

    # weight reference (lifetime and reference weight (OneWeight+flux+osc+lifetime))
    tray.AddModule(
        weight_hnl_lifetime_framewise,
        U_tau4_sq=1e-03,
        streams=[icetray.I3Frame.DAQ],
    )

    ##### end #####

    tray.Add(
        "I3Writer",
        filename = outfile,
        Streams = [icetray.I3Frame.DAQ,],
    )

    if write_hdf5:

        outfile_hdf5 = outfile.replace('.i3.zst',".hdf5")

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
            output = outfile_hdf5,
            keys = HDF5_KEYS,
        )

    tray.Execute()
    tray.Finish()


    t1 = time.perf_counter() # time.clock()

    print("Time it took: {}s".format(t1-t0))
    print('done ...')


if __name__ == '__main__':


    main()

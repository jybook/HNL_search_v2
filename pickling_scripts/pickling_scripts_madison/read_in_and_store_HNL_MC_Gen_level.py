#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/icetray-start
#METAPROJECT /cvmfs/icecube.opensciencegrid.org/users/lfischer/oscnext_hnl_meta/build

#
# Author: Leander Fischer
#

import time

t0 = time.time()

from ProgressBar import print_progress

import sys, os
import numpy as np

import glob

from icecube import (
    icetray,
    dataclasses,
    dataio,
    phys_services,
    simclasses,
    recclasses,
    LeptonInjector,
)

from icecube.millipede import MillipedeFitParams

import icecube.icetray.i3logging as logging

import pandas as pd

import collections

from optparse import OptionParser

usage = "usage: %prog [options]"
parser = OptionParser(usage)
parser.add_option(
    "-i",
    "--inpath",
    default="/data/ana/BSM/HNL/MC/190607/Gen",
    dest="INPATH",
    help="Read input from INPATH subfolders (.i3{.gz/.zst} format)",
)
parser.add_option(
    "-o",
    "--outfile",
    default="/data/ana/BSM/HNL/MC/pickle_files/gen_level",
    dest="OUTPATH",
    help="Write output pickle files to OUTPATH",
)
parser.add_option(
    "-v",
    "--verbose",
    default=False,
    action="store_true",
    dest="VERBOSE",
    help="logging.log_info info level output",
)
parser.add_option(
    "-r",
    "--real",
    default=False,
    action="store_true",
    dest="REAL",
    help="Real Run. If not set, test is done on first 2 files of each set, whithout storing.",
)
parser.add_option(
    "-l",
    "--logtofile",
    default=False,
    action="store_true",
    dest="LOGTOGILE",
    help="Store log into a file with fit name as name.",
)
(options, args) = parser.parse_args()

if not options.INPATH:
    parser.error("Inpath not specified.")
if not options.OUTPATH:
    parser.error("Outpath not specified.")

if options.VERBOSE:
    logging.set_level("INFO")
else:
    logging.set_level("WARN")


if options.LOGTOGILE:
    logstring = "pickle_log_190607_generator_level"
    logpath = os.path.join(options.OUTPATH, logstring)
    logging.log_info("Logpath: {}".format(logpath))
    logging.rotating_files(logpath)


base_path = options.INPATH
logging.log_info("Inpath: {}".format(base_path))

file_names = []


file_names = np.sort(glob.glob(base_path + "/**/*.i3.zst"))
if not len(file_names):
    file_names = np.sort(glob.glob(base_path + "/**/*.i3.bz2"))

num_files = len(file_names)
logging.log_info("Number of files: {}".format(num_files))
if not options.REAL:
    logging.log_info("Example file(first): {}".format(file_names[0]))


"""
Define lists to store the read in True/Reco data
"""

# summary variables
failed_reco = []
empty_files = []

# # true primary variables
# p_true_pos = [[], [], []]
# p_true_energy = []
# p_true_zenith = []
# p_true_azimuth = []
# p_true_time = []

# # true HNL decay daughter variables
# decay_angles = []
# decay_types = []
# decay_shapes = []
# decay_location = []

# # true first (DIS) cascade variables
# casc_0_true_pos = [[], [], []]
# casc_0_true_energy = []
# casc_0_true_zenith = []
# casc_0_true_azimuth = []
# casc_0_true_time = []

# # true second (HNL decay) cascade variables
# casc_1_true_pos = [[], [], []]
# casc_1_true_energy = []
# casc_1_true_zenith = []
# casc_1_true_azimuth = []
# casc_1_true_time = []

# true HNL variables
HNL_true_pos = [[], [], []]
HNL_true_energy = []
HNL_true_zenith = []
HNL_true_azimuth = []
HNL_true_time = []

# true additional variables
# decay_l_true = []
particle_type = []

# event properties parameters
HNL_mass = []
total_energy = []
lifetime = []
mHNL_min = []
mHNL_max = []
distance = []
distanceMin = []
distanceMax = []

file_count = 0
n_events = 0

for filepath in file_names:
    # print progressbar
    print_progress(file_count, num_files)

    # logging.log_info(filepath)

    file_name = os.path.split(filepath)[-1]
    # logging.log_info(file_name)

    infile = dataio.I3File(filepath)

    count = 0

    while infile.more():  # iterate through all physics frames end save values

        break_bool = False

        try:
            f = infile.pop_daq()
        except RuntimeError as err:
            # logging.log_info('{}'.format(err))
            logging.log_info(
                "Skipping remaining frames of file {}. (No physics frames present anymore.)".format(
                    file_name
                )
            )
            failed_reco.append(file_name)
            continue

        """
        extract true variables
        """

        mctree = f["I3MCTree"]
        p_true = mctree.primaries[0]

        p_daughters = mctree.get_daughters(p_true)

        assert len(p_daughters) == 2

        for p_daughter in p_daughters:
            if p_daughter.type == dataclasses.I3Particle.Hadrons:
                casc_0_true = p_daughter
                # print(casc_0_true)
            else:
                hnl_true = p_daughter

        hnl_daughters = mctree.get_daughters(hnl_true)

        if len(hnl_daughters) == 0:
            logging.log_info(
                "Event has no (visible) decay products. Skipping this event for the moment."
            )
            continue

        # # true HNL decay daughter variables (create temp lists, so the overall lists have the proper shape)
        # decay_angles_temp = []
        # decay_types_temp = []
        # decay_shapes_temp = []
        # decay_location_temp = []

        # assert len(hnl_daughters) > 0
        # for count_hnl_daughters, hnl_daughter in enumerate(hnl_daughters):

        #     # decay_angles_temp.append(
        #     #     phys_services.I3Calculator.angle(hnl_daughter, casc_0_true)
        #     #     / icetray.I3Units.deg
        #     # )
        #     # decay_types_temp.append(hnl_daughter.type.real)
        #     # decay_shapes_temp.append(hnl_daughter.shape.real)
        #     # decay_location_temp.append(hnl_daughter.location_type.real)

        #     if not count_hnl_daughters:
        #         casc_1_true = hnl_daughter
        #     else:
        #         assert casc_1_true.pos == hnl_daughter.pos
        #         if np.isnan(hnl_daughter.energy):
        #             break_bool = True
        #         casc_1_true.energy = casc_1_true.energy + hnl_daughter.energy

        # if break_bool:
        #     # logging.log_info('There is a NaN energy in a hnl daugther particle. Skipping this event.')
        #     continue

        # # true primary variables
        # p_true_pos[0].append(p_true.pos.x)
        # p_true_pos[1].append(p_true.pos.y)
        # p_true_pos[2].append(p_true.pos.z)
        # p_true_energy.append(p_true.energy)
        # p_true_zenith.append(p_true.dir.zenith)
        # p_true_azimuth.append(p_true.dir.azimuth)
        # p_true_time.append(p_true.time)

        # # true HNL decay daughter variables
        # decay_angles.append(decay_angles_temp)
        # decay_types.append(decay_types_temp)
        # decay_shapes.append(decay_shapes_temp)
        # decay_location.append(decay_location_temp)

        # # true first (DIS) cascade variables
        # casc_0_true_pos[0].append(casc_0_true.pos.x)
        # casc_0_true_pos[1].append(casc_0_true.pos.y)
        # casc_0_true_pos[2].append(casc_0_true.pos.z)
        # casc_0_true_energy.append(casc_0_true.energy)
        # casc_0_true_zenith.append(casc_0_true.dir.zenith)
        # casc_0_true_azimuth.append(casc_0_true.dir.azimuth)
        # casc_0_true_time.append(casc_0_true.time)

        # # true second (HNL decay) cascade variables
        # casc_1_true_pos[0].append(casc_1_true.pos.x)
        # casc_1_true_pos[1].append(casc_1_true.pos.y)
        # casc_1_true_pos[2].append(casc_1_true.pos.z)
        # casc_1_true_energy.append(casc_1_true.energy)
        # casc_1_true_zenith.append(casc_1_true.dir.zenith)
        # casc_1_true_azimuth.append(casc_1_true.dir.azimuth)
        # casc_1_true_time.append(casc_1_true.time)

        # true HNL variables
        HNL_true_pos[0].append(hnl_true.pos.x)
        HNL_true_pos[1].append(hnl_true.pos.y)
        HNL_true_pos[2].append(hnl_true.pos.z)
        HNL_true_energy.append(hnl_true.energy)
        HNL_true_zenith.append(hnl_true.dir.zenith)
        HNL_true_azimuth.append(hnl_true.dir.azimuth)
        HNL_true_time.append(hnl_true.time)

        # true additional variables
        # decay_l_true.append(
        #     phys_services.I3Calculator.distance(casc_0_true, casc_1_true)
        #     / icetray.I3Units.m
        # )
        particle_type.append(p_true.pdg_encoding)

        # event properties parameters
        HNL_mass.append(f["EventProperties"].mHNL)
        total_energy.append(f["EventProperties"].totalEnergy)
        lifetime.append(f["EventProperties"].lifetime)
        mHNL_min.append(f["EventProperties"].mHNL_min)
        mHNL_max.append(f["EventProperties"].mHNL_max)
        distance.append(f["EventProperties"].distance)
        distanceMin.append(f["EventProperties"].distanceMin)
        distanceMax.append(f["EventProperties"].distanceMax)

        n_events += 1
        count += 1
        # if(count == 2):break
        # break

    infile.close()
    if not count:
        empty_files.append(file_name)
    file_count += 1

    if not options.REAL:
        if file_count == 1:
            break

    # # test new features.
    # if(file_count == 1): break

logging.log_info("All events: {}".format(n_events))

logging.log_info("Number of empty files: {}".format(len(empty_files)))

# create data dict with variables
data_dict = {
    # # true primary variables
    # "true_x": np.array(p_true_pos[0]),
    # "true_y": np.array(p_true_pos[1]),
    # "true_z": np.array(p_true_pos[2]),
    # "true_energy": np.array(p_true_energy),
    # "true_zenith": np.array(p_true_zenith),
    # "true_azimuth": np.array(p_true_azimuth),
    # "true_time": np.array(p_true_time),
    # # true first (DIS) cascade variables
    # "casc0_true_x": np.array(casc_0_true_pos[0]),
    # "casc0_true_y": np.array(casc_0_true_pos[1]),
    # "casc0_true_z": np.array(casc_0_true_pos[2]),
    # "casc0_true_energy": np.array(casc_0_true_energy),
    # "casc0_true_zenith": np.array(casc_0_true_zenith),
    # "casc0_true_azimuth": np.array(casc_0_true_azimuth),
    # "casc0_true_time": np.array(casc_0_true_time),
    # #  true second (HNL decay) cascade variables
    # "casc1_true_x": np.array(casc_1_true_pos[0]),
    # "casc1_true_y": np.array(casc_1_true_pos[1]),
    # "casc1_true_z": np.array(casc_1_true_pos[2]),
    # "casc1_true_energy": np.array(casc_1_true_energy),
    # "casc1_true_zenith": np.array(casc_1_true_zenith),
    # "casc1_true_azimuth": np.array(casc_1_true_azimuth),
    # "casc1_true_time": np.array(casc_1_true_time),
    # true HNL variables
    "HNL_true_x": np.array(HNL_true_pos[0]),
    "HNL_true_y": np.array(HNL_true_pos[1]),
    "HNL_true_z": np.array(HNL_true_pos[2]),
    "HNL_true_energy": np.array(HNL_true_energy),
    "HNL_true_zenith": np.array(HNL_true_zenith),
    "HNL_true_azimuth": np.array(HNL_true_azimuth),
    "HNL_true_time": np.array(HNL_true_time),
    # # true HNL decay daughter variables
    # "decay_angles": np.array(decay_angles),
    # "decay_types": np.array(decay_types),
    # "decay_shapes": np.array(decay_shapes),
    # "decay_location": np.array(decay_location),
    # true additional variables
    # "true_decayL": np.array(decay_l_true),
    "p_type": np.array(particle_type),
    # event properties parameters
    "hnl_mass": np.array(HNL_mass),
    "total_energy": np.array(total_energy),
    "lifetime": np.array(lifetime),
    "mHNL_min": np.array(mHNL_min),
    "mHNL_max": np.array(mHNL_max),
    "distance": np.array(distance),
    "distanceMin": np.array(distanceMin),
    "distanceMax": np.array(distanceMax),
}
print_progress(file_count, num_files)

data = pd.DataFrame(data_dict)

# if not options.REAL:
#     data.info()

# store extracted data as pickle file
store_name = "190607_generator_level_reconstructed_data.pckl"
logging.log_info("Storename: {}".format(store_name))

filepath = os.path.join(options.OUTPATH, store_name)
logging.log_info("Storepath: {}".format(filepath))

if options.REAL:
    # data.to_pickle(path=filepath)
    # in case this extraction breaks i should try with the following lines
    try:
        data.to_pickle(path=filepath)
    except:
        logging.log_info('Storing as pickle did not work, try hdf5.')
        try:
            data.to_hdf(filepath.replace('pckl', 'hdf5'))
        except:
            logging.log_info('Storing as hdf5 did not work, try json.')
            try:
                data.to_json(filepath.replace('pckl', 'json'))
            except:
                logging.log_info('Nothing works these days.')

t1 = time.time()
logging.log_info("Total time it took: {:.3f} s".format(t1 - t0))

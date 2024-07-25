#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/icetray-start
#METAPROJECT /cvmfs/icecube.opensciencegrid.org/users/lfischer/oscnext_hnl_meta/build

#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# Author: Leander Fischer
#

import time
t0 = time.time()

import sys, os
import numpy as np

import glob

from icecube import icetray, dataclasses, dataio, phys_services, simclasses, recclasses, LeptonInjector

from icecube.millipede import MillipedeFitParams

import icecube.icetray.i3logging as logging

import pandas as pd

import collections

from optparse import OptionParser

usage = "usage: %prog [options]"
parser = OptionParser(usage)
parser.add_option("-o", "--outpath",
                 default="/data/ana/BSM/HNL/MC/pickle_files/all_levels",
                 dest="OUTPATH", help="Write output pickle files to OUTPATH")
parser.add_option("--settype",
                 dest="SETTYPE", help="which type of set to extract")
parser.add_option("--setnumber",
                 dest="SETNUMBER", help="which set to extract")
parser.add_option("--level",
                #  default="Millipede",
                 dest="LEVEL", help="level to extract")
parser.add_option("-v", "--verbose", default=False, action="store_true",
                 dest="VERBOSE", help="logging.log_info info level output")
parser.add_option("-r", "--real", default=False, action="store_true",
                 dest="REAL", help="Real Run. If not set, test is done on first 2 files of each set, whithout storing.")
parser.add_option("-l", "--logtofile", default=False, action="store_true",
                 dest="LOGTOGILE", help="Store log into a file with fit name as name.")
(options,args) = parser.parse_args()

if not options.OUTPATH:
    parser.error('Outpath not specified.')

if options.VERBOSE:
    logging.set_level('INFO')
else:
    logging.set_level('WARN')

level_folder_name = options.LEVEL

if options.LEVEL in ['L3', 'L4', 'L5', 'L6', 'L7']:
    level_folder_name = options.LEVEL.replace('L','level')

if options.LOGTOGILE:
    logstring = 'pickle_log_{}_level_{}'.format(options.SETNUMBER, options.LEVEL)
    logpath = os.path.join(options.OUTPATH, logstring)
    logging.log_info('Logpath: {}'.format(logpath))
    logging.rotating_files(logpath)

if options.LEVEL in ['Taupede', 'Millipede']:
    in_path = '/data/ana/BSM/HNL/MC/oscNext/'
    base_path = os.path.join(os.path.join(os.path.join(in_path, options.SETTYPE), options.SETNUMBER), level_folder_name)
else:
    in_path = '/data/ana/LE/oscNext/pass2/'
    base_path = os.path.join(os.path.join(os.path.join(in_path, options.SETTYPE), level_folder_name), options.SETNUMBER)


logging.log_info('Inpath: {}'.format(base_path))

file_names = []

file_names = np.sort(glob.glob(base_path + '/*.i3.zst'))
if not len(file_names):
    file_names = np.sort(glob.glob(base_path + '/*.i3.bz2'))

num_files = len(file_names)
logging.log_info('Number of files: {}'.format(num_files))
# if not options.REAL:
#     logging.log_info('Example file(first): {}'.format(file_names[0]))

if options.SETTYPE == 'muongun':
    logging.log_info('Using muongun processing chain')
else:
    logging.log_info('Using genie/noise processing chain')

'''
Define lists to store the read in True/Reco data
'''

# summary variables
empty_files = []

# # true primary variables
# p_true_pos = [[], [], []]
# p_true_energy = []
# p_true_zenith = []
# p_true_azimuth = []
# p_true_time = []

# true additional variables
weight = []
int_type = []

# # reco primary variables
# p_reco_pos = [[], [], []]
# p_reco_energy = []
# p_reco_zenith = []
# p_reco_azimuth = []
# p_reco_time = []
# p_reco_length = []

# # reco first (DIS) cascade variables
# casc_0_reco_pos = [[], [], []]
# casc_0_reco_energy = []
# casc_0_reco_zenith = []
# casc_0_reco_azimuth = []
# casc_0_reco_time = []

# # reco second (HNL decay) cascade variables
# casc_1_reco_pos = [[], [], []]
# casc_1_reco_energy = []
# casc_1_reco_zenith = []
# casc_1_reco_azimuth = []
# casc_1_reco_time = []

# # reco additional variables (taupede)
# fitstatus = []
# decay_l_reco = []

oscNext_bool = []
if options.LEVEL in ['L3', 'L4', 'L5', 'L6', 'L7']:
    oscNext_bool_key = '{}_oscNext_bool'.format(options.LEVEL)
else:
    oscNext_bool_key = 'L7_oscNext_bool'

# logl = []
# rlogl = []
# ndof = []

# qtotal = []
# predicted_qtotal = []
# squared_residuals = []

# chi_squared = []
# chi_squared_dof = []

# # reco additional variables (Retro)
# retro_crs_prefit__max_llh = []
# retro_crs_prefit__zero_dllh_mean = []

# # millipede energies
# millipede_total_energy = []
# millipede_casc0_energy_40m = []
# millipede_casc1_energy_40m = []
# millipede_casc0_energy_15m = []
# millipede_casc1_energy_15m = []

file_count = 0
n_events = 0

for filepath in file_names:

    # logging.log_info(filepath)

    file_name = os.path.split(filepath)[-1]
    # logging.log_info(file_name)  
    
    infile = dataio.I3File(filepath)
    
    count = 0

    while(infile.more()):  # iterate through all physics frames end save values
        try:
            f = infile.pop_physics()
        except RuntimeError as err:
            # logging.log_info('{}'.format(err))
            logging.log_info('Skipping remaining frames of file {}. (No physics frames present anymore.)'.format(file_name))
            continue
        
        '''
        extract true variables
        '''

        mctree = f['I3MCTree']
        if(len(mctree.primaries)):
            p_true = mctree.primaries[0]
        else:
            logging.log_info('MC Tree is empty, using default initialised particle.')
            p_true = dataclasses.I3Particle()
            # print(p_true)

        if options.SETTYPE == 'muongun':
            p_true = mctree.get_daughters(p_true)[0]
            try:
                assert(np.abs(p_true.type) == 13)
            except AssertionError:
                logging.log_info('Store particle is not a MuMinus/MuPlus (as it should be)')

        # # true primary variables
        # p_true_pos[0].append(p_true.pos.x)
        # p_true_pos[1].append(p_true.pos.y)
        # p_true_pos[2].append(p_true.pos.z)
        # p_true_energy.append(p_true.energy)
        # p_true_zenith.append(p_true.dir.zenith)
        # p_true_azimuth.append(p_true.dir.azimuth)
        # p_true_time.append(p_true.time)

        # true additional variables
        weight.append(f['I3MCWeightDict']['weight'])
        if f['I3MCWeightDict'].has_key('InteractionType'):
            int_type.append(f['I3MCWeightDict']['InteractionType'])
        else:
            int_type.append(-1)

        # '''
        # extract reco variables
        # '''

        # if f.Has('TaupedeFitManualFitStatus'):
        #     manual_fit_status = f['TaupedeFitManualFitStatus'].value

        # if manual_fit_status in [0,50]:
        #     fitname = 'TaupedeFit'
        #     p_reco = f[fitname]
            
        #     daughters_reco = f[fitname+'Particles']
        #     assert(len(daughters_reco) == 2)
            
        #     casc_0_reco = daughters_reco[0]
        #     casc_1_reco = daughters_reco[1]
            
        #     # extract final LLH for the fit
        #     logl.append(f['TaupedeFitFitParams'].logl)
        #     rlogl.append(f['TaupedeFitFitParams'].rlogl)
        #     ndof.append(f['TaupedeFitFitParams'].ndof)

        #     qtotal.append(f['TaupedeFitFitParams'].qtotal)
        #     predicted_qtotal.append(f['TaupedeFitFitParams'].predicted_qtotal)

        #     squared_residuals.append(f['TaupedeFitFitParams'].squared_residuals)

        #     chi_squared.append(f['TaupedeFitFitParams'].chi_squared)
        #     chi_squared_dof.append(f['TaupedeFitFitParams'].chi_squared_dof)


        # else:
        #     # logging.log_info("Taupede fit of event {0} failed with: {1}\nSetting daughter particles to MC truth.".format(n_events, dataclasses.I3Particle.FitStatus.values[manual_fit_status]))
            
        #     # set llh to -1 and reco particles to true particles if reconstruction did not work
        #     logl.append(-1)
        #     rlogl.append(-1)
        #     ndof.append(-1)

        #     qtotal.append(-1)
        #     predicted_qtotal.append(-1)

        #     squared_residuals.append(-1)

        #     chi_squared.append(-1)
        #     chi_squared_dof.append(-1)
    
        #     p_reco = dataclasses.I3Particle()
            
        #     casc_0_reco = dataclasses.I3Particle()
        #     casc_1_reco = dataclasses.I3Particle()

        # # reco primary variables
        # p_reco_pos[0].append(p_reco.pos.x)
        # p_reco_pos[1].append(p_reco.pos.y)
        # p_reco_pos[2].append(p_reco.pos.z)
        # p_reco_energy.append(p_reco.energy)
        # p_reco_zenith.append(p_reco.dir.zenith)
        # p_reco_azimuth.append(p_reco.dir.azimuth)
        # p_reco_time.append(p_reco.time)
        # p_reco_length.append(p_reco.length)

        # # reco first (DIS) cascade variables    
        # casc_0_reco_pos[0].append(casc_0_reco.pos.x)
        # casc_0_reco_pos[1].append(casc_0_reco.pos.y)
        # casc_0_reco_pos[2].append(casc_0_reco.pos.z)
        # casc_0_reco_energy.append(casc_0_reco.energy)
        # casc_0_reco_zenith.append(casc_0_reco.dir.zenith)
        # casc_0_reco_azimuth.append(casc_0_reco.dir.azimuth)
        # casc_0_reco_time.append(casc_0_reco.time)

        # # reco second (HNL decay) cascade variables
        # casc_1_reco_pos[0].append(casc_1_reco.pos.x)
        # casc_1_reco_pos[1].append(casc_1_reco.pos.y)
        # casc_1_reco_pos[2].append(casc_1_reco.pos.z)
        # casc_1_reco_energy.append(casc_1_reco.energy)
        # casc_1_reco_zenith.append(casc_1_reco.dir.zenith)
        # casc_1_reco_azimuth.append(casc_1_reco.dir.azimuth)
        # casc_1_reco_time.append(casc_1_reco.time)


        # # reco additional variables
        # decay_l_reco.append(phys_services.I3Calculator.distance(casc_0_reco,casc_1_reco)/icetray.I3Units.m)
        # fitstatus.append(manual_fit_status)
        oscNext_bool.append(f[oscNext_bool_key].value)

        # # reco additional variables (Retro)
        # retro_crs_prefit__max_llh.append(f['retro_crs_prefit__max_llh'].value)
        # retro_crs_prefit__zero_dllh_mean.append(f['retro_crs_prefit__zero_dllh']['mean'])

        # # millipede energies
        # if f.Has('MuMillipedeFit'):
        #     trackvector = f['MuMillipedeFit']
        #     taupedeCascade0 = casc_0_reco
        #     taupedeCascade1 = casc_1_reco
        #     E_tot = 0.0
        #     E_cascade0_40m = 0.0
        #     E_cascade1_40m = 0.0
        #     E_cascade0_15m = 0.0
        #     E_cascade1_15m = 0.0
        #     for sec in trackvector:
        #         E_tot += sec.energy
        #         d0 = np.sqrt(np.power(sec.pos.x-taupedeCascade0.pos.x,2)+np.power(sec.pos.y-taupedeCascade0.pos.y,2)+np.power(sec.pos.z-taupedeCascade0.pos.z,2))
        #         d1 = np.sqrt(np.power(sec.pos.x-taupedeCascade1.pos.x,2)+np.power(sec.pos.y-taupedeCascade1.pos.y,2)+np.power(sec.pos.z-taupedeCascade1.pos.z,2))
        #         # 40 m
        #         if d0 < 40 and d1 > 40: E_cascade0_40m += sec.energy
        #         elif d0 > 40 and d1 < 40: E_cascade1_40m += sec.energy
        #         elif d0 < 40 and d1 < 40:
        #             if d0 < d1: E_cascade0_40m += sec.energy
        #             else: E_cascade1_40m += sec.energy
        #         # 15 m
        #         if d0 < 15 and d1 > 15: E_cascade0_15m += sec.energy
        #         elif d0 > 15 and d1 < 15: E_cascade1_15m += sec.energy
        #         elif d0 < 15 and d1 < 15:
        #             if d0 < d1: E_cascade0_15m += sec.energy
        #             else: E_cascade1_15m += sec.energy
        # else:
        #     E_tot = -1.0
        #     E_cascade0_40m = -1.0
        #     E_cascade1_40m = -1.0
        #     E_cascade0_15m = -1.0
        #     E_cascade1_15m = -1.0

        # millipede_total_energy.append(E_tot)
        # millipede_casc0_energy_40m.append(E_cascade0_40m)
        # millipede_casc1_energy_40m.append(E_cascade1_40m)
        # millipede_casc0_energy_15m.append(E_cascade0_15m)
        # millipede_casc1_energy_15m.append(E_cascade1_15m)

        n_events += 1
        count += 1
        # if(count == 2):break
        # break

    infile.close()
    if not count:
        empty_files.append(file_name)
    file_count +=1

    if not options.REAL:
        if(file_count == 1): break
    
    # # test new features.
    # if(file_count == 1): break

logging.log_info('All events: {}'.format(n_events))

logging.log_info('Number of empty files: {}'.format(len(empty_files)))

# logging.log_info('Count type of fitstatusssss: {}'.format(collections.Counter(fitstatus)))

logging.log_info('Total number of events in simulation set: {}'.format(len(weight)))
# logging.log_info('Total rate of events: {:.2f} mHz'.format(1e03 * np.sum(weight)))


# create data dict with variables
data_dict = {
    # # true primary variables
    # 'true_x':np.array(p_true_pos[0]),
    # 'true_y':np.array(p_true_pos[1]),
    # 'true_z':np.array(p_true_pos[2]),
    # 'true_energy':np.array(p_true_energy),
    # 'true_zenith':np.array(p_true_zenith),
    # 'true_azimuth':np.array(p_true_azimuth),
    # 'true_time':np.array(p_true_time),

    # true additional variables
    'weight':np.array(weight),
    'int_type':np.array(int_type),

    # # reco primary variables
    # 'reco_x':np.array(p_reco_pos[0]),
    # 'reco_y':np.array(p_reco_pos[1]),
    # 'reco_z':np.array(p_reco_pos[2]),
    # 'reco_energy':np.array(p_reco_energy),
    # 'reco_zenith':np.array(p_reco_zenith),
    # 'reco_azimuth':np.array(p_reco_azimuth),
    # 'reco_time':np.array(p_reco_time),
    # 'reco_length':np.array(p_reco_length),

    # # reco first (DIS) cascade variables
    # 'casc0_reco_x':np.array(casc_0_reco_pos[0]),
    # 'casc0_reco_y':np.array(casc_0_reco_pos[1]),
    # 'casc0_reco_z':np.array(casc_0_reco_pos[2]),
    # 'casc0_reco_energy':np.array(casc_0_reco_energy),
    # 'casc0_reco_zenith':np.array(casc_0_reco_zenith),
    # 'casc0_reco_azimuth':np.array(casc_0_reco_azimuth),
    # 'casc0_reco_time':np.array(casc_0_reco_time),

    # # reco second (HNL decay) cascade variables
    # 'casc1_reco_x':np.array(casc_1_reco_pos[0]),
    # 'casc1_reco_y':np.array(casc_1_reco_pos[1]),
    # 'casc1_reco_z':np.array(casc_1_reco_pos[2]),
    # 'casc1_reco_energy':np.array(casc_1_reco_energy),
    # 'casc1_reco_zenith':np.array(casc_1_reco_zenith),
    # 'casc1_reco_azimuth':np.array(casc_1_reco_azimuth),
    # 'casc1_reco_time':np.array(casc_1_reco_time),

    # # reco additional variables
    # 'fitstatus':np.array(fitstatus),
    # 'reco_decayL':np.array(decay_l_reco),
    'oscNext_bool':np.array(oscNext_bool),
    
    # 'logl':np.array(logl),
    # 'rlogl':np.array(rlogl),
    # 'ndof':np.array(ndof),

    # 'qtotal':np.array(qtotal),
    # 'predicted_qtotal':np.array(predicted_qtotal),

    # 'squared_residuals':np.array(squared_residuals),

    # 'chi_squared':np.array(chi_squared),
    # 'chi_squared_dof':np.array(chi_squared_dof),

    # # reco additional variables (Retro)
    # 'retro_crs_prefit__max_llh':np.array(retro_crs_prefit__max_llh),
    # 'retro_crs_prefit__zero_dllh_mean':np.array(retro_crs_prefit__zero_dllh_mean),

    # # millipede energies
    # 'millipede_total_energy':np.array(millipede_total_energy),
    # 'millipede_casc0_energy_40m':np.array(millipede_casc0_energy_40m),
    # 'millipede_casc1_energy_40m':np.array(millipede_casc1_energy_40m),
    # 'millipede_casc0_energy_15m':np.array(millipede_casc0_energy_15m),
    # 'millipede_casc1_energy_15m':np.array(millipede_casc1_energy_15m),
    }

data = pd.DataFrame(data_dict)

# if not options.REAL:
#     data.info()

# store extracted data as pickle file
store_name = 'data_{}_level_{}.pckl'.format(options.SETNUMBER, options.LEVEL)
logging.log_info('Storename: {}'.format(store_name))

filepath = os.path.join(options.OUTPATH, store_name)
logging.log_info('Storepath: {}'.format(filepath))

# if options.REAL:
#     data.to_pickle(path=filepath)

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

# data.to_pickle(path=filepath)
# data.to_hdf(filepath.replace('pckl', 'hdf5'), key='df')
# data.to_json(filepath.replace('pckl', 'json'))
        
t1 = time.time()
logging.log_info('Total time it took: {:.3f} s'.format(t1-t0))

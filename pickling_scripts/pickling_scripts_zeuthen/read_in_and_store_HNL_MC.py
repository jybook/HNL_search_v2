#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Use with OscNext meta-project (e.g. Leander's):
# unset OS_ARCH; eval `/cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/setup.sh
# export PYTHONPATH=$PYTHONPATH:/afs/ifh.de/group/amanda/scratch/lfischer/software/fridge
# /afs/ifh.de/user/l/lfischer/scratch/ICsoft/meta-projects/oscnext_meta/build/env-shell.sh

#
# Author: Leander Fischer
#

import time
t0 = time.time()


import sys, os
import numpy as np

from icecube import icetray, dataclasses, dataio, phys_services, simclasses, recclasses, LeptonInjector

from icecube.millipede import MillipedeFitParams

import icecube.icetray.i3logging as logging

import pandas as pd

import collections

from optparse import OptionParser

usage = "usage: %prog [options]"
parser = OptionParser(usage)
parser.add_option("-i", "--inpath",
                #  default="/lustre/fs22/group/icecube/lfischer/data/HNL/190605",
                #  default="/lustre/fs22/group/icecube/lfischer/data/HNL/190606",
                 default="/afs/ifh.de/user/l/lfischer/scratch/data/190605",
                 dest="INPATH", help="Read input from INPATH (.i3{.gz/.zst} format)")
parser.add_option("-o", "--outfile",
                 default="/afs/ifh.de/user/l/lfischer/scratch/analysis_phd/plots_all/2021/01_taupede_signal",
                #  default="/afs/ifh.de/user/l/lfischer/scratch/analysis_phd/plots_all/2021/03_preliminary_taupede_190606",
                 dest="OUTPATH", help="Write output pickle files to OUTPATH")
parser.add_option("-v", "--verbose", default=False, action="store_true",
                 dest="VERBOSE", help="logging.log_info info level output")
parser.add_option("-r", "--real", default=False, action="store_true",
                 dest="REAL", help="Real Run. If not set, test is done on first 2 files of each set, whithout storing.")
parser.add_option("-l", "--logtofile", default=False, action="store_true",
                 dest="LOGTOGILE", help="Store log into a file with fit name as name.")
(options,args) = parser.parse_args()

if not options.INPATH:
    parser.error('Inpath not specified.')
if not options.OUTPATH:
    parser.error('Outpath not specified.')

if options.VERBOSE:
    logging.set_level('INFO')
else:
    logging.set_level('WARN')

# logstring = options.INPATH.split('/')[-1] + '_extract_reconstructed_data_log'
# logstring = options.INPATH.split('/')[-1] + '_minuit2_versions_extract_reconstructed_data_log'
# logstring = options.INPATH.split('/')[-1] + '_length_seeds_versions_extract_reconstructed_data_log'
# logstring = options.INPATH.split('/')[-1] + '_clean_build_versions_extract_reconstructed_data_log'
# logstring = options.INPATH.split('/')[-1] + '_effdist_splitinice_versions_extract_reconstructed_data_log'
# logstring = options.INPATH.split('/')[-1] + '_effdist_splitinice_clean_build_versions_extract_reconstructed_data_log'
# logstring = options.INPATH.split('/')[-1] + '_effdist_extract_reconstructed_data_log'
# logstring = options.INPATH.split('/')[-1] + '_effdist_splitinice_pb1_extract_reconstructed_data_log'
# logstring = options.INPATH.split('/')[-1] + '_all_eleven_extract_reconstructed_data_log'

logstring = options.INPATH.split('/')[-1] + '_taupede_reconstructed_data_log'
logging.log_info('Logpath: {}'.format(logstring))

if options.LOGTOGILE:
    logpath = os.path.join(options.OUTPATH, logstring)
    logging.rotating_files(logpath)

taupede_versions = [
    # 'full_fit_v0',
    # 'split_fit_v2',
    # 'split_fit_v3',
    # 'split_fit_v4',

    # 'split_fit_v2_minuit2_simplex',
    # 'split_fit_v2_minuit2_migrad',

    # 'split_fit_v2_length_seeds',

    # 'split_fit_v2_clean_build',
    # 'split_fit_v2_effdist_splitinice',

    # 'split_fit_v2_effdist',

    # 'split_fit_v2_effdist_splitinice_pb1',

    'taupede',
]

for taupede_version in taupede_versions:

    t2 = time.time()

    base_path = os.path.join( os.path.join(options.INPATH, taupede_version) , 'files')
    logging.log_info('######################## NEW MC SET ########################'.format(base_path))
    logging.log_info('Inpath: {}'.format(base_path))

    file_names = []

    for f in os.listdir(base_path):
        if f.endswith('.i3.zst'):
            file_names.append(os.path.join(base_path,f))
    file_names.sort()
    logging.log_info('Number of files: {}'.format(len(file_names)))
    # logging.log_info('Example file(first): {}'.format(file_names[0]))


    '''
    Define lists to store the read in True/Reco data
    '''

    # summary variables
    failed_reco = []
    empty_files = []

    # true primary variables
    p_true_pos = [[], [], []]
    p_true_energy = []
    p_true_zenith = []
    p_true_azimuth = []
    p_true_time = []

    # # true HNL decay daughter variables
    # decay_angles = []
    # decay_types = []
    # decay_shapes = []
    # decay_location = []

    # true first (DIS) cascade variables
    casc_0_true_pos = [[], [], []]
    casc_0_true_energy = []
    casc_0_true_zenith = []
    casc_0_true_azimuth = []
    casc_0_true_time = []

    # true second (HNL decay) cascade variables
    casc_1_true_pos = [[], [], []]
    casc_1_true_energy = []
    casc_1_true_zenith = []
    casc_1_true_azimuth = []
    casc_1_true_time = []

    # true HNL variables
    HNL_true_pos = [[], [], []]
    HNL_true_energy = []
    HNL_true_zenith = []
    HNL_true_azimuth = []
    HNL_true_time = []

    # true additional variables
    weight = []
    decay_l_true = []
    particle_type = []

    # event properties parameters
    HNL_mass = []
    lifetime = []
    mHNL_max = []
    mHNL_min = []
    distance = []
    distanceMin = []
    distanceMax = []

    # reco primary variables
    p_reco_pos = [[], [], []]
    p_reco_energy = []
    p_reco_zenith = []
    p_reco_azimuth = []
    p_reco_time = []

    # reco first (DIS) cascade variables
    casc_0_reco_pos = [[], [], []]
    casc_0_reco_energy = []
    casc_0_reco_zenith = []
    casc_0_reco_azimuth = []
    casc_0_reco_time = []

    # reco second (HNL decay) cascade variables
    casc_1_reco_pos = [[], [], []]
    casc_1_reco_energy = []
    casc_1_reco_zenith = []
    casc_1_reco_azimuth = []
    casc_1_reco_time = []

    # reco additional variables (taupede)
    fitstatus = []
    decay_l_reco = []

    logl = []
    rlogl = []
    ndof = []

    qtotal = []
    predicted_qtotal = []
    squared_residuals = []

    chi_squared = []
    chi_squared_dof = []

    # reco additional variables (Retro)
    retro_crs_prefit__max_llh = []
    retro_crs_prefit__zero_dllh_mean = []

    L7_oscNext_bool = []


    file_count = 0
    n_events = 0

    for filepath in file_names:
        # logging.log_info(filepath)

        file_name = os.path.split(filepath)[-1]
        # logging.log_info(file_name)  
        
        infile = dataio.I3File(filepath)
        
        count = 0

        while(infile.more()):         # iterate through all frames end save values
            try:
                f = infile.pop_physics()
            except RuntimeError as err:
                # logging.log_info('{}'.format(err))
                logging.log_info('Skipping remaining frames of file {}. (No physics frames present anymore.)'.format(file_name))
                failed_reco.append(file_name)
                continue
            
            '''
            extract true variables
            '''

            mctree = f['I3MCTree']        
            p_true = mctree.primaries[0]

            # true primary variables
            p_true_pos[0].append(p_true.pos.x)
            p_true_pos[1].append(p_true.pos.y)
            p_true_pos[2].append(p_true.pos.z)
            p_true_energy.append(p_true.energy)
            p_true_zenith.append(p_true.dir.zenith)
            p_true_azimuth.append(p_true.dir.azimuth)
            p_true_time.append(p_true.time)


            p_daughters = mctree.get_daughters(p_true)

            assert(len(p_daughters) == 2)

            for p_daughter in p_daughters:
                if p_daughter.type == dataclasses.I3Particle.Hadrons:
                    casc_0_true = p_daughter
                    # print(casc_0_true)
            # MODIFICATION (for 190605/190606 set, which has two hadron objects and intermediate HNL particle)
                else:
                    HNL_true = p_daughter
                    casc_1_true = mctree.get_daughters(p_daughter)[0]
                    # print(casc_1_true)
                    
            # # OG Version - START
            #     else:
            #         hnl_true = p_daughter
                    
            # hnl_daughters = mctree.get_daughters(hnl_true)
            
            # if(len(hnl_daughters) == 0):
            #     logging.log_info('Event has no decay products. Skipping this event for the moment.')
            #     continue
            
            # # # true HNL decay daughter variables (create temp lists, so the overall lists have the proper shape)
            # # decay_angles_temp = []
            # # decay_types_temp = []
            # # decay_shapes_temp = []
            # # decay_location_temp = []

            # assert(len(hnl_daughters) > 0)    
            # for count_hnl_daughters, hnl_daughter in enumerate(hnl_daughters):

            #     # decay_angles_temp.append( phys_services.I3Calculator.angle(hnl_daughter,casc_0_true) / icetray.I3Units.deg)
            #     # decay_types_temp.append(hnl_daughter.type)
            #     # decay_shapes_temp.append(hnl_daughter.shape)
            #     # decay_location_temp.append(hnl_daughter.location_type)

            #     if not count_hnl_daughters:
            #         casc_1_true = hnl_daughter
            #     else:
            #         assert(casc_1_true.pos == hnl_daughter.pos)
            #         casc_1_true.energy = casc_1_true.energy + hnl_daughter.energy
            # END
            
            # # true HNL decay daughter variables
            # decay_angles.append(decay_angles_temp)
            # decay_types.append(decay_types_temp)
            # decay_shapes.append(decay_shapes_temp)
            # decay_location.append(decay_location_temp)

            # true first (DIS) cascade variables            
            casc_0_true_pos[0].append(casc_0_true.pos.x)
            casc_0_true_pos[1].append(casc_0_true.pos.y)
            casc_0_true_pos[2].append(casc_0_true.pos.z)
            casc_0_true_energy.append(casc_0_true.energy)
            casc_0_true_zenith.append(casc_0_true.dir.zenith)
            casc_0_true_azimuth.append(casc_0_true.dir.azimuth)
            casc_0_true_time.append(casc_0_true.time)
            
            #  true second (HNL decay) cascade variables
            casc_1_true_pos[0].append(casc_1_true.pos.x)
            casc_1_true_pos[1].append(casc_1_true.pos.y)
            casc_1_true_pos[2].append(casc_1_true.pos.z)
            casc_1_true_energy.append(casc_1_true.energy)
            casc_1_true_zenith.append(casc_1_true.dir.zenith)
            casc_1_true_azimuth.append(casc_1_true.dir.azimuth)
            casc_1_true_time.append(casc_1_true.time)

            # true HNL variables
            HNL_true_pos[0].append(HNL_true.pos.x)
            HNL_true_pos[1].append(HNL_true.pos.y)
            HNL_true_pos[2].append(HNL_true.pos.z)
            HNL_true_energy.append(HNL_true.energy)
            HNL_true_zenith.append(HNL_true.dir.zenith)
            HNL_true_azimuth.append(HNL_true.dir.azimuth)
            HNL_true_time.append(HNL_true.time)

            # true additional variables
            weight.append(f['I3MCWeightDict']['weight'])
            decay_l_true.append(phys_services.I3Calculator.distance(casc_0_true,casc_1_true)/icetray.I3Units.m)
            particle_type.append(p_true.pdg_encoding)

            HNL_mass.append( f['EventProperties'].mHNL )
            lifetime.append( f['EventProperties'].lifetime )
            mHNL_min.append( f['EventProperties'].mHNL_min )
            mHNL_max.append( f['EventProperties'].mHNL_max )
            distance.append( f['EventProperties'].distance )
            distanceMin.append( f['EventProperties'].distanceMin )
            distanceMax.append( f['EventProperties'].distanceMax )

            '''
            extract reco variables
            '''

            if f.Has('TaupedeFitManualFitStatus'):
                manual_fit_status = f['TaupedeFitManualFitStatus'].value
            elif f.Has('RetroSeedParticleManualFitStatus'):
                manual_fit_status = f['RetroSeedParticleManualFitStatus'].value


            if manual_fit_status in [0,50]:
                fitname = 'TaupedeFit'
                p_reco = f[fitname]
                
                daughters_reco = f[fitname+'Particles']
                assert(len(daughters_reco) == 2)
                
                casc_0_reco = daughters_reco[0]
                casc_1_reco = daughters_reco[1]
                
                # extract final LLH for the fit
                logl.append(f['TaupedeFitFitParams'].logl)
                rlogl.append(f['TaupedeFitFitParams'].rlogl)
                ndof.append(f['TaupedeFitFitParams'].ndof)

                qtotal.append(f['TaupedeFitFitParams'].qtotal)
                predicted_qtotal.append(f['TaupedeFitFitParams'].predicted_qtotal)

                squared_residuals.append(f['TaupedeFitFitParams'].squared_residuals)

                chi_squared.append(f['TaupedeFitFitParams'].chi_squared)
                chi_squared_dof.append(f['TaupedeFitFitParams'].chi_squared_dof)


            else:
                # logging.log_info("Taupede fit of event {0} failed with: {1}\nSetting daughter particles to MC truth.".format(n_events, dataclasses.I3Particle.FitStatus.values[manual_fit_status]))
                
                # set llh to -1 and reco particles to true particles if reconstruction did not work
                logl.append(-1)
                rlogl.append(-1)
                ndof.append(-1)

                qtotal.append(-1)
                predicted_qtotal.append(-1)

                squared_residuals.append(-1)

                chi_squared.append(-1)
                chi_squared_dof.append(-1)

                
                p_reco = dataclasses.I3Particle()
                
                casc_0_reco = dataclasses.I3Particle()
                casc_1_reco = dataclasses.I3Particle()

            # reco primary variables
            p_reco_pos[0].append(p_reco.pos.x)
            p_reco_pos[1].append(p_reco.pos.y)
            p_reco_pos[2].append(p_reco.pos.z)
            p_reco_energy.append(p_reco.energy)
            p_reco_zenith.append(p_reco.dir.zenith)
            p_reco_azimuth.append(p_reco.dir.azimuth)
            p_reco_time.append(p_reco.time)

            # reco first (DIS) cascade variables    
            casc_0_reco_pos[0].append(casc_0_reco.pos.x)
            casc_0_reco_pos[1].append(casc_0_reco.pos.y)
            casc_0_reco_pos[2].append(casc_0_reco.pos.z)
            casc_0_reco_energy.append(casc_0_reco.energy)
            casc_0_reco_zenith.append(casc_0_reco.dir.zenith)
            casc_0_reco_azimuth.append(casc_0_reco.dir.azimuth)
            casc_0_reco_time.append(casc_0_reco.time)

            # reco second (HNL decay) cascade variables
            casc_1_reco_pos[0].append(casc_1_reco.pos.x)
            casc_1_reco_pos[1].append(casc_1_reco.pos.y)
            casc_1_reco_pos[2].append(casc_1_reco.pos.z)
            casc_1_reco_energy.append(casc_1_reco.energy)
            casc_1_reco_zenith.append(casc_1_reco.dir.zenith)
            casc_1_reco_azimuth.append(casc_1_reco.dir.azimuth)
            casc_1_reco_time.append(casc_1_reco.time)

            # reco additional variables
            decay_l_reco.append(phys_services.I3Calculator.distance(casc_0_reco,casc_1_reco)/icetray.I3Units.m)
            fitstatus.append(manual_fit_status)
            L7_oscNext_bool.append(f['L7_oscNext_bool'].value)

            # reco additional variables (Retro)
            retro_crs_prefit__max_llh.append(f['retro_crs_prefit__max_llh'].value)
            retro_crs_prefit__zero_dllh_mean.append(f['retro_crs_prefit__zero_dllh']['mean'])

            n_events += 1
            count += 1
            # if(count == 2):break
            # break

        infile.close()
        if not count:
            empty_files.append(file_name)
        file_count +=1

        if not options.REAL:
            if(file_count == 2): break
        
        # # test new features.
        # if(file_count == 1): break

    logging.log_info('All events: {}'.format(n_events))

    logging.log_info('Number of empty files: {}'.format(len(empty_files)))

    logging.log_info('Count type of fitstatusssss: {}'.format(collections.Counter(fitstatus)))

    logging.log_info('Total number of events in simulation set: {}'.format(len(weight)))
    # logging.log_info('Total rate of events: {:.2f} mHz'.format(1e03 * np.sum(weight)))


    # create data dict with variables
    data_dict = {
        # true primary variables
        'true_x':np.array(p_true_pos[0]),
        'true_y':np.array(p_true_pos[1]),
        'true_z':np.array(p_true_pos[2]),
        'true_energy':np.array(p_true_energy),
        'true_zenith':np.array(p_true_zenith),
        'true_azimuth':np.array(p_true_azimuth),
        'true_time':np.array(p_true_time),
  
        # # true first (DIS) cascade variables  
        'casc0_true_x':np.array(casc_0_true_pos[0]),
        'casc0_true_y':np.array(casc_0_true_pos[1]),
        'casc0_true_z':np.array(casc_0_true_pos[2]),
        'casc0_true_energy':np.array(casc_0_true_energy),
        'casc0_true_zenith':np.array(casc_0_true_zenith),
        'casc0_true_azimuth':np.array(casc_0_true_azimuth),
        'casc0_true_time':np.array(casc_0_true_time),
  
        # #  true second (HNL decay) cascade variables
        'casc1_true_x':np.array(casc_1_true_pos[0]),
        'casc1_true_y':np.array(casc_1_true_pos[1]),
        'casc1_true_z':np.array(casc_1_true_pos[2]),
        'casc1_true_energy':np.array(casc_1_true_energy),
        'casc1_true_zenith':np.array(casc_1_true_zenith),
        'casc1_true_azimuth':np.array(casc_1_true_azimuth),
        'casc1_true_time':np.array(casc_1_true_time),

        # true HNL variables
        'HNL_true_x':np.array(HNL_true_pos[0]),
        'HNL_true_y':np.array(HNL_true_pos[1]),
        'HNL_true_z':np.array(HNL_true_pos[2]),
        'HNL_true_energy':np.array(HNL_true_energy),
        'HNL_true_zenith':np.array(HNL_true_zenith),
        'HNL_true_azimuth':np.array(HNL_true_azimuth),
        'HNL_true_time':np.array(HNL_true_time),

        # # # true HNL decay daughter variables
        # 'decay_angles':np.array(decay_angles),
        # 'decay_types':np.array(decay_types),
        # 'decay_shapes':np.array(decay_shapes),
        # 'decay_location':np.array(decay_location),

        # true additional variables
        'weight':np.array(weight),
        'true_decayL':np.array(decay_l_true),
        'p_type':np.array(particle_type),

        # event properties parameters
        'hnl_mass':np.array(HNL_mass), 
        'lifetime':np.array(lifetime),
        'mHNL_min':np.array(mHNL_min),
        'mHNL_max':np.array(mHNL_max),
        'distance':np.array(distance),
        'distanceMin':np.array(distanceMin),
        'distanceMax':np.array(distanceMax),

        # reco primary variables
        'reco_x':np.array(p_reco_pos[0]),
        'reco_y':np.array(p_reco_pos[1]),
        'reco_z':np.array(p_reco_pos[2]),
        'reco_energy':np.array(p_reco_energy),
        'reco_zenith':np.array(p_reco_zenith),
        'reco_azimuth':np.array(p_reco_azimuth),
        'reco_time':np.array(p_reco_time),
  
        # reco first (DIS) cascade variables
        'casc0_reco_x':np.array(casc_0_reco_pos[0]),
        'casc0_reco_y':np.array(casc_0_reco_pos[1]),
        'casc0_reco_z':np.array(casc_0_reco_pos[2]),
        'casc0_reco_energy':np.array(casc_0_reco_energy),
        'casc0_reco_zenith':np.array(casc_0_reco_zenith),
        'casc0_reco_azimuth':np.array(casc_0_reco_azimuth),
        'casc0_reco_time':np.array(casc_0_reco_time),

        # reco second (HNL decay) cascade variables
        'casc1_reco_x':np.array(casc_1_reco_pos[0]),
        'casc1_reco_y':np.array(casc_1_reco_pos[1]),
        'casc1_reco_z':np.array(casc_1_reco_pos[2]),
        'casc1_reco_energy':np.array(casc_1_reco_energy),
        'casc1_reco_zenith':np.array(casc_1_reco_zenith),
        'casc1_reco_azimuth':np.array(casc_1_reco_azimuth),
        'casc1_reco_time':np.array(casc_1_reco_time),
  
        # reco additional variables
        'fitstatus':np.array(fitstatus),
        'reco_decayL':np.array(decay_l_reco),
        
        'logl':np.array(logl),
        'rlogl':np.array(rlogl),
        'ndof':np.array(ndof),

        'qtotal':np.array(qtotal),
        'predicted_qtotal':np.array(predicted_qtotal),

        'squared_residuals':np.array(squared_residuals),

        'chi_squared':np.array(chi_squared),
        'chi_squared_dof':np.array(chi_squared_dof),

        'L7_oscNext_bool':np.array(L7_oscNext_bool),

        # reco additional variables (Retro)
        'retro_crs_prefit__max_llh':np.array(retro_crs_prefit__max_llh),
        'retro_crs_prefit__zero_dllh_mean':np.array(retro_crs_prefit__zero_dllh_mean),
        }

    data = pd.DataFrame(data_dict)
    # if not options.REAL:
    #     data.info()

    # store extracted data as pickle file
    store_name = options.INPATH.split('/')[-1] + '_' + taupede_version + '_reconstructed_data.pckl'
    logging.log_info('Storename: {}'.format(store_name))

    filepath = os.path.join(options.OUTPATH, store_name)
    logging.log_info('Storepath: {}'.format(filepath))

    if options.REAL:
        data.to_pickle(path=filepath)

    t3 = time.time()
    logging.log_info('Time it took: {:.3f} s'.format(t3-t2))

    # # test (only first taupede version)
    # break

t1 = time.time()
logging.log_info('Total time it took: {:.3f} s'.format(t1-t0))

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

import os
import numpy as np

import icecube.icetray.i3logging as logging

import pandas as pd

from optparse import OptionParser

usage = "usage: %prog [options]"
parser = OptionParser(usage)
parser.add_option("-i", "--inpath",
                #  default="/lustre/fs22/group/icecube/lfischer/data/HNL/190605",
                #  default="/lustre/fs22/group/icecube/lfischer/data/HNL/190606",
                 default="/afs/ifh.de/user/l/lfischer/scratch/data/190605",
                 dest="INPATH", help="Read input from INPATH")
parser.add_option("-o", "--outfile",
                 default="/afs/ifh.de/user/l/lfischer/scratch/analysis_phd/plots_all/2021/01_taupede_signal",
                #  default="/afs/ifh.de/user/l/lfischer/scratch/analysis_phd/plots_all/2021/03_preliminary_taupede_190606",
                 dest="OUTPATH", help="Write output pickle files to OUTPATH")
parser.add_option("-v", "--verbose", default=False, action="store_true",
                 dest="VERBOSE", help="logging.log_info info level output")
parser.add_option("-r", "--real", default=False, action="store_true",
                 dest="REAL", help="Real Run. If not set, test is done on first file.")
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

# logstring = options.INPATH.split('/')[-1] + '_extract_runtime_data_log'
# logstring = options.INPATH.split('/')[-1] + '_minuit2_versions_extract_runtime_data_log'
# logstring = options.INPATH.split('/')[-1] + '_lenght_seeds_versions_extract_runtime_data_log'
# logstring = options.INPATH.split('/')[-1] + '_clean_build_versions_extract_runtime_data_log'
# logstring = options.INPATH.split('/')[-1] + '_effdist_splitinice_versions_extract_runtime_data_log'
# logstring = options.INPATH.split('/')[-1] + '_effdist_splitinice_clean_build_versions_extract_runtime_data_log'
# logstring = options.INPATH.split('/')[-1] + '_effdist_extract_runtime_data_log'
# logstring = options.INPATH.split('/')[-1] + '_effdist_splitinice_pb1_extract_runtime_data_log'

logstring = options.INPATH.split('/')[-1] + '_taupede_runtime_data_log'
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

    mc_folder = os.path.join(taupede_version, 'LOGS')
    file_base = os.path.join(options.INPATH, mc_folder)
    logging.log_info('File base path: {}'.format(file_base))

    store_path = options.OUTPATH

    log_files = []

    for f in os.listdir(file_base):
        # print(f)
        if os.path.isfile(os.path.join(file_base, f)):
            log_files.append(os.path.join(file_base,f))

    log_files.sort()
    logging.log_info('Number of files to read in: {}'.format(len(log_files)))
    # logging.log_info('Example file(first): {}'.format(log_files[0]))


    runtime_per_file = []
    files = []

    failed_files = []

    for file_path in log_files:
        # logging.log_info(file_path)

        output_file_bool = False

        inf = open(file_path, 'r')
        previous_line = ''
        for line in inf:
            if output_file_bool:
                outfile = line.split('\n')[0].split('/')[-1]
                # logging.log_info(outfile)
                output_file_bool = False
            if line == 'output file:\n':
                output_file_bool = True

            # if(line.startswith(' *** Break *** ')):
            #     break
                
            if(line.startswith('RuntimeError: problems opening file')):
                # logging.log_info('File {} failed.'.format(outfile))
                failed_files.append(outfile)
                break
            
            if(line.startswith('It took')):
                runtime_per_file.append(int(line.split(' ')[2]))
                files.append(outfile)

        if not options.REAL:
            break
        # break

    logging.log_info('{}'.format(len(runtime_per_file)))
    logging.log_info('{}'.format(len(files)))

    if not options.REAL:
        logging.log_info('{}'.format((runtime_per_file)))
        logging.log_info('{}'.format((files)))

    data_dict = {
        'names':np.array(files),
        'runtime':np.array(runtime_per_file),
    }

    data = pd.DataFrame(data_dict)
    data.sort_values('names', inplace=True)
    # logging.log_info('Dataframe head: {}'.format(data.to_string()))
    # data.info()

    # store extracted data as pickle file
    store_name = options.INPATH.split('/')[-1] + '_' + taupede_version + '_runtime_data.pckl'
    logging.log_info('{}'.format(store_name))

    filepath = os.path.join(options.OUTPATH, store_name)
    logging.log_info('{}'.format(filepath))

    if options.REAL:
        data.to_pickle(path=filepath)

    t3 = time.time()
    logging.log_info('Time it took: {:.3f} s'.format(t3-t2))

    # # for testing purposes
    # break

t1 = time.time()
logging.log_info('Total time it took: {:.3f} s'.format(t1-t0))

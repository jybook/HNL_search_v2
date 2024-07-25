#!/usr/bin/env python
# -*- coding: utf-8 -*- 

from optparse import OptionParser

usage = "usage: %prog [options]"
parser = OptionParser(usage)
parser.add_option("-i", "--recopath",
                 dest="RECOPATH", help="Directory of the taupede version to run.")
parser.add_option("-v", "--verbose", default=False, action="store_true",
                 dest="VERBOSE", help="Print info level output.")
parser.add_option("-r", "--real", default=False, action="store_true",
                 dest="REAL", help="Real run.")
(options,args) = parser.parse_args()

if not options.RECOPATH:
    parser.error('Infile not specified.')

import os, shutil
import logging
import subprocess

logging.info('Options: {}'.format(options))

if options.VERBOSE:
    logging.getLogger().setLevel(level=logging.INFO)
else:
    logging.getLogger().setLevel(level=logging.WARNING)

RECOPATH = options.RECOPATH
logging.info('{}'.format(RECOPATH))
if not os.path.exists(RECOPATH):
    raise RuntimeError('Reco Path does not exist.')

RECOSCRIPT = os.path.join(RECOPATH, 'monopod_and_taupede.py')
logging.info('{}'.format(RECOSCRIPT))
if not os.path.isfile(RECOSCRIPT):
    raise RuntimeError('Reco Script does not exist.')

simulation_infos_dicts = []

n_files = {
    12:'645',
    14:'1550',
    16:'350',
}

max_runtimes = {
    12:'23:59:59',
    14:'23:59:59',
    16:'47:59:59',
}

# add info about genie nominal sets
for genie_set in [12,14,16]:
    simulation_info = dict()
    simulation_info['simulation_set'] = '{}0000'.format(genie_set)
    simulation_info['inpath'] = '/lustre/fs22/group/icecube/lfischer/data/OscNext/pass2/genie/level7_v02.00/{}0000'.format(genie_set)
    simulation_info['filetype'] = 'genie'
    simulation_info['n_jobs'] = n_files[genie_set]
    simulation_info['max_runtime'] = max_runtimes[genie_set]
    simulation_infos_dicts.append(simulation_info)

# add info about muongun nominal sets
simulation_info = dict()
simulation_info['simulation_set'] = '130000'
simulation_info['inpath'] = '/lustre/fs22/group/icecube/lfischer/data/OscNext/pass2/muongun/level7_v02.00/130000'
simulation_info['filetype'] = 'muongun'
simulation_info['n_jobs'] = '19999'
simulation_info['max_runtime'] = '00:29:59'
simulation_infos_dicts.append(simulation_info)

# add info about noise nominal sets

simulation_info = dict()
simulation_info['simulation_set'] = '888003'
simulation_info['inpath'] = '/lustre/fs22/group/icecube/lfischer/data/OscNext/pass2/noise/level7_v02.00/888003'
simulation_info['filetype'] = 'noise'
simulation_info['n_jobs'] = '10000'
simulation_info['max_runtime'] = '00:29:59'
simulation_infos_dicts.append(simulation_info)

path_of_this_script = os.path.split(os.path.abspath(__file__))[0]
CLUSTER_SCRIPT = os.path.join(path_of_this_script, 'cluster_bash_script.sh')


for simulation_info in simulation_infos_dicts:

    INPATH = simulation_info['inpath']

    SIMULATION_SET = simulation_info['simulation_set']

    FILETYPE = simulation_info['filetype']

    OUTPATH = os.path.join(RECOPATH, SIMULATION_SET)

    LOGPATH = os.path.join(OUTPATH, 'LOGS')

    # N_JOBS = 1 # for testing
    N_JOBS = simulation_info['n_jobs']

    MAX_RUNTIME = simulation_info['max_runtime']

    logging.info('Running taupede on {1} file(s) of OscNext {0} from {3} using the script in {2}'.format(SIMULATION_SET, N_JOBS, RECOPATH, INPATH))

    qsub_line = ['qsub', '-t', '1-{}'.format(N_JOBS), '-l', 'h_rt={}'.format(MAX_RUNTIME), '-o', '{}'.format(LOGPATH), '{}'.format(CLUSTER_SCRIPT), '{}'.format(INPATH), '{}'.format(RECOPATH), '{}'.format(SIMULATION_SET), '{}'.format(FILETYPE)]
    logging.info(qsub_line)

    if options.REAL:
        if not os.path.exists(OUTPATH):
            logging.info('Creating destination directory: {}'.format(OUTPATH))
            os.mkdir(OUTPATH)

        if not os.path.exists(LOGPATH):
            logging.info('Creating destination directory: {}'.format(LOGPATH))
            os.mkdir(LOGPATH)

        logging.info('Submitting jobs, for real now.')
        subprocess.call(qsub_line)

    # break

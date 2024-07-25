#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/icetray-start
#METAPROJECT /cvmfs/icecube.opensciencegrid.org/users/lfischer/oscnext_hnl_meta/build

#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# Author: Leander Fischer
#

from optparse import OptionParser

usage = "usage: %prog [options]"
parser = OptionParser(usage)
parser.add_option("-s", "--script", default='/data/user/lfischer/hnl_analysis/reconstruction/taupede_directional_seeds.py',
                 dest="SCRIPT", help="script location")
parser.add_option("-i", "--inpath", default='/data/ana/BSM/HNL/MC/190607/Ares/IC86.AVG/L7/domeff_0.97/',
                 dest="INPATH", help="input base location of the Level 7 Simulation")
parser.add_option("-o", "--outpath", default='/data/ana/BSM/HNL/MC/190607/Ares/IC86.AVG/Taupde_Directional_Seeding/domeff_0.97/',
                 dest="OUTPATH", help="iutput base location of the Taupede Level Simulation")
parser.add_option("-r", "--real", default=False, action="store_true",
                 dest="REAL", help="do the real run of all files")
parser.add_option("-v", "--verbose", default=False, action="store_true",
                 dest="VERBOSE", help="print info level output")
(options,args) = parser.parse_args()

import os
import numpy as np
import subprocess

import icecube.icetray.i3logging as logging

if options.VERBOSE:
    logging.set_level('INFO')
else:
    logging.set_level('WARN')

logging.log_info('Options: {}'.format(options))

executable = options.SCRIPT
in_basedir = options.INPATH

out_basedir = options.OUTPATH

if not os.path.exists(out_basedir):
    logging.log_info('Creating output base directory: {}'.format(out_basedir))
    os.makedirs(out_basedir)

dagman_location = '/scratch/lfischer/taupede/dagman'

"""
Write submission script
"""

# write submission script
submission_script_name = 'dagman_submit_taupede_directional_seeding_hnl_signal.sub'
submission_script_path = os.path.join(dagman_location, submission_script_name)

logging.log_info('Writing submission script.')
submission_script = open(submission_script_path, 'w')

logfolder = '/scratch/lfischer/condor_logs/taupede/'
if not os.path.isdir(logfolder):os.makedirs(logfolder)

outfolder = '/scratch/lfischer/condor_output/taupede/'
if not os.path.isdir(outfolder):os.makedirs(outfolder)

errfolder = '/scratch/lfischer/condor_error/taupede/hnl_signal/'
if not os.path.isdir(errfolder):os.makedirs(errfolder)

submission_script.write('# dagman submit taupede directional seeding reco for all hnl signal files in one go')
submission_script.write('\n')
submission_script.write('\n')
submission_script.write('executable = {}'.format(executable))
submission_script.write('\n')
submission_script.write('\n')
submission_script.write('arguments = -i $(i3_file_in_path) -o $(i3_file_out_path) -v')
submission_script.write('\n')
submission_script.write('\n')
submission_script.write('log = {}'.format(os.path.join(logfolder, 'hnl_signal.log')))
submission_script.write('\n')
submission_script.write('output = {}'.format(os.path.join(outfolder, 'hnl_signal.out')))
submission_script.write('\n')
submission_script.write('error = {}'.format(os.path.join(errfolder, '$(Cluster)_$(Process).err')))
submission_script.write('\n')
submission_script.write('\n')
submission_script.write('should_transfer_files = true')
submission_script.write('\n')
submission_script.write('\n')
submission_script.write('request_memory = 2GB')
submission_script.write('\n')
submission_script.write('\n')
submission_script.write('queue')
submission_script.write('\n')

submission_script.close()
logging.log_info('Submission Script written. ({})'.format(submission_script_path))

"""
Loop over directories and files and write the dagman script
"""

# write job script
job_script_name = 'dagman_job_taupede_directional_seeding_hnl_signal.dag'
job_script_path = os.path.join(dagman_location, job_script_name)

logging.log_info('Writing job script.')
job_script = open(job_script_path, 'w')

job_script.write('# dagman job file to submit hnl signal set')
job_script.write('\n')

filenumber = 0

for directory in sorted(os.listdir(in_basedir)):
    in_dir = os.path.join(in_basedir, directory)
    logging.log_info('{}'.format(in_dir))

    out_dir = os.path.join(out_basedir, directory)
    if not os.path.exists(out_dir):
        logging.log_info('Creating output directory: {}'.format(out_dir))
        os.makedirs(out_dir)

    file_list = os.listdir(in_dir)
    for i3_file_in in sorted(file_list):

        filenumber += 1

        if not i3_file_in.endswith('.i3.zst'):continue

        i3_file_in_path = os.path.join(in_dir, i3_file_in)
        if not options.REAL:
            logging.log_info('{}'.format(i3_file_in_path))

        i3_file_out = "Taupede_Directional_Seeding" + i3_file_in.split('L7')[-1]
        i3_file_out_path = os.path.join(out_dir, i3_file_out)

        job_script.write('JOB\tjob{}\t{}'.format(filenumber, submission_script_path))
        job_script.write('\n')
        job_script.write('VARS\tjob{}\ti3_file_in_path="{}"\ti3_file_out_path="{}"'.format(filenumber, i3_file_in_path, i3_file_out_path))
        job_script.write('\n')


        if not options.REAL:
            break
    if not options.REAL:
        break

logging.log_info('Written {} jobs for this set'.format(filenumber))

job_script.close()
logging.log_info('Job Script written. ({})'.format(job_script_path))
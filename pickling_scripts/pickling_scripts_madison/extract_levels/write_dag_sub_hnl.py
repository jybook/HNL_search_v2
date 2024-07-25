#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# Author: Leander Fischer
#

import os
import logging

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

levels = [
    # 'Gen',
    # 'Phot',
    # 'Det',
    # 'L1',
    # 'L2',
    # 'L3',
    # 'L4',
    # 'L5',
    # 'L6',
    'L7',
    'Taupede',
    'Millipede',
]

dagman_location = '/scratch/lfischer/dagman'
submission_script = '/data/user/lfischer/hnl_analysis/pickling_scripts/pickling_scripts_madison/extract_levels/pickle_hnl.sub'

# write job script
job_script_name = 'dagman_submit_190607_all_levels.dag'
job_script_path = os.path.join(dagman_location, job_script_name)

log.info('Writing job script.')
job_script = open(job_script_path, 'w')

job_script.write('# dagman job file to pickle all 190607 levels')
job_script.write('\n')

jobnumber = 1

for level in levels:
    log.info('Level: {}'.format(level))

    job_script.write('JOB\tjob{}\t{}'.format(jobnumber, submission_script))
    job_script.write('\n')
    job_script.write('VARS\tjob{}\tlevel="{}"'.format(jobnumber, level))
    job_script.write('\n')

    jobnumber += 1

job_script.close()
log.info('Job Script written: {}'.format(job_script_path))
log.info('To submit enter: "condor_submit_dag {}'.format(job_script_path))

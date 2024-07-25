#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# Author: Leander Fischer
#

import os
import logging

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

mc_sets = {
    # 'genie':['120000', '140000', '160000'],
    'genie':['160000'],
    # 'muongun':['130000'],
    # 'noise':['888003']
}

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
    # 'L7',
    # 'Taupede',
    'Millipede',
]

dagman_location = '/scratch/lfischer/dagman'
submission_script = '/data/user/lfischer/hnl_analysis/pickling_scripts/pickling_scripts_madison/extract_levels/pickle_oscnext.sub'

# write job script
# job_script_name = 'dagman_submit_oscnext_all_levels.dag'
# job_script_name = 'dagman_submit_oscnext_all_levels_Tau_Mil.dag'
# job_script_name = 'dagman_submit_oscnext_130000_L3_reduced.dag'
# job_script_name = 'dagman_submit_oscnext_130000_Millipede_reduced.dag'
# job_script_name = 'dagman_submit_oscnext_888003_Millipede_reduced.dag'
job_script_name = 'dagman_submit_oscnext_160000_Millipede_reduced.dag'
job_script_path = os.path.join(dagman_location, job_script_name)

log.info('Writing job script.')
job_script = open(job_script_path, 'w')

job_script.write('# dagman job file to pickle all oscnext sets and levels')
job_script.write('\n')

jobnumber = 1

for settype, setnumbers in mc_sets.iteritems():
    for setnumber in setnumbers:
        for level in levels:
            log.info('Level: {}'.format(level))

            job_script.write('JOB\tjob{}\t{}'.format(jobnumber, submission_script))
            job_script.write('\n')
            job_script.write('VARS\tjob{}\tlevel="{}"\ttype="{}"\tnumber="{}"'.format(jobnumber, level, settype, setnumber))
            job_script.write('\n')

            jobnumber += 1

job_script.close()
log.info('Job Script written: {}'.format(job_script_path))
log.info('To submit enter: "condor_submit_dag {}'.format(job_script_path))

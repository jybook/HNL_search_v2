#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# Author: Leander Fischer
#

import  os
import logging

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

mc_sets = {
    'genie':['120000', '140000', '160000'],
    'muongun':['130000'],
    'noise':['888003']
}

executable = '/data/user/lfischer/hnl_analysis/reconstruction/millipede.py'
dagman_location = '/scratch/lfischer/millipede/dagman'
input_base = '/data/ana/BSM/HNL/MC/oscNext'
output_base = input_base

for mc_set_name, mc_set_numbers in mc_sets.iteritems():
    for mc_set_number in mc_set_numbers:
        log.info('MC Set: {} {}'.format(mc_set_name, mc_set_number))

        input_files_location = os.path.join(input_base, mc_set_name)
        input_files_location = os.path.join(input_files_location, 'taupede')
        input_files_location = os.path.join(input_files_location, mc_set_number)
        log.info('Input Files Location: {}'.format(input_files_location))

        output_files_location = os.path.join(output_base, mc_set_name)
        output_files_location = os.path.join(output_files_location, 'millipede')
        output_files_location = os.path.join(output_files_location, mc_set_number)
        log.info('Output Files Location: {}'.format(output_files_location))


        # write submission script
        submission_script_name = 'dagman_submit_millipede_oscnext_{}_{}.sub'.format(mc_set_name, mc_set_number)
        submission_script_path = os.path.join(dagman_location, submission_script_name)
        
        log.info('Writing submission script.')
        submission_script = open(submission_script_path, 'w')

        submission_script.write('# dagman submit millipede reco for all {} {} files in one go'.format(mc_set_name, mc_set_number))
        submission_script.write('\n')
        submission_script.write('\n')
        submission_script.write('executable = {}'.format(executable))
        submission_script.write('\n')
        submission_script.write('\n')
        submission_script.write('arguments = -i $(file) -o {} -v'.format(os.path.join(output_files_location, '$Fnx(file)')))
        submission_script.write('\n')
        submission_script.write('\n')
        submission_script.write('log = /scratch/lfischer/condor_logs/millipede/oscnext_{}_{}.log'.format(mc_set_name, mc_set_number))
        submission_script.write('\n')
        submission_script.write('output = /scratch/lfischer/condor_output/millipede/oscnext_{}_{}.out'.format(mc_set_name, mc_set_number))
        submission_script.write('\n')
        submission_script.write('error = /scratch/lfischer/condor_error/millipede/{}/$(Cluster)_$(Process).err'.format(mc_set_number))
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
        log.info('Submission Script written. ({})'.format(submission_script_path))


        # write job script
        job_script_name = 'dagman_job_millipede_oscnext_{}_{}.dag'.format(mc_set_name, mc_set_number)
        job_script_path = os.path.join(dagman_location, job_script_name)
        
        log.info('Writing job script.')
        job_script = open(job_script_path, 'w')

        job_script.write('# dagman job file to submit {} {} set'.format(mc_set_name, mc_set_number))
        job_script.write('\n')

        filenumber = 0
        for filename in sorted(os.listdir(input_files_location)):
            if filename.endswith('i3.zst'):
                filenumber += 1
                # log.info('File: {}'.format(filename))

                job_script.write('JOB\tjob{}\t{}'.format(filenumber, submission_script_path))
                job_script.write('\n')
                job_script.write('VARS\tjob{}\tfile="{}"'.format(filenumber, os.path.join(input_files_location, filename)))
                job_script.write('\n')

                # break
            
        log.info('Written {} jobs for this set'.format(filenumber))

        job_script.close()
        log.info('Job Script written. ({})'.format(job_script_path))

    #     break
    # break

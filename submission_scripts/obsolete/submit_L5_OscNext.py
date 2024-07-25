#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Script to create dag/sub files to process background simulation with bugfixed L5

Author: Leander Fischer
'''

import  os
import logging

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

mc_sets = {
    # 'genie':['120000', '140000', '160000'],
    'muongun':['130000'],
    # 'noise':['888003']
}

executable = '/data/user/lfischer/hnl_analysis/submission_scripts/jobs/run_oscNext.py'
# dagman_location = '/scratch/lfischer/oscnext_L5/dagman'
dagman_location = '/scratch/lfischer/oscnext_L5_bugfix/dagman'
input_base = '/data/ana/LE/oscNext/pass2/'
output_base = '/data/ana/BSM/HNL/MC/oscNext/'
output_tmp = '/data/ana/BSM/HNL/MC/oscNext/tmp/'
gcd_file = '/cvmfs/icecube.opensciencegrid.org/data/GCD/GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE_withScaledNoise.i3.gz'


for mc_set_name, mc_set_numbers in mc_sets.iteritems():
    for mc_set_number in mc_set_numbers:
        log.info('MC Set: {} {}'.format(mc_set_name, mc_set_number))

        input_files_location = os.path.join(input_base, mc_set_name)
        input_files_location = os.path.join(input_files_location, 'level4')
        input_files_location = os.path.join(input_files_location, mc_set_number)
        log.info('Input Files Location: {}'.format(input_files_location))

        output_files_location = os.path.join(output_base, mc_set_name)
        output_files_location = os.path.join(output_files_location, 'L5_bugfix')
        if not os.path.isdir(output_files_location):os.makedirs(output_files_location)
        output_files_location = os.path.join(output_files_location, mc_set_number)
        if not os.path.isdir(output_files_location):os.makedirs(output_files_location)
        log.info('Output Files Location: {}'.format(output_files_location))


        # write submission script
        submission_script_name = 'dagman_submit_oscnext_L5_bugfix_{}_{}.sub'.format(mc_set_name, mc_set_number)
        submission_script_path = os.path.join(dagman_location, submission_script_name)
        
        log.info('Writing submission script.')
        submission_script = open(submission_script_path, 'w')

        logfolder = '/scratch/lfischer/condor_logs/l5_bugfix/'
        if not os.path.isdir(logfolder):os.makedirs(logfolder)

        outfolder = '/scratch/lfischer/condor_output/l5_bugfix/'
        if not os.path.isdir(outfolder):os.makedirs(outfolder)

        errfolder = '/scratch/lfischer/condor_error/l5_bugfix/{}'.format(mc_set_number)
        if not os.path.isdir(errfolder):os.makedirs(errfolder)

        submission_script.write('# dagman submit oscnext l5 bugfix for all {} {} files in one go'.format(mc_set_name, mc_set_number))
        submission_script.write('\n')
        submission_script.write('\n')
        submission_script.write('executable = {}'.format(executable))
        submission_script.write('\n')
        submission_script.write('\n')
        submission_script.write('arguments = -t {} -i $(file) -o {} -g {} -l 5 --tmp-dir {}'.format(mc_set_name, os.path.join(output_files_location, '$Fnx(file)'), gcd_file, output_tmp))
        submission_script.write('\n')
        submission_script.write('\n')
        submission_script.write('log = {}'.format(os.path.join(logfolder, 'oscnext_{}_{}.log'.format(mc_set_name, mc_set_number))))
        submission_script.write('\n')
        submission_script.write('output = {}'.format(os.path.join(outfolder, 'oscnext_{}_{}.out'.format(mc_set_name, mc_set_number))))
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
        log.info('Submission Script written. ({})'.format(submission_script_path))

        # write job script
        job_script_name = 'dagman_job_oscnext_L5_bugfix_{}_{}.dag'.format(mc_set_name, mc_set_number)
        job_script_path = os.path.join(dagman_location, job_script_name)
        
        log.info('Writing job script.')
        job_script = open(job_script_path, 'w')

        job_script.write('# dagman job file to submit {} {} set'.format(mc_set_name, mc_set_number))
        job_script.write('\n')

        n_files = 200
        filenumber = 1
        for filename in sorted(os.listdir(input_files_location)):
            if filename.endswith('i3.zst'):
                # log.info('File: {}'.format(filename))

                job_script.write('JOB\tjob{}\t{}'.format(filenumber, submission_script_path))
                job_script.write('\n')
                job_script.write('VARS\tjob{}\tfile="{}"'.format(filenumber, os.path.join(input_files_location, filename)))
                job_script.write('\n')

                if filenumber ==n_files: break
                filenumber += 1
            
        log.info('Written {} jobs for this set'.format(filenumber))

        job_script.close()
        log.info('Job Script written. ({})'.format(job_script_path))

    #     break
    # break

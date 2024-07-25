#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/icetray-start
#METAPROJECT /data/user/lfischer/software/oscnext_hnl/build/

'''
Author: Leander Fischer

Script to produce dagman file to process HNL simulation with Taupede.
'''

import glob, os, sys, random
import numpy as np
from shutil import copyfile
from optparse import OptionParser
from os.path import expandvars

def main():
    '''
    Main function: Write dag file to submit Taupede level jobs to condor.
    '''

    usage = "usage: %prog [options] inputfile"
    parser = OptionParser(usage)
    parser.add_option("--limit",type="int",
                    dest="limit", help="How many files do you want?. [default: %default]")
    parser.add_option("-i", "--identifier",type="string",
                    dest="identifier", help="Set name (infiles). [default: %default]")
    parser.add_option("-o", "--identifier_out",type="string",
                    dest="identifier_out", help="Set name (outfiles). [default: %default]")
    parser.add_option("-l", "--processing_level",type="string", default='L8',
                    dest="level", help="Which processing level to run on [default: %default]")
    parser.add_option("-g", "--gcdfile",
                    default='/cvmfs/icecube.opensciencegrid.org/data/GCD/GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE_withScaledNoise.i3.gz',
                    help="Read in GCD file [default: %default]")
    parser.add_option("-s", "--seedparticle", default='FlercnnSeedParticle',
                    dest="SEEDPARTICLE", help="Choose particle to use as seed for taupede, either real I3Particle or from ('RetroSeedParticle', 'TruthSeedParticle' or 'FiniteRecoFit+OnlineL2_SplineMPE_DirectHitsC_DirTrackLength') [default: %default]")
    parser.add_option("-b", "--BrightDOMThreshold",dest="BrightDOMThreshold",type="int",
                    help="DOMs with charge larger than this factor times the median will be excluded from fit [default: %default]", default=30)
    parser.add_option("-m","--icemodel", default="spice_lea",
                    dest="ICEMODEL",help="Ice model (spice_3.2.1, spice_lea etc.) [default: %default]")
    parser.add_option("-p", "--pulseseries", default='SplitInIcePulses',
                    dest="PULSESERIES", help="Choose pulse series to use [default: %default]")
    parser.add_option("-r", "--REALRUN", default=False, action="store_true",
                    dest="REALRUN", help="Do real run. (otherwise just test file) [default: %default]")
    parser.add_option("--mc_location",type="string", default='/data/ana/BSM/HNL/MC/',
                    dest="mc_location", help="The location to write the MC files to [default: %default]")
    parser.add_option("--dag_location",type="string", default='/scratch/lfischer/dagman/',
                    dest="dag_location", help="The location to write the dag files to [default: %default]")
    parser.add_option("--real_sim", type="int",
                    dest="real_sim", default=1, help="Is this the real (model dependent) simulation? [default: %default]")
    (options,args) = parser.parse_args()


    if not options.identifier_out:
        options.identifier_out = options.identifier

    limit 		          = options.limit
    retries 	          = '2'
    gcdfile               = options.gcdfile
    brightdomthreshold    = options.BrightDOMThreshold
    icemodel		      = options.ICEMODEL
    seedparticle          = options.SEEDPARTICLE
    pulseseries           = options.PULSESERIES
    identifier 	          = options.identifier
    identifier_out	      = options.identifier_out
    level	              = options.level
    realrun 		      = options.REALRUN
    mc_location 	      = options.mc_location
    dag_location	      = options.dag_location
    real_sim	          = options.real_sim

    mother_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    submit_location = os.path.join(mother_path, 'submit')

    if not os.path.exists(dag_location):
        os.makedirs(dag_location)

    if not os.path.exists(gcdfile):
        print('Missing GCD file: ' + gcdfile)
        sys.exit()
    print('GCD file: '+gcdfile)

    outdir_base = mc_location + identifier_out
    if not os.path.exists(outdir_base):
        os.makedirs(outdir_base)

    ###### Write taupede information file ######
    counter = 0
    taupede_data_file = outdir_base+"/simulation_log_Taupede_{}_v{}.txt".format(identifier_out, counter)
    while(True):
        if os.path.isfile(taupede_data_file):
            counter +=1
            taupede_data_file = taupede_data_file.split('_v')[0] + '_v{}.txt'.format(counter)
        else:
            break

    text_file = open(taupede_data_file, "w")
    text_file.write("Identifier: " + str(identifier) + '\n')
    text_file.write("Identifier Out: " + str(identifier_out) + '\n')
    text_file.write("Running on level: " + str(level) + '\n')
    text_file.write("Limit: " + str(limit) + '\n')
    text_file.write("GCD file: " + str(gcdfile) + '\n')
    text_file.write("Brightdomthreshold: " + str(brightdomthreshold) + '\n')
    text_file.write("Icemodel: " + str(icemodel) + '\n')
    text_file.write("Seedparticle: " + str(seedparticle) + '\n')
    text_file.write("Pulseseries: " + str(pulseseries) + '\n')
    text_file.write('\n')
    text_file.close()


    ###### Write dag file ######
    counter = 0
    subf = os.path.join(dag_location, 'dagman-Taupede_{}_v{}.dag'.format(identifier_out, counter))

    while(True):
        if os.path.isfile(subf):
            counter +=1
            subf = subf.split('_v')[0] + '_v{}.dag'.format(counter)
        else:
            break

    subtarget = open(subf, 'w')
    os.chmod(subf, 0o775)

    processing_dict = {}
    processing_dict['GCD'] = gcdfile

    processing_dict['Taupede'] = {}
    processing_dict['Taupede']['infiles'] = []
    processing_dict['Taupede']['outfiles'] = []
    processing_dict['Taupede']['submit_location'] = submit_location + '/submit-Taupede.condor'

    infile_list = []

    folders = sorted(glob.glob(mc_location+identifier.split('/')[0]+'/{}/domeff**'.format(level)+'/[0-9]*'))
    # print(folders)

    for folder in folders:
        folder_name = os.path.basename(folder)
        # print(folder_name)
        files = sorted(glob.glob(folder + '/*.i3.zst'))
        if len(files) == 0:
            files = sorted(glob.glob(folder + '/*.i3.bz2'))
        for file in files:
            infile_list.append(file)

        domeff = infile_list[0].split('domeff')[-1].split('/')[0]

        folder_out_dir = os.path.join(os.path.join(os.path.join(outdir_base, 'Taupede'), 'domeff{}'.format(domeff)), folder_name)
        # print(folder_out_dir)
        if not os.path.exists(folder_out_dir):
            os.makedirs(folder_out_dir)
        os.chmod(folder_out_dir, 0o775)

    print('Number of {} files: '.format(level) + str(len(infile_list)))

    if not limit:
        limit = len(infile_list)

    for infile in infile_list:
        # print(infile)
        processing_dict['Taupede']['infiles'].append(infile)
        infile_base = os.path.basename(infile)
        folder = infile.split('/')[-2]
        # print(folder)

        Taupede_outfile = os.path.join(os.path.join(os.path.join(os.path.join(outdir_base, 'Taupede'), 'domeff{}'.format(domeff)), folder), infile_base.replace(level, 'Taupede'))

        processing_dict['Taupede']['outfiles'].append(Taupede_outfile)

        if not realrun:break

    def print_Taupede():
        subtarget.write('JOB\tjob'+str(jobName)+'\t' + processing_dict['Taupede']['submit_location'])
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(jobName)+'\tgcdfile="' + processing_dict['GCD'] + '"')
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(jobName)+'\tinfile="' + processing_dict['Taupede']['infiles'][job_counter] + '"')
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(jobName)+'\toutfile="' + processing_dict['Taupede']['outfiles'][job_counter] +'"')
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(jobName)+'\tidentifier_out="'+identifier_out+'"')
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(jobName)+'\tbrightdomthreshold="'+str(brightdomthreshold)+'"')
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(jobName)+'\ticemodel="'+str(icemodel)+'"')
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(jobName)+'\tseedparticle="'+seedparticle+'"')
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(jobName)+'\tpulseseries="'+pulseseries+'"')
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(jobName)+'\treal_sim="'+str(real_sim)+'"')
        subtarget.write("\n")
        subtarget.write('RETRY\tjob'+str(jobName)+'\t'+retries)
        subtarget.write("\n")

    job_counter = 0
    job_list = []

    while job_counter < int(limit) and job_counter < len(processing_dict['Taupede']['outfiles']):
        if not os.path.isfile(processing_dict['Taupede']['outfiles'][job_counter]) or os.path.getsize(processing_dict['Taupede']['outfiles'][job_counter]) < 100000:
            jobName = str(job_counter+1) + 'a'
            print_Taupede()
            subtarget.write("\n")
            job_list.append(processing_dict['Taupede']['outfiles'][job_counter])
            job_counter +=1
        else:
            print('File already processed: '+ processing_dict['Taupede']['outfiles'][job_counter])
            # print(job_counter)
            job_counter +=1

    subtarget.close()
    ###### End ######

    print('Jobs setup: {}'.format(len(job_list)))
    if job_counter == 0:
        os.remove(subf)
        os.remove(taupede_data_file)
        print('No jobs, removing dag/info file')
    else:
        print('Dag file: {}'.format(subf))
        print('Submit with:\ncondor_submit_dag {}'.format(subf))
    print('done...')


if __name__ == '__main__':
    main()

#!eval $(/cvmfs/icecube.opensciencegrid.org/py3-v4.2.1/setup.sh)
#METAPROJECT /data/user/jbook/i3/build

#=================================
# S. Axani
#
# Modified by D Vanerrom, L Fischer, J Book
#=================================

import os, sys
import numpy as np

from optparse import OptionParser

def main():
    '''
    Main function: Write dag file to submit Generation level jobs to condor.
    '''

    usage = "usage: %prog [options] inputfile"
    parser = OptionParser(usage)
    parser.add_option("-i", "--identifier",type="string",
            dest="identifier", help="HNL signal set number (19XXXX) [default: %default]")
    parser.add_option("--ram",type="string", default="1000",
            help="The amount of ram to use (in MB) for each job [default: %default], minimum 1GB")
    parser.add_option("--Emin",type="float", default="2",
            dest="Emin", help="The minimum energy to simulate [default: %default]")
    parser.add_option("--Emax",type="int", default="10000",
            dest="Emax", help="The maximum energy to simulate [default: %default]")
    parser.add_option("--index",type="float", default="2.",
            dest="index", help="The index to simulate [default: %default]")
    parser.add_option("-m", "--HNL_mass",type="float",
            dest="HNL_mass", help="The mass of the HNL [default: %default]")
    parser.add_option("--nFiles",type="int", default="10000",
            dest="nFiles", help="The total number of files you want to generate [default: %default]")
    parser.add_option("--nEvents",type="int", default="100000",
            dest="nEvents", help="The number of events in each file [default: %default]")
    parser.add_option("--capacity",type="int", default="1000",
            dest="capacity", help="The number of files in each folder [default: %default]")
    parser.add_option("--limit",type="string",
            dest="limit", help="set to 1-10000 to setup files from 1-10000 [default: %default]") 
    parser.add_option("--radius",type="string", default="600",
            dest="radius", help="The injection radius [default: %default]") 
    parser.add_option("--length",type="string", default="600",
            dest="length", help="The cylinder length of event genreration [default: %default]")
    parser.add_option("--mc_location",type="string", default='/data/ana/BSM/HNL/MC/',
            dest="mc_location", help="The location to write the MC files to [default: %default]")
    parser.add_option("--dag_location",type="string", default='/scratch/jbook/dagman/',
            dest="dag_location", help="The location to write the dag files to [default: %default]")
    parser.add_option("--username",type="string",
                    dest="username", help="Set username (needed for OSG support) [default: %default]")
    parser.add_option("--osg", default=False, action="store_true",
                     dest="osg", help="Is this job running on the OSG? [default: %default]")
    (options,args) = parser.parse_args()

    if not options.username:
        parser.error('Username must be specified with --username option.')
    if not options.identifier:
        parser.error("Identifier not given")
    if not options.HNL_mass:
        parser.error("HNL mass not given")
                

    identifier		= options.identifier
    ram 			= options.ram
    Emin 			= options.Emin
    Emax 			= options.Emax
    index 			= options.index
    HNL_mass 		= options.HNL_mass
    nFiles 			= options.nFiles
    nEvents 		= options.nEvents
    capacity 		= options.capacity
    limit 			= options.limit
    radius 			= options.radius
    length 			= options.length
    mc_location 	= options.mc_location
    dag_location	= options.dag_location
    osg             = options.osg
    username        = options.username

    if limit == None:
        limit = '1-'+str(nFiles)

    mother_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    submit_location = os.path.join(mother_path, 'submit')
    submit_location_osg = '/scratch/{}/submission/'.format(username)
    process_location = os.path.join(mother_path, 'process')

    if not os.path.exists(dag_location):
        os.makedirs(dag_location)

    outdir_base = mc_location +'/' + identifier+'/Gen/'
    if not os.path.exists(outdir_base):
        os.makedirs(outdir_base)
        os.chmod(outdir_base, 0o775)


    ###### Write Gen information file ######
    counter = 0
    gen_data_file = outdir_base+"/Generation_data_{}_v{}.txt".format(identifier, counter)
    while(True):
        if os.path.isfile(gen_data_file):
            counter +=1
            gen_data_file = gen_data_file.split('_v')[0] + '_v{}.txt'.format(counter)
        else:
            break

    with open(gen_data_file, "w") as text_file:
        text_file.write("{0:<25} {1} \n".format('Identifier:', identifier))
        text_file.write("{0:<25} {1} \n".format('ram:', ram))
        text_file.write("{0:<25} {1} \n".format('Emin:', Emin))
        text_file.write("{0:<25} {1} \n".format('Emax:', Emax))
        text_file.write("{0:<25} {1} \n".format('index:', index))
        text_file.write("{0:<25} {1} \n".format('HNL_mass:', HNL_mass))
        text_file.write("{0:<25} {1} \n".format('nFiles:', nFiles))
        text_file.write("{0:<25} {1} \n".format('limit:', limit))
        text_file.write("{0:<25} {1} \n".format('nEvents per file:', nEvents))
        text_file.write("{0:<25} {1} \n".format('capacity:', capacity))
        text_file.write("{0:<25} {1} \n".format('radius:', radius))
        text_file.write("{0:<25} {1} \n".format('length:', length))
        text_file.write("{0:<25} {1} \n".format('Using OSG?', osg))
        text_file.write("{0:<25} {1} \n".format('Who ran this?', username))

    ###### Write dag file ######
    counter = 0
    subf_base = 'dagman-Gen_{}_v{}.dag'.format(identifier, counter)
    subf = os.path.join(dag_location, subf_base)

    while(True):
        if os.path.isfile(subf):
            counter +=1
            subf = subf.split('_v')[0] + '_v{}.dag'.format(counter)
        else:
            break

    subf_base = os.path.basename(subf)

    subtarget = open(subf, 'w')
    os.chmod(subf, 0o775)

    smallest = limit.split('-')[0]
    biggest = limit.split('-')[1]

    eventlist = np.arange(int(smallest), int(biggest)+1)

    job_number = 1

    ###### Start writing dag file ######
    for seed in eventlist:
        # create folder
        if seed<0:
            print("Seed is negative, this really should not be.")
            sys.exit()
        low = str(int(np.floor((seed-1)/capacity))*capacity+1).zfill(5)
        high = str(int(np.floor((seed-1)/capacity))*capacity+capacity).zfill(5)
        folder = low+'-'+high
        outlocation = os.path.join(outdir_base,folder)
        if not os.path.exists(outlocation):
            os.makedirs(outlocation)
            os.chmod(outlocation, 0o775)

        outfile = os.path.join(outlocation, 'Gen_{:05d}.i3.zst'.format(seed))
        if os.path.exists(outfile):
            # print("Generation file {} already exists, skipping this file.".format(outfile))
            continue
        if osg:
            subtarget.write('JOB\tjob'+str(job_number)+'\t'+ '/scratch/{}/submission/submit-Gen-OSG.condor'.format(username))
        else:
            subtarget.write('JOB\tjob'+str(job_number)+'\t'+ submit_location + '/submit-Gen.condor')
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(job_number)+'\tseed="' + str(seed) + '"')  # og version
        # subtarget.write('VARS\tjob'+str(job_number)+'\tseed="' + str(seed+10000) + '"')  # add to seed to avoid duplication with long lengths sets
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(job_number)+'\tram="' + str(ram) + '"')
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(job_number)+'\tidentifier_out="' + str(identifier) + '"')
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(job_number)+'\tindex="' + str(index) + '"')
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(job_number)+'\tHNL_mass="' + str(HNL_mass) + '"')
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(job_number)+'\tEmin="' + str(Emin) + '"')
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(job_number)+'\tEmax="' + str(Emax) + '"')
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(job_number)+'\tnEvents="' + str(nEvents) + '"')
        subtarget.write("\n")
        if osg:
            subtarget.write('VARS\tjob'+str(job_number)+'\toutfile="' + 'gsiftp://gridftp-users.icecube.wisc.edu{}'.format(outfile) + '"')
        else:
            subtarget.write('VARS\tjob'+str(job_number)+'\toutfile="' + str(outfile) + '"')
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(job_number)+'\tradius="' + str(radius) + '"')
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(job_number)+'\tlength="' + str(length) + '"')
        subtarget.write("\n")
        job_number += 1
    subtarget.close()
    ###### End ######

    print('We setup: {} jobs'.format(job_number-1))
    if job_number - 1 == 0:
        os.remove(subf)
        os.remove(gen_data_file)
        print('No jobs, removing dag/info file')
    else:
        if osg:
            subf_on_sub1 = os.path.join('/scratch/{}/dagman'.format(username), subf_base)
            print("Copying dag/sub/process scripts to sub-2.")
            # print('scp {} sub-2:{}'.format(subf, subf_on_sub1))
            os.system('scp {} sub-2:{}'.format(subf, subf_on_sub1))
            
            # print('scp {} sub-2:{}'.format(submit_location + '/submit-Gen-OSG.condor', '/scratch/{}/submission/submit-Gen-OSG.condor'.format(username)))
            os.system('scp {} sub-2:{}'.format(submit_location + '/submit-Gen-OSG.condor', '/scratch/{}/submission/submit-Gen-OSG.condor'.format(username)))

            # print('scp {} sub-2:{}'.format(process_location + '/process_Gen.py', '/scratch/{}/submission/process_Gen-OSG.py'.format(username)))
            os.system('scp {} sub-2:{}'.format(process_location + '/process_Gen.py', '/scratch/{}/submission/process_Gen-OSG.py'.format(username)))
            print('Submit (from sub-2!) with: "condor_submit_dag {}"'.format(subf_on_sub1))
        else:
            print('Dag file: {}'.format(subf))
            print('Submit with: "condor_submit_dag {}"'.format(subf))
    print('done...')


if __name__ == '__main__':
    main()

#!/bin/sh /cvmfs/icecube.opensciencegrid.org/users/BeyondStandardModel/bsm-py2-v3.1.1-1/icetray-start-standard
#METAPROJECT /cvmfs/icecube.opensciencegrid.org/users/lfischer/oscnext_hnl_meta/build/

#=================================
# L. Fischer
# Questions? leander.fischer@desy.de
#=================================

import glob, os, sys
from optparse import OptionParser

def main():
    '''
    Main function: Write dag file to submit L9 level jobs to condor.
    '''

    usage = "usage: %prog [options] inputfile"
    parser = OptionParser(usage)
    parser.add_option("-e", "--domeff",type="string", default='1.00',
                    dest="domeff", help="Which dom efficiency are we using. (ie. '0.95','0.97','1.00','1.03','1.05','1.10','0.90', '2.00'). Just used for folder naming [default: %default]")
    parser.add_option("--limit",type="int",
                    dest="limit", help="How many files do you want to have processed (absolute number?). [default: %default]")
    parser.add_option("-i", "--identifier",type="string",
                    dest="identifier", help="Set name (infiles) [default: %default]")
    parser.add_option("-o", "--identifier_out",type="string",
                    dest="identifier_out", help="Set name (outfiles) [default: %default]")
    parser.add_option("-g", "--GCD",
                    default='/cvmfs/icecube.opensciencegrid.org/data/GCD/GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE_withScaledNoise.i3.gz',
                    help="Read in GCD file [default: %default]")
    parser.add_option("-r","--REALRUN", default=False, action="store_true",
                    dest="REALRUN", help="Do real run. (otherwise just test file) [default: %default]")
    parser.add_option("--mc_location",type="string", default='/data/ana/BSM/HNL/MC/',
                    dest="mc_location", help="The location to write the MC files to [default: %default]")
    parser.add_option("--dag_location",type="string", default='/scratch/lfischer/dagman/',
                    dest="dag_location", help="The location to write the dag files to [default: %default]")
    parser.add_option("--username",type="string",
                    dest="username", help="Set username (needed for OSG support) [default: %default]")
    parser.add_option("--osg", default=False, action="store_true",
                     dest="osg", help="Is this job running on the OSG, if so, you need to run this script on cobalt, not the submit node! [default: %default]")
    parser.add_option("--real_sim", type="int",
                    dest="real_sim", default=1, help="Is this the real (model dependent) simulation? [default: %default]")
    (options,args) = parser.parse_args()

    if not options.username:
        parser.error('Username must be specified with --username option.')
    if not options.identifier:
        parser.error('Set identifier must be specified with --identifier option.')
    if not options.identifier_out:
        options.identifier_out = options.identifier

    domeff          = options.domeff
    limit 		    = options.limit
    retries 	    = '2'
    GCD             = options.GCD
    identifier 	    = options.identifier
    identifier_out	= options.identifier_out
    realrun 		= options.REALRUN
    mc_location 	= options.mc_location
    dag_location	= options.dag_location
    osg             = options.osg
    username        = options.username
    real_sim	    = options.real_sim


    list_domeff = ['0.95', '0.97', '1.00', '1.03', '1.05', '1.10', '0.90', '2.00']

    print("Using dom eficiency: {}".format(domeff))
    if domeff in list_domeff:
        D = str(list_domeff.index(domeff))
    else:
        print('Bad dom eff')
        sys.exit()

    mother_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    submit_location = os.path.join(mother_path, 'submit')
    submit_location_osg = '/scratch/{}/submission/'.format(username)
    process_location = os.path.join(mother_path, 'process')

    if not os.path.exists(dag_location):
        os.makedirs(dag_location)

    if not os.path.exists(GCD):
        print('Missing GCD file: {}'.format(GCD))
        sys.exit()
    print('GCD file: {}'.format(GCD))

    outdir_base = mc_location + identifier_out
    if not os.path.exists(outdir_base):
        os.makedirs(outdir_base)

    ###### Write L9 information file ######
    counter = 0
    L9_data_file = outdir_base+"/simulation_log_L9_{}_v{}.txt".format(identifier, counter)
    while(True):
        if os.path.isfile(L9_data_file):
            counter +=1
            L9_data_file = L9_data_file.split('_v')[0] + '_v{}.txt'.format(counter)
        else:
            break

    text_file = open(L9_data_file, "w")
    text_file.write("Identifier: " + str(identifier) + '\n')
    text_file.write("Identifier Out: " + str(identifier_out) + '\n')
    text_file.write("DOM efficiency: " + str(domeff) + '\n')
    text_file.write("Limit: " + str(limit) + '\n')
    text_file.write("GCD file: " + str(GCD) + '\n')
    text_file.write("{0:<25} {1} \n".format('Using OSG?', osg))
    text_file.write("{0:<25} {1} \n".format('Who ran this?', username))

    text_file.write('\n')
    text_file.close()

    ###### Write dag file ######
    counter = 0
    subf_base = 'dagman-L9_{}_v{}.dag'.format(identifier_out, counter)
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

    processing_dict = {}
    processing_dict['GCD'] = GCD

    processing_dict['L9'] = {}
    processing_dict['L9']['outdirs'] = []
    processing_dict['L9']['outfiles'] = []
    processing_dict['L9']['outfiles_hdf5'] = []
    processing_dict['L9']['submit_location'] = os.path.join(submit_location, 'submit-L9.condor')

    infile_list = []

    folders = sorted(glob.glob(mc_location+identifier.split('/')[0]+'/L8/domeff'+domeff+'/[0-9]*'))
    # print(folders)

    level_list = ['L9']

    for folder in folders:
        folder_name = os.path.basename(folder)
        files = sorted(glob.glob(folder + '/*.i3.zst'))
        if len(files) == 0:
            files = sorted(glob.glob(folder + '/*.i3.bz2'))
        for file in files:
            infile_list.append(file)
        for level in level_list:
            folder_out_dir = os.path.join(os.path.join(os.path.join(outdir_base, level), 'domeff{}'.format(domeff)), folder_name)
            if not os.path.exists(folder_out_dir):
                os.makedirs(folder_out_dir)
            os.chmod(folder_out_dir, 0o775)

    print('Number of L8 files: {}'.format(len(infile_list)))

    if not limit:
        limit = len(infile_list)

    for infile in infile_list:
        infile_base = os.path.basename(infile)
        folder = infile.split('/')[-2]

        L9_outdir = os.path.join(os.path.join(os.path.join(outdir_base, 'L9'), 'domeff{}'.format(domeff)), folder)
        L9_outfile = infile_base.replace('L8', 'L9').replace('.i3.zst', '')

        if osg:
            processing_dict['L9']['outdirs'].append( 'gsiftp://gridftp-users.icecube.wisc.edu{}'.format(L9_outdir) )
        else:
            processing_dict['L9']['outdirs'].append(L9_outdir)

        processing_dict['L9']['outfiles'].append(L9_outfile)

        if not realrun:break

    def print_L9(i):
        if osg:
            subtarget.write('JOB\tjob'+str(jobName)+'\t' + processing_dict['L9']['submit_location_osg'])
            subtarget.write("\n")
            subtarget.write('VARS\tjob'+str(jobName)+'\tINFILE="' + 'gsiftp://gridftp-users.icecube.wisc.edu{}'.format(infile_list[i]) + '"')
        else:
            subtarget.write('JOB\tjob'+str(jobName)+'\t' + processing_dict['L9']['submit_location'])
            subtarget.write("\n")
            subtarget.write('VARS\tjob'+str(jobName)+'\tINFILE="' + infile_list[i] + '"')
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(jobName)+'\tOUTDIR="' + processing_dict['L9']['outdirs'][i] +'/"')
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(jobName)+'\tOUTFILE="' + processing_dict['L9']['outfiles'][i] +'"')
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(jobName)+'\tLOGNAME="' + processing_dict['L9']['outfiles'][i] +'"')
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(jobName)+'\tidentifier_out="' + identifier_out +'"')
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(jobName)+'\tREALSIM="'+str(real_sim)+'"')
        subtarget.write("\n")
        subtarget.write('RETRY\tjob'+str(jobName)+'\t'+retries)
        subtarget.write("\n")

    job_counter = 0
    job_list = []
    while job_counter < int(limit) and job_counter < len(infile_list):
        this_outfile = os.path.join(processing_dict['L9']['outdirs'][job_counter], processing_dict['L9']['outfiles'][job_counter] + '.i3.zst')
        if not os.path.isfile(this_outfile) or os.path.getsize(this_outfile) < 100000:

            jobName = str(job_counter+1) + 'a'
            print_L9(i=job_counter)
            job_list.append(this_outfile)

            job_counter +=1
        else:
            print('File already processed:  {}'.format(this_outfile))
            job_counter +=1

        if not realrun: break

    subtarget.close()
    ###### End ######

    print('Jobs setup: {}'.format(len(job_list)))
    if job_counter == 0:
        os.remove(subf)
        os.remove(L9_data_file)
        print('No jobs, removing dag/info file')
    else:
        if osg:
            subf_on_sub1 = os.path.join('/scratch/{}/dagman'.format(username), subf_base)
            print("Copying dag/sub/process scripts to sub-1.")
            os.system('scp {} sub-1:{}'.format(subf, subf_on_sub1))

            os.system('scp {} sub-1:{}'.format(os.path.join(submit_location, 'submit-L9-OSG.condor'), processing_dict['L9']['submit_location_osg']))

            os.system('scp {} sub-1:{}'.format(os.path.join(process_location, 'process_flercnn_L9-OGS.sh'), os.path.join(submit_location_osg, 'process_flercnn_L9-OSG.sh')))

            print('Submit (from sub-1!) with: "condor_submit_dag {}"'.format(subf_on_sub1))
        else:
            print('Dag file: {}'.format(subf))
            print('Submit with:\ncondor_submit_dag {}'.format(subf))
    print('done...')


if __name__ == '__main__':
    main()

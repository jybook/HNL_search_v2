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
    Main function: Write dag file to submit L6/L7/L8 level jobs to condor.
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

    ###### Write L6/L7/L8 information file ######
    counter = 0
    L6_L7_L8_data_file = outdir_base+"/simulation_log_L6_L7_L8_{}_v{}.txt".format(identifier, counter)
    while(True):
        if os.path.isfile(L6_L7_L8_data_file):
            counter +=1
            L6_L7_L8_data_file = L6_L7_L8_data_file.split('_v')[0] + '_v{}.txt'.format(counter)
        else:
            break

    text_file = open(L6_L7_L8_data_file, "w")
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
    subf_base = 'dagman-L6_L7_L8_{}_v{}.dag'.format(identifier_out, counter)
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

    processing_dict['L6'] = {}
    processing_dict['L6']['infiles'] = []
    processing_dict['L6']['outfiles'] = []
    processing_dict['L6']['lognames'] = []
    processing_dict['L6']['outfiles_hdf5'] = []
    processing_dict['L6']['submit_location_osg'] = os.path.join(submit_location_osg, 'submit-L6-OSG.condor')
    processing_dict['L6']['submit_location'] = os.path.join(submit_location, 'submit-L6.condor')

    processing_dict['L7'] = {}
    processing_dict['L7']['outdirs'] = []
    processing_dict['L7']['outfiles'] = []
    processing_dict['L7']['lognames'] = []
    processing_dict['L7']['outfiles_hdf5'] = []
    processing_dict['L7']['submit_location_osg'] = os.path.join(submit_location_osg, 'submit-L7-OSG.condor')
    processing_dict['L7']['submit_location'] = os.path.join(submit_location, 'submit-L7.condor')

    processing_dict['L8'] = {}
    processing_dict['L8']['outfiles'] = []
    processing_dict['L8']['lognames'] = []
    processing_dict['L8']['outfiles_hdf5'] = []
    processing_dict['L8']['submit_location_osg'] = os.path.join(submit_location_osg, 'submit-L8-OSG.condor')
    processing_dict['L8']['submit_location'] = os.path.join(submit_location, 'submit-L8.condor')


    infile_list = []

    folders = sorted(glob.glob(mc_location+identifier.split('/')[0]+'/L5/domeff'+domeff+'/[0-9]*'))
    # print(folders)

    level_list = ['L6','L7','L8']

    for folder in folders:
        folder_name = os.path.basename(folder)
        # print(folder_name)
        files = sorted(glob.glob(folder + '/*.i3.zst'))
        if len(files) == 0:
            files = sorted(glob.glob(folder + '/*.i3.bz2'))
        for file in files:
            infile_list.append(file)
        for level in level_list:
            folder_out_dir = os.path.join(os.path.join(os.path.join(outdir_base, level), 'domeff{}'.format(domeff)), folder_name)
            # print(folder_out_dir)
            if not os.path.exists(folder_out_dir):
                os.makedirs(folder_out_dir)
            os.chmod(folder_out_dir, 0o775)

    print('Number of L5 files: {}'.format(len(infile_list)))

    if not limit:
        limit = len(infile_list)

    for infile in infile_list:
        infile_base = os.path.basename(infile)
        folder = infile.split('/')[-2]

        # # modify infile/infile_base/outfile_base to work on the grid
        # if osg:
        #     infile = 'gsiftp://gridftp-users.icecube.wisc.edu{}'.format(infile)
        #     outdir_base = 'gsiftp://gridftp-users.icecube.wisc.edu{}'.format(outdir_base)

        L6_outfile = os.path.join(os.path.join(os.path.join(os.path.join(outdir_base, 'L6'), 'domeff{}'.format(domeff)), folder), infile_base.replace('L5', 'L6'))
        L6_logname = infile_base.replace('L5', 'L6').replace('.i3.zst', '')

        L7_outdir = os.path.join(os.path.join(os.path.join(outdir_base, 'L7'), 'domeff{}'.format(domeff)), folder)
        L7_outfile = infile_base.replace('L5', 'L7').replace('.i3.zst', '')
        L7_logname = L7_outfile

        L8_outfile = os.path.join(os.path.join(os.path.join(os.path.join(outdir_base, 'L8'), 'domeff{}'.format(domeff)), folder), infile_base.replace('L5', 'L8'))
        L8_logname = infile_base.replace('L5', 'L8').replace('.i3.zst', '')

        processing_dict['L6']['infiles'].append(infile)
        processing_dict['L6']['outfiles'].append(L6_outfile)
        processing_dict['L6']['lognames'].append(L6_logname)
        processing_dict['L6']['outfiles_hdf5'].append(L6_outfile.replace('.i3.zst', '.hdf5'))

        if osg:
            processing_dict['L7']['outdirs'].append( 'gsiftp://gridftp-users.icecube.wisc.edu{}'.format(L7_outdir) )
        else:
            processing_dict['L7']['outdirs'].append(L7_outdir)

        processing_dict['L7']['outfiles'].append(L7_outfile)
        processing_dict['L7']['lognames'].append(L7_logname)

        processing_dict['L8']['outfiles'].append(L8_outfile)
        processing_dict['L8']['lognames'].append(L8_logname)
        processing_dict['L8']['outfiles_hdf5'].append(L8_outfile.replace('.i3.zst', '.hdf5'))

        if not realrun:break

    def print_L6():
        if osg:
            subtarget.write('JOB\tjob'+str(jobName)+'\t' + processing_dict['L6']['submit_location_osg'])
        else:
            subtarget.write('JOB\tjob'+str(jobName)+'\t' + processing_dict['L6']['submit_location'])
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(jobName)+'\tGCD="' + processing_dict['GCD'] + '"')
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(jobName)+'\tINFILE="' + processing_dict['L6']['infiles'][i] + '"')
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(jobName)+'\tOUTFILE="' + processing_dict['L6']['outfiles'][i] +'"')
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(jobName)+'\tLOGNAME="' + processing_dict['L6']['lognames'][i] +'"')
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(jobName)+'\tidentifier_out="' + identifier_out +'"')
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(jobName)+'\tHDF5FILE="' + processing_dict['L6']['outfiles_hdf5'][i] +'"')
        subtarget.write("\n")
        subtarget.write('RETRY\tjob'+str(jobName)+'\t'+retries)
        subtarget.write("\n")

    def print_L7():
        if osg:
            subtarget.write('JOB\tjob'+str(jobName)+'\t' + processing_dict['L7']['submit_location_osg'])
            subtarget.write("\n")
            subtarget.write('VARS\tjob'+str(jobName)+'\tINFILE="' + 'gsiftp://gridftp-users.icecube.wisc.edu{}'.format(processing_dict['L6']['outfiles'][i]) + '"')
        else:
            subtarget.write('JOB\tjob'+str(jobName)+'\t' + processing_dict['L7']['submit_location'])
            subtarget.write("\n")
            subtarget.write('VARS\tjob'+str(jobName)+'\tINFILE="' + processing_dict['L6']['outfiles'][i] + '"')
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(jobName)+'\tOUTDIR="' + processing_dict['L7']['outdirs'][i] +'/"')
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(jobName)+'\tOUTFILE="' + processing_dict['L7']['outfiles'][i] +'"')
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(jobName)+'\tLOGNAME="' + processing_dict['L7']['lognames'][i] +'"')
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(jobName)+'\tidentifier_out="' + identifier_out +'"')
        subtarget.write("\n")
        subtarget.write('RETRY\tjob'+str(jobName)+'\t'+retries)
        subtarget.write("\n")

    def print_L8():
        if osg:
            subtarget.write('JOB\tjob'+str(jobName)+'\t' + processing_dict['L8']['submit_location_osg'])
        else:
            subtarget.write('JOB\tjob'+str(jobName)+'\t' + processing_dict['L8']['submit_location'])
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(jobName)+'\tGCD="' + processing_dict['GCD'] + '"')
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(jobName)+'\tINFILE="' + os.path.join(processing_dict['L7']['outdirs'][i].split('gsiftp://gridftp-users.icecube.wisc.edu')[-1], processing_dict['L7']['outfiles'][i]) + '.i3.zst"')
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(jobName)+'\tOUTFILE="' + processing_dict['L8']['outfiles'][i] +'"')
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(jobName)+'\tLOGNAME="' + processing_dict['L8']['lognames'][i] +'"')
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(jobName)+'\tidentifier_out="' + identifier_out +'"')
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(jobName)+'\tHDF5FILE="' + processing_dict['L8']['outfiles_hdf5'][i] +'"')
        subtarget.write("\n")
        subtarget.write('RETRY\tjob'+str(jobName)+'\t'+retries)
        subtarget.write("\n")

    job_counter = 0

    i = 0
    while job_counter < int(limit) and job_counter < len(processing_dict['L8']['outfiles']) and i < len(processing_dict['L8']['outfiles']):
        processing_level = ""
        if not os.path.isfile(processing_dict['L8']['outfiles'][i]) or os.path.getsize(processing_dict['L8']['outfiles'][i]) < 200000:
            processing_level = 'L8'
            if not os.path.isfile(processing_dict['L7']['outfiles'][i]) or os.path.getsize(processing_dict['L7']['outfiles'][i]) < 200000:
                processing_level = 'L7'
                if not os.path.isfile(processing_dict['L6']['outfiles'][i]):
                    processing_level = 'L6'
        # else:
        #     print('File already processed: '+ processing_dict['L8']['outfiles'][i])
        # print(job_counter)
        
        if processing_level == 'L8':	
            jobName = str(job_counter+1) + 'c'
            print_L8()
            job_counter +=1
        elif processing_level == 'L7':	
            jobName = str(job_counter+1) + 'b'
            print_L7()
            jobName = str(job_counter+1) + 'c'
            print_L8()
            subtarget.write('PARENT\tjob'+str(job_counter+1)+'b'+'\tCHILD '+'job'+str(job_counter+1)+'c')
            subtarget.write("\n")
            job_counter +=1
        elif processing_level == 'L6':	
            jobName = str(job_counter+1) + 'a'
            print_L6()
            jobName = str(job_counter+1) + 'b'
            print_L7()
            jobName = str(job_counter+1) + 'c'
            print_L8()
            subtarget.write('PARENT\tjob'+str(job_counter+1)+'a'+'\tCHILD '+'job'+str(job_counter+1)+'b')
            subtarget.write("\n")
            subtarget.write('PARENT\tjob'+str(job_counter+1)+'b'+'\tCHILD '+'job'+str(job_counter+1)+'c')
            subtarget.write("\n")
            job_counter +=1
        i += 1

    subtarget.close()
    ###### End ######

    print('Jobs setup: {}'.format(job_counter))
    if job_counter == 0:
        os.remove(subf)
        os.remove(L6_L7_L8_data_file)
        print('No jobs, removing dag/info file')
    else:
        if osg:
            subf_on_sub1 = os.path.join('/scratch/{}/dagman'.format(username), subf_base)
            print("Copying dag/sub/process scripts to sub-1.")
            os.system('scp {} sub-1:{}'.format(subf, subf_on_sub1))

            os.system('scp {} sub-1:{}'.format(os.path.join(submit_location, 'submit-L6-OSG.condor'), processing_dict['L6']['submit_location_osg']))
            os.system('scp {} sub-1:{}'.format(os.path.join(submit_location, 'submit-L7-OSG.condor'), processing_dict['L7']['submit_location_osg']))
            os.system('scp {} sub-1:{}'.format(os.path.join(submit_location, 'submit-L8-OSG.condor'), processing_dict['L8']['submit_location_osg']))

            os.system('scp {} sub-1:{}'.format(os.path.join(process_location, 'process_flercnn_L6.sh'), os.path.join(submit_location_osg, 'process_flercnn_L6-OSG.sh')))
            os.system('scp {} sub-1:{}'.format(os.path.join(process_location, 'process_flercnn_L7-OGS.sh'), os.path.join(submit_location_osg, 'process_flercnn_L7-OSG.sh')))
            os.system('scp {} sub-1:{}'.format(os.path.join(process_location, 'process_flercnn_L8.sh'), os.path.join(submit_location_osg, 'process_flercnn_L8-OSG.sh')))

            print('Submit (from sub-1!) with: "condor_submit_dag {}"'.format(subf_on_sub1))
        else:
            print('Dag file: {}'.format(subf))
            print('Submit with:\ncondor_submit_dag {}'.format(subf))
    print('done...')


if __name__ == '__main__':
    main()

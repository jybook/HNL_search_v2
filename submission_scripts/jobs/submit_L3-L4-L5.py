#!/bin/sh /cvmfs/icecube.opensciencegrid.org/users/BeyondStandardModel/bsm-py2-v3.1.1-1/icetray-start-standard
#METAPROJECT /data/user/lfischer/software/oscnext_hnl/build/

#=================================
# Spencer N. Axani
# Questions? saxani@mit.edu

# Modified by L. Fischer
#=================================

import glob, os, sys
from optparse import OptionParser

def main():
    '''
    Main function: Write dag file to submit L3/L4/L5 level jobs to condor.
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
    parser.add_option("-g", "--gcdfile",
                    default='/cvmfs/icecube.opensciencegrid.org/data/GCD/GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE_withScaledNoise.i3.gz',
                    help="Read in GCD file [default: %default]")
    parser.add_option("-d","--datatype",type="string", default="LeptonInjector",
                    dest="datatype", help="Generator data type [default: %default]")
    parser.add_option("-r","--REALRUN", default=False, action="store_true",
                    dest="REALRUN", help="Do real run. (otherwise just test file) [default: %default]")
    parser.add_option("--mc_location",type="string", default='/data/ana/BSM/HNL/MC/',
            dest="mc_location", help="The location to write the MC files to [default: %default]")
    parser.add_option("--dag_location",type="string", default='/scratch/lfischer/dagman/',
            dest="dag_location", help="The location to write the dag files to [default: %default]")
    (options,args) = parser.parse_args()

    if not options.identifier:
        options.identifier_out = parser.error('Set identifier not given.')
    if not options.identifier_out:
        options.identifier_out = options.identifier

    domeff          = options.domeff
    limit 		    = options.limit
    retries 	    = '2'
    gcdfile         = options.gcdfile
    datatype        = options.datatype
    identifier 	    = options.identifier
    identifier_out	= options.identifier_out
    realrun 		= options.REALRUN
    mc_location 	= options.mc_location
    dag_location	= options.dag_location

    catchup = 'False'

    list_domeff = ['0.95', '0.97', '1.00', '1.03', '1.05', '1.10', '0.90', '2.00']

    print("Using dom eficiency: {}".format(domeff))
    if domeff in list_domeff:
        D = str(list_domeff.index(domeff))
    else:
        print('Bad dom eff')
        sys.exit()

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

    ###### Write L3/L4/L5 information file ######
    counter = 0
    L3_L4_L5_data_file = outdir_base+"/simulation_log_L3_L4_L5_{}_v{}.txt".format(identifier, counter)
    while(True):
        if os.path.isfile(L3_L4_L5_data_file):
            counter +=1
            L3_L4_L5_data_file = L3_L4_L5_data_file.split('_v')[0] + '_v{}.txt'.format(counter)
        else:
            break

    text_file = open(L3_L4_L5_data_file, "w")
    text_file.write("Identifier: " + str(identifier) + '\n')
    text_file.write("Identifier Out: " + str(identifier_out) + '\n')
    text_file.write("DOM efficiency: " + str(domeff) + '\n')
    text_file.write("Limit: " + str(limit) + '\n')
    text_file.write("GCD file: " + str(gcdfile) + '\n')
    text_file.write('\n')
    text_file.close()


    ###### Write dag file ######
    counter = 0
    subf = os.path.join(dag_location, 'dagman-L3_L4_L5_{}_v{}.dag'.format(identifier_out, counter))

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
    processing_dict['domeff'] = domeff
    processing_dict['datatype'] = datatype

    processing_dict['L3'] = {}
    processing_dict['L3']['infiles'] = []
    processing_dict['L3']['outfiles'] = []
    processing_dict['L3']['outfiles_hdf5'] = []
    processing_dict['L3']['submit_location'] = submit_location + '/submit-L3.condor'

    processing_dict['L4'] = {}
    processing_dict['L4']['outfiles'] = []
    processing_dict['L4']['outfiles_hdf5'] = []
    processing_dict['L4']['submit_location'] = submit_location + '/submit-L4.condor'

    processing_dict['L5'] = {}
    processing_dict['L5']['outfiles'] = []
    processing_dict['L5']['outfiles_hdf5'] = []
    processing_dict['L5']['submit_location'] = submit_location + '/submit-L5.condor'


    infile_list = []

    folders = sorted(glob.glob(mc_location+identifier.split('/')[0]+'/L2/domeff'+domeff+'/[0-9]*'))
    # print(folders)

    level_list = ['L3','L4','L5']

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

    print('Number of L2 files: '+str(len(infile_list)))

    if not limit:
        limit = len(infile_list)

    for infile in infile_list:
        # print(infile)
        processing_dict['L3']['infiles'].append(infile)
        infile_base = os.path.basename(infile)
        folder = infile.split('/')[-2]
        # print(folder)

        L3_outfile = os.path.join(os.path.join(os.path.join(os.path.join(outdir_base, 'L3'), 'domeff{}'.format(domeff)), folder), infile_base.replace('L2', 'L3'))
        L4_outfile = os.path.join(os.path.join(os.path.join(os.path.join(outdir_base, 'L4'), 'domeff{}'.format(domeff)), folder), infile_base.replace('L2', 'L4'))
        L5_outfile = os.path.join(os.path.join(os.path.join(os.path.join(outdir_base, 'L5'), 'domeff{}'.format(domeff)), folder), infile_base.replace('L2', 'L5'))

        processing_dict['L3']['outfiles'].append(L3_outfile)
        processing_dict['L3']['outfiles_hdf5'].append(L3_outfile.replace('.i3.zst', '.hdf5'))
        processing_dict['L4']['outfiles'].append(L4_outfile)
        processing_dict['L4']['outfiles_hdf5'].append(L4_outfile.replace('.i3.zst', '.hdf5'))
        processing_dict['L5']['outfiles'].append(L5_outfile)
        processing_dict['L5']['outfiles_hdf5'].append(L5_outfile.replace('.i3.zst', '.hdf5'))

        if not realrun:break

    def print_L3():
        subtarget.write('JOB\tjob'+str(jobName)+'\t' + processing_dict['L3']['submit_location'])
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(jobName)+'\tdatatype="' + processing_dict['datatype'] + '"')
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(jobName)+'\tgcdfile="' + processing_dict['GCD'] + '"')
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(jobName)+'\tinfile="' + processing_dict['L3']['infiles'][i] + '"')
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(jobName)+'\toutfile="' + processing_dict['L3']['outfiles'][i] +'"')
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(jobName)+'\thdf5file="' + processing_dict['L3']['outfiles_hdf5'][i] +'"')
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(jobName)+'\tidentifier_out="'+identifier_out+'"')
        subtarget.write("\n")
        subtarget.write('RETRY\tjob'+str(jobName)+'\t'+retries)
        subtarget.write("\n")

    def print_L4():
        subtarget.write('JOB\tjob'+str(jobName)+'\t' + processing_dict['L4']['submit_location'])
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(jobName)+'\tdatatype="' + processing_dict['datatype'] + '"')
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(jobName)+'\tgcdfile="' + processing_dict['GCD'] + '"')
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(jobName)+'\tinfile="' + processing_dict['L3']['outfiles'][i] + '"')
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(jobName)+'\toutfile="' + processing_dict['L4']['outfiles'][i] +'"')
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(jobName)+'\thdf5file="' + processing_dict['L4']['outfiles_hdf5'][i] +'"')
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(jobName)+'\tidentifier_out="'+identifier_out+'"')
        subtarget.write("\n")
        subtarget.write('RETRY\tjob'+str(jobName)+'\t'+retries)
        subtarget.write("\n")

    def print_L5():
        subtarget.write('JOB\tjob'+str(jobName)+'\t' + processing_dict['L5']['submit_location'])
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(jobName)+'\tdatatype="' + processing_dict['datatype'] + '"')
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(jobName)+'\tgcdfile="' + processing_dict['GCD'] + '"')
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(jobName)+'\tinfile="' + processing_dict['L4']['outfiles'][i] + '"')
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(jobName)+'\toutfile="' + processing_dict['L5']['outfiles'][i] +'"')
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(jobName)+'\thdf5file="' + processing_dict['L5']['outfiles_hdf5'][i] +'"')
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(jobName)+'\tidentifier_out="'+identifier_out+'"')
        subtarget.write("\n")
        subtarget.write('RETRY\tjob'+str(jobName)+'\t'+retries)
        subtarget.write("\n")

    job_counter = 0

    i = 0
    while job_counter < int(limit) and job_counter < len(processing_dict['L5']['outfiles']) and i < len(processing_dict['L5']['outfiles']):
        processing_level = ""
        if not os.path.isfile(processing_dict['L5']['outfiles'][i]) or os.path.getsize(processing_dict['L5']['outfiles'][i]) < 200000:
            processing_level = 'L5'
            if not os.path.isfile(processing_dict['L4']['outfiles'][i]) or os.path.getsize(processing_dict['L4']['outfiles'][i]) < 200000:
                processing_level = 'L4'
                if catchup == 'True':
                    if not os.path.isfile(processing_dict['L3']['outfiles'][i])  or os.path.getsize(processing_dict['L3']['outfiles'][i]) < 200000:
                        i += 1
                        continue
                    else:
                        print('year')
                else:
                    if not os.path.isfile(processing_dict['L3']['outfiles'][i]):
                        processing_level = 'L3'
        # else:
        #     print('File already processed: '+ processing_dict['L5']['outfiles'][i])
        # print(job_counter)
        
        if processing_level == 'L5':	
            jobName = str(job_counter+1) + 'c'
            print_L5()
            job_counter +=1
        elif processing_level == 'L4':	
            jobName = str(job_counter+1) + 'b'
            print_L4()
            jobName = str(job_counter+1) + 'c'
            print_L5()
            subtarget.write('PARENT\tjob'+str(job_counter+1)+'b'+'\tCHILD '+'job'+str(job_counter+1)+'c')
            subtarget.write("\n")
            job_counter +=1
        elif processing_level == 'L3':	
            jobName = str(job_counter+1) + 'a'
            print_L3()
            jobName = str(job_counter+1) + 'b'
            print_L4()
            jobName = str(job_counter+1) + 'c'
            print_L5()
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
        os.remove(L3_L4_L5_data_file)
        print('No jobs, removing dag/info file')
    else:
        print('Dag file: {}'.format(subf))
        print('Submit with:\ncondor_submit_dag {}'.format(subf))
    print('done...')


if __name__ == '__main__':
    main()

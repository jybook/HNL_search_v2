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
    Main function: Write dag file to submit Det/L1/L2 level jobs to condor.
    '''

    usage = "usage: %prog [options] inputfile"
    parser = OptionParser(usage)
    parser.add_option("-i", "--identifier",type="string",
                    dest="identifier", help="Set name (infiles) [default: %default]")
    parser.add_option("-o", "--identifier_out",type="string",
                    dest="identifier_out", help="Set name (outfiles), will use infile identifier if not set [default: %default]")
    parser.add_option("-e", "--domeff",type="string", default='1.00',
                    dest="domeff", help="Which dom efficiency are we using.ie. '0.95','0.97','1.00','1.03','1.05','1.10','0.90', '2.00' [default: %default]")
    parser.add_option("-l", "--limit",type="int",
                    dest="limit", help="How many files do you want to have processed (absolute number?), will use all if unset [default: %default]")
    parser.add_option("-m","--icemodel", default="spice_3.2.1",
                    dest="ICEMODEL",help="Ice model (spice_3.2.1, spice_mie, spice_lea, etc) [default: %default]")
    parser.add_option("-n","--noise",default="vuvuzela",
                    dest="NOISE",help="Noise model (vuvuzela/poisson/none) [default: %default]")
    parser.add_option("-g", "--gcdfile",
                    default='/cvmfs/icecube.opensciencegrid.org/data/GCD/GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE_withScaledNoise.i3.gz',
                    help="Read in GCD file [default: %default]")
    parser.add_option("--holeice",type="string", default="angsens/as.flasher_p1_0.30_p2_-1",
                    dest="holeice", help="What hole ice model do you want? [default: %default]")
    parser.add_option("--REALRUN", default=False, action="store_true",
                    dest="REALRUN", help="Do real run. (otherwise just test file) [default: %default]")
    parser.add_option("--mc_location",type="string", default='/data/ana/BSM/HNL/MC/',
            dest="mc_location", help="The location to write the MC files to [default: %default]")
    parser.add_option("--dag_location",type="string", default='/scratch/lfischer/dagman/',
            dest="dag_location", help="The location to write the dag files to [default: %default]")
    (options,args) = parser.parse_args()

    if not options.identifier_out:
        options.identifier_out = options.identifier

    domeff			= options.domeff
    holeice 	    = options.holeice
    limit 		    = options.limit
    retries 	    = '2'
    gcdfile			= options.gcdfile
    identifier 	    = options.identifier
    identifier_out	= options.identifier_out
    realrun 		= options.REALRUN
    icemodel		= options.ICEMODEL
    noise           = options.NOISE
    mc_location 	= options.mc_location
    dag_location	= options.dag_location

    catchup = 'False'

    print('Using hole ice: {}'.format(holeice))

    list_domeff = ['0.90', '0.95', '0.97', '1.00', '1.03', '1.05', '1.10', '2.00']

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

    ###### Write Det/L1/L2 information file ######
    counter = 0
    det_l1_l2_data_file = outdir_base+"/simulation_log_Det_L1_L2_{}_v{}.txt".format(identifier, counter)
    while(True):
        if os.path.isfile(det_l1_l2_data_file):
            counter +=1
            det_l1_l2_data_file = det_l1_l2_data_file.split('_v')[0] + '_v{}.txt'.format(counter)
        else:
            break

    text_file = open(det_l1_l2_data_file, "w")
    text_file.write("Identifier: " + str(identifier) + '\n')
    text_file.write("Identifier Out: " + str(identifier_out) + '\n')
    text_file.write("DOM efficiency: " + str(domeff) + '\n')
    text_file.write("Noise model: " + str(noise) + '\n')
    text_file.write("Limit: " + str(limit) + '\n')
    text_file.write("Hole Ice: " + str(holeice) + '\n')
    text_file.write("GCD file: " + str(gcdfile) + '\n')
    text_file.write('\n')
    text_file.close()
    if not realrun:os.remove(det_l1_l2_data_file)


    ###### Write dag file ######
    counter = 0
    subf = os.path.join(dag_location, 'dagman-Det_L1_L2_{}_v{}.dag'.format(identifier_out, counter))

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
    processing_dict['holeice'] = holeice
    processing_dict['icemodel'] = icemodel
    processing_dict['noise'] = noise
    processing_dict['Det'] = {}
    processing_dict['Det']['outfiles'] = []
    processing_dict['Det']['infiles'] = []
    processing_dict['Det']['submit_location'] = submit_location + '/submit-Det.condor'

    processing_dict['L1'] = {}
    processing_dict['L1']['outfiles'] = []
    processing_dict['L1']['submit_location'] = submit_location + '/submit-L1.condor'

    processing_dict['L2'] = {}
    processing_dict['L2']['outfiles'] = []
    processing_dict['L2']['submit_location'] = submit_location + '/submit-L2.condor'

    infile_list = []

    folders = sorted(glob.glob(mc_location+identifier.split('/')[0]+'/Phot/[0-9]*'))

    level_list = ['Det', 'L1', 'L2']

    for folder in folders:
        folder_name = os.path.basename(folder)
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

    print('Number of Photon files: '+str(len(infile_list)))

    if not limit:
        limit = len(infile_list)

    for infile in infile_list:
        processing_dict['Det']['infiles'].append(infile)
        infile_base = os.path.basename(infile)
        folder = infile.split('/')[-2]

        det_outfile 	= os.path.join(os.path.join(os.path.join(os.path.join(outdir_base, 'Det'), 'domeff{}'.format(domeff)), folder), infile_base.replace('Phot', 'Det'))
        L1_outfile 		= det_outfile.replace('Det','L1')
        L2_outfile 		= det_outfile.replace('Det','L2')

        processing_dict['Det']['outfiles'].append(det_outfile)
        processing_dict['L1']['outfiles'].append(L1_outfile)
        processing_dict['L2']['outfiles'].append(L2_outfile)

        if not realrun:break

    def print_det():
        subtarget.write('JOB\tjob'+str(jobName)+'\t' + processing_dict['Det']['submit_location'])
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(jobName)+'\tgcdfile="' + processing_dict['GCD'] + '"')
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(jobName)+'\tinfile="' + processing_dict['Det']['infiles'][i] + '"')
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(jobName)+'\toutfile="' + processing_dict['Det']['outfiles'][i] +'"')
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(jobName)+'\tdomeff="'+processing_dict['domeff']+'"')
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(jobName)+'\tholeice="'+processing_dict['holeice']+'"')
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(jobName)+'\ticemodel="'+processing_dict['icemodel']+'"')
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(jobName)+'\tnoise="'+processing_dict['noise']+'"')
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(jobName)+'\tidentifier_out="'+identifier_out+'"')
        subtarget.write("\n")
        subtarget.write('RETRY\tjob'+str(jobName)+'\t'+retries)
        subtarget.write("\n")

    def print_L1():
        subtarget.write('JOB\tjob'+str(jobName)+'\t' + processing_dict['L1']['submit_location'])
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(jobName)+'\tgcdfile="' + processing_dict['GCD'] + '"')
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(jobName)+'\tinfile="' + processing_dict['Det']['outfiles'][i] + '"')
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(jobName)+'\toutfile="' + processing_dict['L1']['outfiles'][i] +'"')
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(jobName)+'\tidentifier_out="'+identifier_out+'"')
        subtarget.write("\n")
        subtarget.write('RETRY\tjob'+str(jobName)+'\t'+retries)
        subtarget.write("\n")

    def print_L2():
        subtarget.write('JOB\tjob'+str(jobName)+'\t' + processing_dict['L2']['submit_location'])
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(jobName)+'\tgcdfile="' + processing_dict['GCD'] + '"')
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(jobName)+'\tinfile="' + processing_dict['L1']['outfiles'][i] + '"')
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(jobName)+'\toutfile="' + processing_dict['L2']['outfiles'][i] +'"')
        subtarget.write("\n")
        subtarget.write('VARS\tjob'+str(jobName)+'\tidentifier_out="'+identifier_out+'"')
        subtarget.write("\n")
        subtarget.write('RETRY\tjob'+str(jobName)+'\t'+retries)
        subtarget.write("\n")

    jobNum = 1
    job_counter = 0
    job_num_list = []

    i = 0
    while i < int(limit) and job_counter < len(processing_dict['L2']['outfiles']) and i < len(processing_dict['L2']['outfiles']):
        processing_level = ""
        if not os.path.isfile(processing_dict['L2']['outfiles'][i]) or os.path.getsize(processing_dict['L2']['outfiles'][i]) < 200000:
            processing_level = 'L2'
            if not os.path.isfile(processing_dict['L1']['outfiles'][i]) or os.path.getsize(processing_dict['L1']['outfiles'][i]) < 200000:
                processing_level = 'L1'
                if catchup == 'True':
                    if not os.path.isfile(processing_dict['Det']['outfiles'][i])  or os.path.getsize(processing_dict['Det']['outfiles'][i]) < 200000:
                        i += 1
                        continue
                    else:
                        print('year')
                else:
                    if not os.path.isfile(processing_dict['Det']['outfiles'][i]):
                        processing_level = 'Det'
        # else:
        #     print('File already processed: '+ processing_dict['L2']['outfiles'][i])
        # print(job_counter)
        
        if processing_level == 'L2':
            job_num_list.append(jobNum)
            jobName = str(jobNum) + 'c'
            print_L2()
            job_counter +=1
            jobNum +=1
        elif processing_level == 'L1':	
            job_num_list.append(jobNum)
            jobName = str(jobNum) + 'b'
            print_L1()
            jobName = str(jobNum) + 'c'
            print_L2()
            subtarget.write('PARENT\tjob'+str(jobNum)+'b'+'\tCHILD '+'job'+str(jobNum)+'c')
            subtarget.write("\n")
            job_counter +=1
            jobNum +=1
        elif processing_level == 'Det':	
            job_num_list.append(jobNum)
            jobName = str(jobNum) + 'a'
            print_det()
            subtarget.write('VARS\tjob'+str(jobName)+'\tfilenum="' + str(jobNum) + '"')
            subtarget.write("\n")
            jobName = str(jobNum) + 'b'
            print_L1()
            jobName = str(jobNum) + 'c'
            print_L2()
            subtarget.write('PARENT\tjob'+str(jobNum)+'a'+'\tCHILD '+'job'+str(jobNum)+'b')
            subtarget.write("\n")
            subtarget.write('PARENT\tjob'+str(jobNum)+'b'+'\tCHILD '+'job'+str(jobNum)+'c')
            subtarget.write("\n")
            job_counter +=1
            jobNum +=1
        i += 1

    subtarget.close()
    ###### End ######

    print('Jobs setup: {}'.format(job_counter))
    if job_counter == 0:
        os.remove(subf)
        os.remove(det_l1_l2_data_file)
        print('No jobs, removing dag/info file')
    else:
        print('Dag file: {}'.format(subf))
        print('Submit with:\ncondor_submit_dag {}'.format(subf))
    print('done...')


if __name__ == '__main__':
    main()

#!/bin/sh /cvmfs/icecube.opensciencegrid.org/users/BeyondStandardModel/bsm-py2-v3.1.1-1/icetray-start-standard
#METAPROJECT /cvmfs/icecube.opensciencegrid.org/users/lfischer/oscnext_hnl_meta/build/

#=================================
# Spencer N. Axani
# Questions? saxani@mit.edu

# Modified by D. Vannerom/L. Fischer
#=================================

import glob, os

from optparse import OptionParser

def main():
    '''
    Main function: Write dag file to submit Photon level jobs to condor.
    '''

    usage = "usage: %prog [options] inputfile"
    parser = OptionParser(usage)
    parser.add_option("-i", "--identifier",type="string",
                    dest="identifier", help="Set name (infiles) [default: %default]")
    parser.add_option("-o", "--identifier_out",type="string",
                    dest="identifier_out", help="Set name (outfiles), will use infile identifier if not set [default: %default]")
    parser.add_option("--ram",type="string",default='4000',help='Ram for each cluster job [default: %default]')
    parser.add_option("-l", "--limit",type="string",
                    dest="limit", help="Set to 1-10000 to setup files from 1-10000") 
    parser.add_option("-r", "--runnumber", type="string", default="1",
                    help="The run/dataset number for this simulation, is used as seed for random generator [default: %default]")
    parser.add_option("-m","--icemodel", default="spice_3.2.1",
                    dest="ICEMODEL",help="Ice model (spice_3.2.1, spice_mie, spice_lea, etc) [default: %default]")
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
    parser.add_option("--username",type="string",
                    dest="username", help="Set username (needed for OSG support) [default: %default]")
    parser.add_option("--osg", default=False, action="store_true",
                     dest="osg", help="Is this job running on the OSG? [default: %default]")
    (options,args) = parser.parse_args()

    if not options.username:
        parser.error('Username must be specified with --username option.')
    if not options.identifier_out:
        options.identifier_out = options.identifier
    holeice 	    = options.holeice

    identifier		= options.identifier
    identifier_out	= options.identifier_out
    limit 			= options.limit
    ram 			= options.ram
    gcdfile 		= options.gcdfile
    RUNNUM 			= options.runnumber
    icemodel		= options.ICEMODEL
    realrun 		= options.REALRUN
    mc_location 	= options.mc_location
    dag_location	= options.dag_location
    osg             = options.osg
    username        = options.username

    print('Using hole ice: {}'.format(holeice))

    mother_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    submit_location = os.path.join(mother_path, 'submit')
    process_location = os.path.join(mother_path, 'process')

    if not os.path.exists(dag_location):
        os.makedirs(dag_location)

    outdir_base = mc_location + identifier_out+'/'+'Phot'
    if not os.path.exists(outdir_base):
        os.makedirs(outdir_base)
        os.chmod(outdir_base, 0o775)


    ###### Write Phot information file ######
    counter = 0
    phot_data_file = outdir_base+"/Phot_data_{}_v{}.txt".format(identifier, counter)
    while(True):
        if os.path.isfile(phot_data_file):
            counter +=1
            phot_data_file = phot_data_file.split('_v')[0] + '_v{}.txt'.format(counter)
        else:
            break

    text_file = open(phot_data_file, "w")
    text_file.write("{0:<25} {1} \n".format('Identifier In:', identifier))
    text_file.write("{0:<25} {1} \n".format('Identifier Out:', identifier_out))
    text_file.write("{0:<25} {1} \n".format('Limit:', limit))
    text_file.write("{0:<25} {1} \n".format('IceModel:', icemodel))
    text_file.write("{0:<25} {1} \n".format('HoleIce:', holeice))
    text_file.write("{0:<25} {1} \n".format('Using OSG?', osg))
    text_file.write("{0:<25} {1} \n".format('Using RAM:', ram))
    text_file.close()

    ###### Write dag file ######
    counter = 0
    subf_base = 'dagman-Phot_{}_v{}.dag'.format(identifier_out, counter)
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


    job_number = 1

    infolders = sorted(glob.glob(mc_location + identifier+ '/Gen/[0-9]*'))

    for folder in infolders:
        this_outdir = os.path.join(outdir_base, os.path.basename(folder))
        if not os.path.exists(this_outdir):
            os.makedirs(this_outdir)
            os.chmod(this_outdir, 0o775)
        infiles = sorted(glob.glob(folder+'/*.i3.zst'))
        if len(infiles) == 0:
            infiles = sorted(glob.glob(folder+'/*.i3.bz2'))

        for infile in infiles:
            file_number = int(infile.split('Gen_')[-1].split('.i3')[0])
            # print(file_number)
            if limit:
                if job_number < int(limit.split('-')[0]) or job_number > int(limit.split('-')[-1]):
                    continue
            outfile = os.path.join(this_outdir, (os.path.basename(infile)).replace('Gen','Phot'))
            # print(outfile)
            if (not os.path.exists(outfile) or os.path.getsize(outfile) < 4000):
                if osg:
                    subtarget.write('JOB\tjob'+str(job_number)+'\t'+ '/scratch/{}/submission/submit-Phot-OSG.condor'.format(username))
                else:
                    subtarget.write('JOB\tjob'+str(job_number)+'\t'+ submit_location + '/submit-Phot.condor')
                subtarget.write("\n")
                subtarget.write('VARS\tjob'+str(job_number)+'\tgcdfile="' + str(gcdfile) + '"')
                subtarget.write("\n")
                subtarget.write('VARS\tjob'+str(job_number)+'\tinfile="' + str(infile) + '"')
                subtarget.write("\n")
                subtarget.write('VARS\tjob'+str(job_number)+'\toutfile="' + str(outfile) + '"')
                subtarget.write("\n")
                subtarget.write('VARS\tjob'+str(job_number)+'\tRUNNUM="' + str(RUNNUM) + '"')
                subtarget.write("\n")
                subtarget.write('VARS\tjob'+str(job_number)+'\tfile_number="' + str(file_number) + '"')
                subtarget.write("\n")
                subtarget.write('VARS\tjob'+str(job_number)+'\tram="' + ram +'"')
                subtarget.write("\n")
                subtarget.write('VARS\tjob'+str(job_number)+'\tidentifier_out="' + identifier_out +'"')
                subtarget.write("\n")
                subtarget.write('VARS\tjob'+str(job_number)+'\ticemodel="' + icemodel +'"')
                subtarget.write("\n")
                subtarget.write('VARS\tjob'+str(job_number)+'\tholeice="' + holeice +'"')
                subtarget.write("\n")
                subtarget.write('RETRY\tjob'+str(job_number)+'\t2')
                subtarget.write("\n")
                job_number += 1
            if not realrun and job_number == 2:break
        if not realrun and job_number == 2:break
    subtarget.close()
    ###### End ######

    print('We setup: {} files'.format(job_number-1))
    if (job_number - 1) == 0:
        os.remove(subf)
        os.remove(phot_data_file)
        print('No jobs, removing dag/info file')
    else:
        if osg:
            subf_on_sub1 = os.path.join('/scratch/{}/dagman'.format(username), subf_base)
            print("Copying dag/sub/process scripts to sub-1.")
            os.system('scp {} sub-1:{}'.format(subf, subf_on_sub1))
            os.system('scp {} sub-1:{}'.format(submit_location + '/submit-Phot-OSG.condor', '/scratch/{}/submission/submit-Phot-OSG.condor'.format(username)))
            os.system('scp {} sub-1:{}'.format(process_location + '/process_Phot.py', '/scratch/{}/submission/process_Phot-OSG.py'.format(username)))
            print('Submit (from sub-1!) with: "condor_submit_dag {}"'.format(subf_on_sub1))
        else:
            print('Dag file: {}'.format(subf))
            print('Submit with: "condor_submit_dag {}"'.format(subf))
    print('done...')


if __name__ == '__main__':
    main()

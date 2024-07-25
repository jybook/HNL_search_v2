#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/icetray-start
#METAPROJECT /cvmfs/icecube.opensciencegrid.org/users/lfischer/oscnext_hnl_meta/build

'''
Author: Leander Fischer

Script to produce dagman file to process oscnext simulation with Taupede and Millipede.
'''

import glob, os, sys, random
import numpy as np
from shutil import copyfile
from optparse import OptionParser
from os.path import expandvars
usage = "usage: %prog [options] inputfile"
parser = OptionParser(usage)

parser.add_option("-l", "--limit",type="int",
                dest="limit", help="How many files do you want?.")
parser.add_option("-i", "--identifier",type="string",
			    dest="identifier", help="Set name (infiles).")
parser.add_option("-o", "--identifier_out",type="string",
			    dest="identifier_out", help="Set name (outfiles).")
parser.add_option("-g", "--gcdfile",
			    default='/cvmfs/icecube.opensciencegrid.org/data/GCD/GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE_withScaledNoise.i3.gz',
			    help="Read in GCD file")
parser.add_option("-r", "--REALRUN", default=False, action="store_true",
			    dest="REALRUN", help="Do real run. (otherwise just test file)")
(options,args) = parser.parse_args()

if not options.identifier_out:
	options.identifier_out = options.identifier

limit 		    = options.limit
retries 	    = '2'
gcdfile         = options.gcdfile
identifier 	    = options.identifier
identifier_out	= options.identifier_out
realrun 		= options.REALRUN


###### Paths to set manually ######
in_location = '/data/ana/LE/oscNext/pass2/'
dag_location_path = '/scratch/lfischer/dagman/'
output_base = '/data/ana/BSM/HNL/MC/oscNext'
###### End ######

mc_sets = {
    'genie':['120000', '120001', '120004', '120100', '120102', '120500', '120503', '140000', '140001', '140004', '140100', '140102', '140500', '140503', '160000', '160001', '160004', '160100', '160102', '160500', '160503'],
    'muongun':['130000', '130001', '130004', '130100', '130102', '130500', '130503'],
    'noise':['888003']
}

mother_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
submit_location = os.path.join(mother_path, 'submit')

if not os.path.exists(dag_location_path):
	os.makedirs(dag_location_path)

dag_location = dag_location_path

if not os.path.exists(gcdfile):
	print('Missing GCD file: ' + gcdfile)
	sys.exit()
print('GCD file: '+gcdfile)

outdir_base = output_base
if not os.path.exists(outdir_base):
	os.makedirs(outdir_base)

try:	
	text_file = open(outdir_base+"/simulation_log_Taupede_Millipede.txt", "w")
	text_file.write("Identifier In: " + str(identifier) + '\n')
	text_file.write("Identifier Out: " + str(identifier_out) + '\n')
	text_file.write("Limit: " + str(limit) + '\n')
	text_file.write("GCD file: " + str(gcdfile) + '\n')
	text_file.write('\n')
	text_file.close()
except:
	pass

counter = 0
subf = os.path.join(dag_location, 'dagman-Taupede_Millipede_{}_v{}.dag'.format(identifier_out, counter))

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

processing_dict['Millipede'] = {}
processing_dict['Millipede']['outfiles'] = []
processing_dict['Millipede']['submit_location'] = submit_location + '/submit-Millipede.condor'


infile_list = []

mc_set_number = identifier
mc_set_name = None

for name, numbers in mc_sets.iteritems():
    if identifier in numbers:mc_set_name = name

if not mc_set_name:
	print('Set Number {} not known.'.format(identifier))
	sys.exit()

input_files_location = os.path.join(in_location, mc_set_name)
input_files_location = os.path.join(input_files_location, 'level7')
input_files_location = os.path.join(input_files_location, mc_set_number)
print('Input Files Location: {}'.format(input_files_location))

outdir_base = os.path.join(output_base, mc_set_name)
outdir_base = os.path.join(outdir_base, mc_set_number)
print('Output Files Location: {}'.format(outdir_base))


level_list = ['Taupede','Millipede']

folder = input_files_location

files = sorted(glob.glob(folder + '/*.i3.zst'))
if len(files) == 0:
    files = sorted(glob.glob(folder + '/*.i3.bz2'))
for file in files:
    infile_list.append(file)

for level in level_list:
    folder_out_dir = os.path.join(outdir_base, level)
    # print(folder_out_dir)
    if not os.path.exists(folder_out_dir):
        os.makedirs(folder_out_dir)
    os.chmod(folder_out_dir, 0o775)

print('Number of L7 files: '+str(len(infile_list)))

if not limit:
	limit = len(infile_list)

for infile in infile_list:
    # print(infile)
    processing_dict['Taupede']['infiles'].append(infile)
    infile_base = os.path.basename(infile)
    folder = infile.split('/')[-2]
	# print(folder)

    Taupede_outfile     = os.path.join(os.path.join(outdir_base, 'Taupede'), infile_base.replace('level7_v02.00', 'Taupede'))
    Millipede_outfile   = os.path.join(os.path.join(outdir_base, 'Millipede'), infile_base.replace('level7_v02.00', 'Millipede'))

    processing_dict['Taupede']['outfiles'].append(Taupede_outfile)
    processing_dict['Millipede']['outfiles'].append(Millipede_outfile)

    if not realrun:break

def print_Taupede():
    subtarget.write('JOB\tjob'+str(jobName)+'\t' + processing_dict['Taupede']['submit_location'])
    subtarget.write("\n")
    subtarget.write('VARS\tjob'+str(jobName)+'\tgcdfile="' + processing_dict['GCD'] + '"')
    subtarget.write("\n")
    subtarget.write('VARS\tjob'+str(jobName)+'\tinfile="' + processing_dict['Taupede']['infiles'][i] + '"')
    subtarget.write("\n")
    subtarget.write('VARS\tjob'+str(jobName)+'\toutfile="' + processing_dict['Taupede']['outfiles'][i] +'"')
    subtarget.write("\n")
    subtarget.write('VARS\tjob'+str(jobName)+'\tidentifier_out="'+identifier_out+'"')
    subtarget.write("\n")
    subtarget.write('RETRY\tjob'+str(jobName)+'\t'+retries)
    subtarget.write("\n")

def print_Millipede():
    subtarget.write('JOB\tjob'+str(jobName)+'\t' + processing_dict['Millipede']['submit_location'])
    subtarget.write("\n")
    subtarget.write('VARS\tjob'+str(jobName)+'\tgcdfile="' + processing_dict['GCD'] + '"')
    subtarget.write("\n")
    subtarget.write('VARS\tjob'+str(jobName)+'\tinfile="' + processing_dict['Taupede']['outfiles'][i] + '"')
    subtarget.write("\n")
    subtarget.write('VARS\tjob'+str(jobName)+'\toutfile="' + processing_dict['Millipede']['outfiles'][i] +'"')
    subtarget.write("\n")
    subtarget.write('VARS\tjob'+str(jobName)+'\tidentifier_out="'+identifier_out+'"')
    subtarget.write("\n")
    subtarget.write('RETRY\tjob'+str(jobName)+'\t'+retries)
    subtarget.write("\n")


jobNum = 1
job_counter = 0
job_num_list = []


i = 0
while job_counter < int(limit) and job_counter < len(processing_dict['Millipede']['outfiles']) and i < len(processing_dict['Millipede']['outfiles']):
    processing_level = ""
    if identifier == '888003':
        print('Not checking for file size to (noise set)')
        if not os.path.isfile(processing_dict['Millipede']['outfiles'][i]):
            processing_level = 'Millipede'
            if not os.path.isfile(processing_dict['Taupede']['outfiles'][i]):
                processing_level = 'Taupede'
        else:
            print('File already processed: '+ processing_dict['Millipede']['outfiles'][i])
    else:
        if not os.path.isfile(processing_dict['Millipede']['outfiles'][i]) or os.path.getsize(processing_dict['Millipede']['outfiles'][i]) < 200000:
            processing_level = 'Millipede'
            if not os.path.isfile(processing_dict['Taupede']['outfiles'][i]):
                processing_level = 'Taupede'
        else:
            print('File already processed: '+ processing_dict['Millipede']['outfiles'][i])
    # print(job_counter)
	
    if processing_level == 'Millipede':	
        job_num_list.append(jobNum)
        jobName = str(jobNum) + 'b'
        print_Millipede()
        job_counter +=1
        jobNum +=1
    elif processing_level == 'Taupede':	
        job_num_list.append(jobNum)
        jobName = str(jobNum) + 'a'
        print_Taupede()
        jobName = str(jobNum) + 'b'
        print_Millipede()
        subtarget.write('PARENT\tjob'+str(jobNum)+'a'+'\tCHILD '+'job'+str(jobNum)+'b')
        subtarget.write("\n")
        job_counter +=1
        jobNum +=1
    i += 1

subtarget.close()


print('Jobs setup: '+str(job_counter))
print('File: '+subf)
print('Submit with: "condor_submit_dag {}"'.format(subf))
print('done...')
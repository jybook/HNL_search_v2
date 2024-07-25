#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/icetray-start
#METAPROJECT /data/user/lfischer/software/oscnext_hnl/build/

#
# Author: Leander Fischer
#

import glob, os, sys, random
import numpy as np
from shutil import copyfile
from optparse import OptionParser
from os.path import expandvars

usage = "usage: %prog [options] inputfile"
parser = OptionParser(usage)
parser.add_option("-i", "--identifier",type="string",
		dest="identifier", help="HNL signal set number (19XXXX)")
parser.add_option("--nFiles",type="int", default="10000",
		dest="nFiles", help="The total number of files you want to generate [default: %default]")
parser.add_option("--nEvents",type="int", default="2000",
		dest="nEvents", help="The number of events in each file [default: %default]")
parser.add_option("--capacity",type="int", default="1000",
		dest="capacity", help="The number of files in each folder [default: %default]")
parser.add_option("--ram",type="string",
		dest="ram", help="ram for each process on the cluster [default: %default]", default="2000")
parser.add_option("--limit",type="string",
        dest="limit", help="set to 1-10000 to setup files from 1-10000") 

(options,args) = parser.parse_args()

identifier			= options.identifier
nFiles 				= options.nFiles
nEvents 			= options.nEvents
capacity 			= options.capacity
limit 				= options.limit
ram 				= options.ram

if limit == None:
    limit = '1-'+str(nFiles)

###### Paths to set manually ######
in_location = '/data/ana/BSM/HNL/MC/'
dag_location_path = '/scratch/lfischer/dagman/'
###### End ######

mother_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
submit_location = os.path.join(mother_path, 'submit')

dag_location = dag_location_path

outdir_base = in_location + identifier+'/Gen/'
if not os.path.exists(outdir_base):
	os.makedirs(outdir_base)
	os.chmod(outdir_base, 0o775)

###### Write gen information file ######
counter = 0
gen_data_file = outdir_base+"/Generation_data_{}_v{}.txt".format(identifier, counter)
while(True):
	if os.path.isfile(gen_data_file):
		counter +=1
		gen_data_file = gen_data_file.split('_v')[0] + '_v{}.txt'.format(counter)
	else:
		break

text_file = open(gen_data_file, "w")
text_file.write("{0:<25} {1} \n".format('identifier:', identifier))
text_file.write("{0:<25} {1} \n".format('nFiles:', nFiles))
text_file.write("{0:<25} {1} \n".format('nEvents', nEvents))
text_file.write("{0:<25} {1} \n".format('capacity:', capacity))
text_file.write("{0:<25} {1} \n".format('limit:', limit))
text_file.write("{0:<25} {1} \n".format('ram:', ram))
text_file.close()


###### Write dag file ######
counter = 0
subf = os.path.join(dag_location, 'dagman-Gen_{}_v{}.dag'.format(identifier, counter))

while(True):
	if os.path.isfile(subf):
		counter +=1
		subf = subf.split('_v')[0] + '_v{}.dag'.format(counter)
	else:
		break

subtarget = open(subf, 'w')
os.chmod(subf, 0o775)

smallest = limit.split('-')[0]
biggest = limit.split('-')[1]

eventlist = np.arange(int(smallest), int(biggest)+1)

job_number = 1

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
		continue

	subtarget.write('JOB\tjob'+str(job_number)+'\t'+ submit_location + '/submit-Gen_simple_double_cascades.condor')
	subtarget.write("\n")
	subtarget.write('VARS\tjob'+str(job_number)+'\toutfile="' + str(outfile) + '"')
	subtarget.write("\n")
	subtarget.write('VARS\tjob'+str(job_number)+'\tseed="' + str(seed) + '"')
	subtarget.write("\n")
	subtarget.write('VARS\tjob'+str(job_number)+'\tnEvents="' + str(nEvents) + '"')
	subtarget.write("\n")
	subtarget.write('VARS\tjob'+str(job_number)+'\tram="' + str(ram) + '"')
	subtarget.write("\n")
	subtarget.write('VARS\tjob'+str(job_number)+'\tidentifier_out="' + str(identifier) + '"')
	subtarget.write("\n")
	job_number += 1

subtarget.close()
###### End ######

import time
time.sleep(1)
print('Jobs setup: '+str(job_number-1)+' jobs')
if job_number -1 == 0:
    os.remove(subf)
    os.remove(gen_data_file)
    print('No jobs, removing dag/info file')
else:
    print('Dag file: '+subf)
    print('Submit with: "condor_submit_dag {}"'.format(subf))
print('done...')
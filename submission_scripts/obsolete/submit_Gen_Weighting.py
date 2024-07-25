#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/icetray-start
#METAPROJECT /cvmfs/icecube.opensciencegrid.org/users/lfischer/oscnext_hnl_meta/build

# Before submission you need to:
# source /cvmfs/icecube.opensciencegrid.org/users/BeyondStandardModel/bsm-py2-v3.1.1-1/setup.sh
# /cvmfs/icecube.opensciencegrid.org/users/lfischer/oscnext_meta_V01-00-05/build/env-shell.sh

import glob, os
from optparse import OptionParser

usage = "usage: %prog [options] inputfile"
parser = OptionParser(usage)
parser.add_option("-i", "--identifier",type="string",
				dest="identifier", help="Set name (infiles).")
parser.add_option("-o", "--identifier_out",type="string",
				dest="identifier_out", help="Set name (outfiles).")
parser.add_option("--ram",type="string")
parser.add_option("-l", "--limit",type="string",
				dest="limit", help="Set to 1-10000 to setup files from 1-10000") 
parser.add_option("-d","--datatype",type="string", default="LeptonInjector",
                  dest="datatype", help="Generator data type (used for weighting).")
parser.add_option("-g", "--gcdfile",
				default='/cvmfs/icecube.opensciencegrid.org/data/GCD/GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE_withScaledNoise.i3.gz',
				help="Read in GCD file")
parser.add_option("--REALRUN", default=False, action="store_true",
				dest="REALRUN", help="Do real run. (otherwise just test file)")
(options,args) = parser.parse_args()

if not options.identifier_out:
	options.identifier_out = options.identifier

if not options.ram:
    options.ram = '2000'

identifier		= options.identifier
identifier_out	= options.identifier_out
limit 			= options.limit
ram 			= options.ram
gcdfile 		= options.gcdfile
realrun 		= options.REALRUN
datatype        = options.datatype


###### Paths to set manually ######
in_location = '/data/ana/BSM/HNL/MC/'
dag_location_path = '/scratch/lfischer/dagman/'
###### End ######

mother_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
submit_location = os.path.join(mother_path, 'submit')

if not os.path.exists(dag_location_path):
	os.makedirs(dag_location_path)

dag_location = dag_location_path
 
outdir_base = in_location + identifier_out+'/'+'Gen_Weighted'
if not os.path.exists(outdir_base):
	os.makedirs(outdir_base)
	os.chmod(outdir_base, 0o775)


###### Write dag file ######
counter = 0
subf = os.path.join(dag_location, 'dagman-Gen_Weighted_{}_v{}.dag'.format(identifier_out, counter))

while(True):
	if os.path.isfile(subf):
		counter +=1
		subf = subf.split('_v')[0] + '_v{}.dag'.format(counter)
	else:
		break

subtarget = open(subf, 'w')
os.chmod(subf, 0o775)


job_number = 1

infolders = sorted(glob.glob(in_location + identifier+ '/Gen/[0-9]*'))

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
		if limit and file_number >= int(limit.split('-')[0]):
			continue
		elif limit and file_number <= int(limit.split('-')[1]):break  # assuming that job_numbers rise monotonic
		else:
			outfile = os.path.join(this_outdir, (os.path.basename(infile)).replace('Gen','Gen_Weighted'))
			if (not os.path.exists(outfile) or os.path.getsize(outfile) < 4000):
				seed = str(int(infile.split('_')[-1].split('.')[0]))
				subtarget.write('JOB\tjob'+str(job_number)+'\t'+ submit_location + '/submit-Gen_Weighting.condor')
				subtarget.write("\n")
				subtarget.write('VARS\tjob'+str(job_number)+'\tgcdfile="' + str(gcdfile) + '"')
				subtarget.write("\n")
				subtarget.write('VARS\tjob'+str(job_number)+'\tinfile="' + str(infile) + '"')
				subtarget.write("\n")
				subtarget.write('VARS\tjob'+str(job_number)+'\toutfile="' + str(outfile) + '"')
				subtarget.write("\n")
				subtarget.write('VARS\tjob'+str(job_number)+'\tdatatype="' + str(datatype) + '"')
				subtarget.write("\n")
				subtarget.write('VARS\tjob'+str(job_number)+'\tram="' + ram +'"')
				subtarget.write("\n")
				subtarget.write('VARS\tjob'+str(job_number)+'\tidentifier_out="' + identifier_out +'"')
				subtarget.write("\n")	
				subtarget.write('RETRY\tjob'+str(job_number)+'\t2')
				subtarget.write("\n")
				job_number += 1
		if not realrun:break
	if not realrun:break
subtarget.close()
###### End ######


print('We setup: '+str(job_number-1)+' files')
if (job_number - 1) == 0:
    os.remove(subf)
    print('No jobs, removing dag/info file')
else:
    print('Dag file: '+subf)
    print('Submit with: "condor_submit_dag {}"'.format(subf))
    print('Before submission you need to:\nsource /cvmfs/icecube.opensciencegrid.org/users/BeyondStandardModel/bsm-py2-v3.1.1-1/setup.sh\n/cvmfs/icecube.opensciencegrid.org/users/lfischer/oscnext_meta_V01-00-05/build/env-shell.sh')
print('done...')

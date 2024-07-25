#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/icetray-start
#METAPROJECT /cvmfs/icecube.opensciencegrid.org/users/lfischer/oscnext_hnl_meta/build


# Before submission you need to (or before submitting dag which takes the current environment):
# source /cvmfs/icecube.opensciencegrid.org/users/BeyondStandardModel/bsm-py2-v3.1.1-1/setup.sh
# /cvmfs/icecube.opensciencegrid.org/users/lfischer/oscnext_meta_V01-00-05/build/env-shell.sh


#
# Author: Leander Fischer
#

import os, glob
from optparse import OptionParser

usage = "usage: %prog [options] inputfile"
parser = OptionParser(usage)
parser.add_option(
    "-i",
    "--infolder",
    type="string",
    default="/data/ana/BSM/HNL/MC/190607/Gen",
    dest="INFOLDER",
    help="Path to top level Gen File directory.",
)
parser.add_option(
    "--REALRUN",
    default=False,
    action="store_true",
    dest="REALRUN",
    help="Do real run. (otherwise just test file)",
)
(options, args) = parser.parse_args()

in_location = options.INFOLDER
realrun = options.REALRUN

###### Paths to set manually ######
dag_location_path = "/scratch/lfischer/dagman/"
###### End ######

mother_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
submit_location = os.path.join(mother_path, "submit")

if not os.path.exists(dag_location_path):
    os.makedirs(dag_location_path)

dag_location = dag_location_path


###### Write dag file ######
counter = 0
subf = os.path.join(dag_location, "dagman-write_hdf5_Gen_v{}.dag".format(counter))

while True:
    if os.path.isfile(subf):
        counter += 1
        subf = subf.split("_v")[0] + "_v{}.dag".format(counter)
    else:
        break

subtarget = open(subf, "w")
os.chmod(subf, 0o775)

job_number = 1

infolders = sorted(glob.glob(in_location + "/[0-9]*"))

for folder in infolders:
    infiles = sorted(glob.glob(folder + "/*.i3.zst"))
    if len(infiles) == 0:
        infiles = sorted(glob.glob(folder + "/*.i3.bz2"))

    for infile in infiles:
        file_number = int(infile.split("Gen_")[-1].split(".i3")[0])
        subtarget.write(
            "JOB\tjob"
            + str(job_number)
            + "\t"
            + submit_location
            + "/submit-write_hdf5_Gen.condor"
        )
        subtarget.write("\n")
        subtarget.write(
            "VARS\tjob" + str(job_number) + '\tinfile="' + str(infile) + '"'
        )
        subtarget.write("\n")
        subtarget.write("RETRY\tjob" + str(job_number) + "\t2")
        subtarget.write("\n")
        job_number += 1
        if not realrun:
            break
    if not realrun:
        break
subtarget.close()
###### End ######

print("We setup: " + str(job_number - 1) + " files")
if (job_number - 1) == 0:
    os.remove(subf)
    os.remove(phot_data_file)
    print("No jobs, removing dag/info file")
else:
    print("Dag file: " + subf)
    print('Submit with: "condor_submit_dag {}"'.format(subf))
    print(
        "Before submission you need to:\nsource /cvmfs/icecube.opensciencegrid.org/users/BeyondStandardModel/bsm-py2-v3.1.1-1/setup.sh\n/cvmfs/icecube.opensciencegrid.org/users/lfischer/oscnext_meta_V01-00-05/build/env-shell.sh"
    )
print("done...")

#!/bin/bash

source /cvmfs/icecube.opensciencegrid.org/users/BeyondStandardModel/bsm-py2-v3.1.1-1/setup.sh
I3_BUILD=/cvmfs/icecube.opensciencegrid.org/users/lfischer/oscnext_hnl_meta/build/
i3env=${I3_BUILD}/env-shell.sh
SCRIPT=$I3_BUILD/oscNext/resources/scripts/run_oscNext.py

FILEIN=$1
FILEOUT=$2
GCD=$3
HDF5FILE=$4

FILE_L5=${FILEIN}
FILE_L6=${FILEOUT}

$i3env $SCRIPT \
--gcdfile ${GCD} \
--inputfile ${FILE_L5} \
--outputfile ${FILE_L6} \
--datatype LeptonInjector --levels 6 \
--hdf5file ${HDF5FILE} \
# --gridftp

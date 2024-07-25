#!/bin/bash

FILE=$1
OUTDIR=$2
OUTFILE=$3
REALSIM=$4

echo singularity exec -B /cvmfs/icecube.opensciencegrid.org:/cvmfs/icecube.opensciencegrid.org -B /data:/data -B /home/lfischer:/home/lfischer --nv /cvmfs/icecube.opensciencegrid.org/users/lfischer/FLERCNN_evaluate/icetray_stable-tensorflow.sif sh -c '/usr/local/icetray/env-shell.sh python /data/user/lfischer/hnl_analysis/submission_scripts/process/process_flercnn_L9.py -i '${FILE}' -o '${OUTDIR}' --name '${OUTFILE}' --model_dir /cvmfs/icecube.opensciencegrid.org/users/lfischer/FLERCNN_evaluate/ --modelname_list XYZ_ending_FLERCNN --variable_list ending --factor_list 1 --cleaned True --real_sim '${REALSIM}''

singularity exec -B /cvmfs/icecube.opensciencegrid.org:/cvmfs/icecube.opensciencegrid.org -B /data:/data -B /home/lfischer:/home/lfischer --nv /cvmfs/icecube.opensciencegrid.org/users/lfischer/FLERCNN_evaluate/icetray_stable-tensorflow.sif sh -c '/usr/local/icetray/env-shell.sh python /data/user/lfischer/hnl_analysis/submission_scripts/process/process_flercnn_L9.py -i '${FILE}' -o '${OUTDIR}' --name '${OUTFILE}' --model_dir /cvmfs/icecube.opensciencegrid.org/users/lfischer/FLERCNN_evaluate/ --modelname_list XYZ_ending_FLERCNN --variable_list ending --factor_list 1 --cleaned True --real_sim '${REALSIM}''

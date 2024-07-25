#!/bin/bash

FILE=$1
OUTDIR=$2
OUTFILE=$3

echo singularity exec -B /cvmfs/icecube.opensciencegrid.org:/cvmfs/icecube.opensciencegrid.org -B /data:/data -B /home/jmicallef:/home/jmicallef --nv /cvmfs/icecube.opensciencegrid.org/users/jmicallef/FLERCNN_evaluate/icetray_stable-tensorflow.sif sh -c '/usr/local/icetray/env-shell.sh python /cvmfs/icecube.opensciencegrid.org/users/jmicallef/FLERCNN_evaluate/CNN_Test_i3.py -i '${FILE}' -o '${OUTDIR}' --name '${OUTFILE}' --model_dir /cvmfs/icecube.opensciencegrid.org/users/jmicallef/FLERCNN_evaluate/ --modelname_list energy_FLERCNN PID_FLERCNN zenith_FLERCNN Vertex_XYZ_FLERCNN  nDOM  muonV3_FLERCNN --variable_list energy class zenith vertex nDOM muonV3 --factor_list 100 1 1 1 1 1 --cleaned True'

singularity exec -B /cvmfs/icecube.opensciencegrid.org:/cvmfs/icecube.opensciencegrid.org -B /data:/data -B /home/jmicallef:/home/jmicallef --nv /cvmfs/icecube.opensciencegrid.org/users/jmicallef/FLERCNN_evaluate/icetray_stable-tensorflow.sif sh -c '/usr/local/icetray/env-shell.sh python /cvmfs/icecube.opensciencegrid.org/users/jmicallef/FLERCNN_evaluate/CNN_Test_i3.py -i '${FILE}' -o '${OUTDIR}' --name '${OUTFILE}' --model_dir /cvmfs/icecube.opensciencegrid.org/users/jmicallef/FLERCNN_evaluate/ --modelname_list energy_FLERCNN PID_FLERCNN zenith_FLERCNN Vertex_XYZ_FLERCNN nDOM muonV3_FLERCNN --variable_list energy class zenith vertex nDOM muonV3 --factor_list 100 1 1 1 1 1 --cleaned True'

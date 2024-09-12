#!/bin/bash
#SBATCH -J HNL_photon_level_test
#SBATCH -c 1 # Number of cores
#SBATCH -p arguelles_delgado_gpu #gpu_requeue #
#SBATCH --gres=gpu:1 #nvidia_a100_1g.10gb:
#SBATCH --mem 10000
#SBATCH -t 0-01:00:00 # Maximum execution time (D-HH:MM:SS)
#SBATCH -o /n/holylfs05/LABS/arguelles_delgado_lab/Lab/HNL_MC/logging/output/hnl_phot_%A_%a_.out # Standard output
#SBATCH -e /n/holylfs05/LABS/arguelles_delgado_lab/Lab/HNL_MC/logging/error/hnl_phot_%A_%a_.err # Standard error
## SBATCH --array ${FILE_N} # Index of file to produce.

source /n/holylfs05/LABS/arguelles_delgado_lab/Lab/common_software/setup.sh

CMD="source /n/holylfs05/LABS/arguelles_delgado_lab/Everyone/jbook/i3/build/env-shell.sh python /n/holylfs05/LABS/arguelles_delgado_lab/Everyone/jbook/HNL_search_v2/submission_scripts/process/process_Phot.py --runnumber ${SEED} --infile=${INFILE_fullpath} --outfile ${OUTFILE_fullpath} -t"
echo $CMD
$CMD

## Uses default GCD file, DOM efficiency, icemodel, holeice model, crossenergy
## REMEMBER! Any parameters added here need to be passed by the snakefile
## The runnumber in this lingo is the same as the random seed - not the run number as defined by job submission scripts

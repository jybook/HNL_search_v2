#!/bin/bash
#SBATCH -J HNL_det_level_test
#SBATCH -c 1 # Number of cores
#SBATCH -p shared
#SBATCH --mem 10000
#SBATCH -t 0-05:00:00 # Maximum execution time (D-HH:MM:SS)
#SBATCH -o /n/holylfs05/LABS/arguelles_delgado_lab/Lab/HNL_MC/logging/output/hnl_det_%A_%a_.out # Standard output
#SBATCH -e /n/holylfs05/LABS/arguelles_delgado_lab/Lab/HNL_MC/logging/error/hnl_det_%A_%a_.err # Standard error
#SBATCH --array 1-5 # Indeces of files to produce. Start at one so seed is never zero.

#SBATCH --mail-user=$USER@g.harvard.edu  
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


source /n/holylfs05/LABS/arguelles_delgado_lab/Lab/common_software/setup.sh

OUTPUT_PREFIX=/n/holylfs05/LABS/arguelles_delgado_lab/Lab/HNL_MC/${SET_N}/Det
INPUT_PREFIX=/n/holylfs05/LABS/arguelles_delgado_lab/Lab/HNL_MC/${SET_N}/Phot

SEED=${RUN_N}${SLURM_ARRAY_TASK_ID}

if [[ ! -d "${OUTPUT_PREFIX}" ]]; then
    mkdir -p ${OUTPUT_PREFIX}
fi

OUTFILE_fullpath=${OUTPUT_PREFIX}/${RUN_N}_${SLURM_ARRAY_TASK_ID}.i3.zst
INFILE_fullpath=${INPUT_PREFIX}/${RUN_N}_${SLURM_ARRAY_TASK_ID}.i3.zst

CMD="source /n/holylfs05/LABS/arguelles_delgado_lab/Everyone/jbook/i3/build/env-shell.sh python /n/holylfs05/LABS/arguelles_delgado_lab/Everyone/jbook/HNL_search_v2/submission_scripts/process/process_Det.py --runnumber=${SEED} --infile=${INFILE_fullpath} --outfile=${OUTFILE_fullpath}"
echo $CMD
$CMD

## Uses default filenumber, GCD file, DOM efficiency, icemodel, noise model (vuvuzela), holeice model
## The runnumber in this lingo is the same as the random seed - not the run number as defined by job submission scripts

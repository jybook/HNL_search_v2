#!/bin/bash
#SBATCH -J HNL_simulation_test
#SBATCH -c 1 # Number of cores
#SBATCH -p arguelles_delgado_gpu
#SBATCH --gres=gpu:nvidia_a100_1g.10gb:1
#SBATCH --mem 10000
#SBATCH -t 0-00:01:00 # Maximum execution time (D-HH:MM:SS)
#SBATCH -o /n/holylfs05/LABS/arguelles_delgado_lab/Lab/HNL_MC/logging/output/hnl_sim_%A_%a_.out # Standard output
#SBATCH -e /n/holylfs05/LABS/arguelles_delgado_lab/Lab/HNL_MC/logging/error/hnl_sim_%A_%a_.err # Standard error
#SBATCH --array 1-5 # Indeces of files to produce. Start at one so seed is never zero.

echo "doing something"
source /n/holylfs05/LABS/arguelles_delgado_lab/Lab/common_software/setup.sh

echo "things sourced; still running"

OUTPUT_PREFIX=/n/holylfs05/LABS/arguelles_delgado_lab/Lab/HNL_MC/${OUTDIR}/

if [[ ! -d "${OUTPUT_PREFIX}" ]]; then
    mkdir -p ${OUTPUT_PREFIX}
fi

OUTFILE_fullpath=${OUTPUT_PREFIX}/${RUN_N}_${SLURM_ARRAY_TASK_ID}.i3.zst #_${OUTFILE}
SEED=${RUN_N}${SLURM_ARRAY_TASK_ID}

echo "about to set cmd"

CMD="source /n/holylfs05/LABS/arguelles_delgado_lab/Everyone/jbook/i3/build/env-shell.sh python /n/holylfs05/LABS/arguelles_delgado_lab/Everyone/jbook/HNL_search_v2/submission_scripts/process/process_Gen.py --seed ${SEED} --outfile ${OUTFILE_fullpath} --Emin ${EMIN} --Emax ${EMAX} --index ${INDEX} --nEvents ${nEVENTS} -m ${HNL_MASS}"
echo $CMD
$CMD
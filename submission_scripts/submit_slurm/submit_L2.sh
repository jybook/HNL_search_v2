RUNS=1       # Number of runs to produce with these parameters
SET_N=190200  # Random seed and set ID
OUTDIR=$SET_N/L2
for ((RUN_N=1; RUN_N<$RUNS+1; RUN_N++ ))
do
        cmd="sbatch --export=RUN_N=$RUN_N,SET_N=$SET_N,OUTDIR=$OUTDIR /n/holylfs05/LABS/arguelles_delgado_lab/Everyone/jbook/HNL_search_v2/submission_scripts/submit_slurm/submit_L2_jobs.sbatch"
        echo $cmd
        $cmd
done
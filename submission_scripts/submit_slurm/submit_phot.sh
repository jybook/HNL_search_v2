RUNS=1       # Number of runs to produce with these parameters
SET_N=11112  # Random seed and set ID
OUTDIR=$SEED/Phot
for ((RUN_N=1; RUN_N<$RUNS+1; RUN_N++ ))
do
        cmd="sbatch --export=RUN_N=$RUN_N,SET_N=$SET_N,EMIN=$EMIN,EMAX=$EMAX,INDEX=$INDEX,nEVENTS=$nEVENTS,HNL_MASS=$HNL_MASS,OUTDIR=$OUTDIR,OUTFILE=$RUN_N.i3.zst /n/holylfs05/LABS/arguelles_delgado_lab/Everyone/jbook/HNL_search_v2/submission_scripts/submit_slurm/submit_phot_jobs.sbatch"
        echo $cmd
        $cmd
done
RUNS=1       # Number of runs to produce with these parameters
SET_N=190200  # Random seed and set ID
EMIN=2 
EMAX=10000 
INDEX=2
nEVENTS=5000 # Number of events per file
HNL_MASS=0.3
OUTDIR=$SET_N/Gen
for ((RUN_N=1; RUN_N<$RUNS+1; RUN_N++ ))
do
        cmd="sbatch --export=RUN_N=$RUN_N,SET_N=$SET_N,EMIN=$EMIN,EMAX=$EMAX,INDEX=$INDEX,nEVENTS=$nEVENTS,HNL_MASS=$HNL_MASS,OUTDIR=$OUTDIR,OUTFILE=$RUN_N.i3.zst /n/holylfs05/LABS/arguelles_delgado_lab/Everyone/jbook/HNL_search_v2/submission_scripts/submit_slurm/submit_gen_jobs.sbatch"
        echo $cmd
        $cmd
done
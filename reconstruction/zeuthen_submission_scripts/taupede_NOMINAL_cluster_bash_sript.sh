#!/bin/zsh

#(otherwise the default shell would be used)
#$ -S /bin/bash

#(the cpu time for this job, less than 30min to start quickly) use -l h_rt instead (as qsub option)
#####$ -l h_cpu=00:29:59
#####$ -l h_cpu=11:59:59
#####$ -l h_cpu=23:59:59
#####$ -l h_cpu=47:59:59

#(the maximum memory usage of this job)
#$ -l h_rss=4G

#(stderr and stdout are merged together to stdout)
#$ -j y

# execute job from current directory and not relative to your home directory
#$ -cwd

# send email on abort
#$ -m a

#$ -l tmpdir_size=1G

#log path is now set as qsub option as well
#####$ -o /lustre/fs22/group/icecube/lfischer/LOGS/taupede/cluster

# set error log path (just in case, you know)
#$ -e /lustre/fs22/group/icecube/lfischer/LOGS/taupede/cluster/error

echo "STARTING"
echo 'Setting up software ...   '
unset OS_ARCH 
eval `/cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/setup.sh`
echo "LOADED CVMFS SOFTWARE"

echo $(($SGE_TASK_ID-1))
FILE_NR=`expr $(($SGE_TASK_ID-1))`
FILE_NR=`printf "%06d" $FILE_NR`

INPATH=$1
RECOPATH=$2
SIMULATION_SET=$3
FILETYPE=$4

# create out path
OUTPATH=${RECOPATH}'/'${SIMULATION_SET}

# create in-/outfilename
INFILE=${INPATH}'/oscNext_'${FILETYPE}'_level7_v02.00_pass2.'${SIMULATION_SET}'.'${FILE_NR}'.i3.zst'
OUTFILE=${OUTPATH}'/oscNext_'${FILETYPE}'_level7_v02.00_pass2.'${SIMULATION_SET}'.'${FILE_NR}'.i3.zst'

# create reco script path
SCRIPT=${RECOPATH}'/monopod_and_taupede.py'

echo 'input file:'
echo $INFILE

echo 'output file:'
echo $OUTFILE

echo 'the job starts at ' `date`

start=$(date +%s)

I3_ENV=/afs/ifh.de/user/l/lfischer/scratch/ICsoft/meta-projects/oscnext_meta/build/env-shell.sh

echo $I3_ENV python ${SCRIPT} -i ${INFILE} -o ${OUTFILE} -v
$I3_ENV python ${SCRIPT} -i ${INFILE} -o ${OUTFILE} -v

end=$(date +%s)
echo 'The job finishes at' `date`

DIFF=$(( $end - $start ))
echo "It took $DIFF seconds"

exit 0
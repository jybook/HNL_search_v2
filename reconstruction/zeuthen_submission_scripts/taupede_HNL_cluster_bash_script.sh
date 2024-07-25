#!/bin/zsh

#(otherwise the default shell would be used)
#$ -S /bin/bash

#(the cpu time for this job, less than 30min to start quickly)
####$ -l h_cpu=11:59:59
#$ -l h_cpu=00:29:59

#(the maximum memory usage of this job)
#$ -l h_rss=2G

#(stderr and stdout are merged together to stdout)
#$ -j y

# execute job from current directory and not relative to your home directory
#$ -cwd

# send email on abort
#$ -m a

#$ -l tmpdir_size=1G

##### EDIT Log Folder For every run (as with nominal)
#$ -o /lustre/fs22/group/icecube/lfischer/data/HNL/SETNUMBER/RECONAME/LOGS


echo "STARTING"
echo 'Setting up software ...   '
unset OS_ARCH 
eval `/cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/setup.sh`
echo "LOADED CVMFS SOFTWARE"


echo $SGE_TASK_ID
FILE_NR=`expr $SGE_TASK_ID`
FILE_NR=`printf "%05d" $FILE_NR`


# input location
INPATH='/lustre/fs22/group/icecube/lfischer/data/HNL/SterileNeutrino/IC86/HighEnergy/HNL/MC/SETNUMBER/Ares/IC86.AVG/L7/domeff_0.97'


# output location
OUTPATH='/lustre/fs22/group/icecube/lfischer/data/HNL/SETNUMBER/RECO_NAME'


# create in-/outfilename
INFILE=${INPATH}'/L7_00_11_'${FILE_NR}'.i3.zst'
OUTFILE=${OUTPATH}'/files/L7_00_11_'${FILE_NR}'.i3.zst'


echo 'input file:'
echo $INFILE

echo 'output file:'
echo $OUTFILE

echo 'the job starts at ' `date`

start=$(date +%s)

I3_ENV=/afs/ifh.de/user/l/lfischer/scratch/ICsoft/meta-projects/oscnext_meta/build/env-shell.sh

SCRIPT=${OUTPATH}'/taupede.py'

echo $I3_ENV python ${SCRIPT} -i ${INFILE} -o ${OUTFILE} -v
$I3_ENV python ${SCRIPT} -i ${INFILE} -o ${OUTFILE} -v

end=$(date +%s)
echo 'The job finishes at' `date`

DIFF=$(( $end - $start ))
echo "It took $DIFF seconds"

exit 0
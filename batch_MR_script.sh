#!/bin/bash
# Tell SGE that we are using the bash shell
#$ -S /bin/bash

# Example file to create a task-farm of identical jobs on Apollo

# Do not lmit the stacksize (the maximum memory a job can use)
ulimit -s unlimited
# Do not limit the number of open files a job can have
#ulimit -n unlimited
# Run the job from the following directory
#$ -cwd
# Created files will have the fs-virgo group
# This feature seems to be disabled on Apollo, so this does not work
newgrp fs-virgo
# Created files will have rw permission for the group and r for the world
umask 002
module load mps/software/
# Set pathnames below relative to the current working directory
#$ -cwd
# Say which queue you want to submit to
#$ -q mps.q@@compute_amd_c6145_mps
## #$ -q eng-inf_himem.q

#$ -l h=!(node108|node106|node105)
# Define a task farm of jobs
#$ -t 1-512
# Limit to 50 concurrent jobs
#$ -tc 50
# Join standard error to standard out
#$ -j y
# Give the job a name
#$ -N _Dust_MR_all
# Name and location of the output file
# SGE will only substitute certain variables here
#$ -o ./logs/$JOB_NAME_$TASK_ID.log

# The parentheses here allow one to do algebra with shell variables
i=$(($SGE_TASK_ID -1))
echo Running on file $i
ff=$(($i))
lf=$(($i))

# Create personalised input parameter files
echo FirstFile $ff | cat >  input/input_batch/input.MR_$i
echo LastFile $lf | cat >>  input/input_batch/input.MR_$i
echo MaxMemSize 8000 >> input/input_batch/input.MR_$i
cat input/input_MR_W1_PLANCK_apollo.par >> input/input_batch/input.MR_$i

# Run jobs
./L-Galaxies input/input_batch/input.MR_$i

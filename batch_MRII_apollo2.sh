#!/bin/bash
# Tell SGE that we are using the bash shell
#$ -S /bin/bash

# Example file to create a task-farm of identical jobs on Apollo

# Do not lmit the stacksize (the maximum memory a job can use)
ulimit -s unlimited
# Do not limit the number of open files a job can have
#ulimit -n unlimited
# Run the job from the following directory
cd /home/sc558/lgal_clay17/
# Created files will have the fs-virgo group
# This feature seems to be disabled on Apollo, so this does not work
# newgrp fs-virgo
# Created files will have rw permission for the group and r for the world
umask 002

# Set pathnames below relative to the current working directory
#$ -cwd
# Say which queue you want to submit to
#$ -q mps.q@@compute_amd_c6145_mps
# Define a task farm of jobs
#$ -t 3-512
## #$ -t 1

# Limit to 50 concurrent jobs
#$ -tc 50
# Join standard error to standard out
#$ -j y
# Give the job a name
#$ -N Clay17_Dust_MRII
# Name and location of the output file
# SGE will only substitute certain variables here
#$ -o ./logs/$JOB_NAME_$TASK_ID.log

module add sge
module add gcc/4.8.1
#module add openmpi/gcc/64/1.7.3
module add intel-mpi/64/4.1.1/036
module add gsl/gcc/1.15


# The parentheses here allow one to do algebra with shell variables
i=$(($SGE_TASK_ID - 1))
echo Running on file $i
ff=$(($i))
lf=$(($i))

# Create personalised input parameter files
echo FirstFile $ff | cat >  input/input_batch/input.MRII_$i
echo LastFile $lf | cat >>  input/input_batch/input.MRII_$i
echo MaxMemSize 10000 >> input/input_batch/input.MRII_$i
cat input/input_MRII_W1_PLANCK_apollo.par >> input/input_batch/input.MRII_$i

# Run jobs
./L-Galaxies input/input_batch/input.MRII_$i

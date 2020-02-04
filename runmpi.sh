#!/bin/bash

#

#PBS -V

#PBS -N myjob

#PBS -q batch

#PBS -t 1-3

#PBS -m ae

#PBS -M “your email”

 

### Specify the number of cpus for your job.  This example will allocate 8 cores

### using 4 nodes with 8 processes per node.

#PBS -l nodes=4:ppn=8

 

# pmem is amount of memory per processor

#PBS -l pmem=2gb

 

#PBS -l walltime=1:00:00

 

### Switch to the working directory; by default TORQUE launches processes

### from your home directory.
$PBS_O_WORKDIR=/home/yourhome/mpi_files

echo Working directory is $PBS_O_WORKDIR

cd $PBS_O_WORKDIR

 

# The number of processors allocated to this run, NPROCS is used by mpirun.
NPROCS=`wc -l < $PBS_NODEFILE`

 

# The number of nodes used in this run.

NNODES=`uniq $PBS_NODEFILE | wc -l`

 

### Display the job context

echo Running on host `hostname`

echo Time is `date`

echo Directory is `pwd`

echo The following processors are allocated to this job:

echo `cat $PBS_NODEFILE`

 

module load openmpi/version

RANDOM=$(date +%N) 

MPIRUN=`which mpirun`

counter=1
 while [ $counter -le 3 ]
 	do
    PartOne=$(( ${PBS_ARRAYID}*3 ))
    PartTwo=$(( PartOne-3 ))
    MaesJob=$(( PartTwo+$counter ))

	${MPIRUN} -v -machinefile $PBS_NODEFILE -np $NPROCS   ABCSMC Seed=$RANDOM Particles=1 Final=0.1 Job=$MaesJob Alpha=0.8 Nparameters=7 Time=${PBS_ARRAYID}
	((counter++))
	done
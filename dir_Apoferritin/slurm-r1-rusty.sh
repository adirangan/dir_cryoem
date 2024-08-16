#!/bin/bash -l

# Standard output and error:
# Initial working directory:
#SBATCH -D ./
# Job Name:
#SBATCH -J r1
# Queue (Partition):
#SBATCH --constraint=rome
#SBATCH --partition=ccm
#
# Request 1 node(s)
#SBATCH --nodes=1
# Set the number of tasks per node (=MPI ranks)
#SBATCH --ntasks-per-node=16
# Set the number of threads per rank (=OpenMP threads)
#SBATCH --cpus-per-task=8

# Wall clock limit:
#SBATCH --time=48:00:00

module load openmpi/4.0.7

#export KMP_AFFINITY=compact,granularity=core,1
export OMP_STACKSIZE=120M
export OMP_NUM_THREADS=8
export GPU=0
export BIOEM_DEBUG_OUTPUT=0
export BIOEM_ALGO=1

listmod="list_mod"
listpart="list_part"
paramfile="Param_BioEM_Apofer-R1"

mkdir outputs 

for mod in `cat ${listmod}` 
do

mpirun --map-by socket:pe=$OMP_NUM_THREADS  -np 16 ~/BenchmarkBioEM/BioEM_Mod_Gera/bioEM --Particlesfile $listpart --ReadMRC  --ReadMultipleMRC --Modelfile models/$mod --ReadPDB --OutputFile Output_algo1_$mod --Inputfile ${paramfile}  --ReadOrientation Quat_36864 

mv Output_algo1_$mod outputs

done


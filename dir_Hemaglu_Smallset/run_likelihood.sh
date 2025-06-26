#!/bin/sh
#SBATCH -p gpu
#SBATCH -C h100
#SBATCH -n 1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --ntasks-per-node=1
#SBATCH --mem=64G
#SBATCH --time=10:00:00
#SBATCH --job-name=test
#SBATCH --output=slurm_out/%x-%j.out
#SBATCH --error=slurm_err/%x-%j.err
#SBATCH --array=0-1

echo "Starting job $SLURM_JOB_ID at `date` on `hostname`"
source ./load_env_h100.sh
tagtemplate=$1

python3 likelihood_1mod_images.py $SLURM_ARRAY_TASK_ID $tagtemplate

echo "Job finished with exit code $? at: `date`"

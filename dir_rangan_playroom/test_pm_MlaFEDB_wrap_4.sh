#!/bin/bash
module load matlab; cd /mnt/home/rangan/dir_cryoem/dir_rangan_playroom ; matlab -batch "test_pm_MlaFEDB_4;";

## run with: ;
## module load slurm disBatch; chmod +x *.sh; sbatch -p ccm -n 20 -c 2 -t 11:59 disBatch test_pm_MlaFEDB_task_4.txt

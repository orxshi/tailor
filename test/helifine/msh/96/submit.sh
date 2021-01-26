#!/bin/bash
#SBATCH -p single
#SBATCH -J helimesh
#SBATCH -t 01:00:00
#SBATCH --ntasks=1
#SBATCH --array=0-5

bash mesh${SLURM_ARRAY_TASK_ID}.sh 96

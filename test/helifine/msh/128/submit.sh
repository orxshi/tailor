#!/bin/bash
#SBATCH -p short
#SBATCH -J helimesh
#SBATCH -t 04:00:00
#SBATCH --ntasks=8
#SBATCH --array=0-5

bash mesh${SLURM_ARRAY_TASK_ID}.sh 128

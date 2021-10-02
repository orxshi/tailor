#!/bin/bash
#SBATCH -p single
#SBATCH -t 00:10:00
#SBATCH --ntasks=1

bash mesh.sh 32

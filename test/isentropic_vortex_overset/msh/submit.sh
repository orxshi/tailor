#!/bin/bash
#SBATCH -p short
#SBATCH -J helimesh
#SBATCH -t 00:10:00
#SBATCH --ntasks=4

bash mesh.sh 1

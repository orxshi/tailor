#!/bin/bash
#SBATCH -p short
#SBATCH -J helimesh
#SBATCH -t 04:00:00
#SBATCH --ntasks=8

bash mesh.sh 1

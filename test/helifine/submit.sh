#!/bin/bash
#SBATCH -p single
#SBATCH -J helifine
#SBATCH --ntasks=1
#SBATCH --hint=nomultithread
#SBATCH -t 00-01:00:00
#SBATCH --output=slurm-%j.log
#SBATCH --error=slurm-%j.err
mpirun --tag-output --report-bindings /usr/bin/time -f '%e %S %U %P %M' -o "timing.dat" --append ./out

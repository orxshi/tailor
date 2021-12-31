#!/bin/bash
#SBATCH -p mid2
#SBATCH --ntasks=64
#SBATCH -t 01-00:00:00
#SBATCH --hint=nomultithread
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.log
mpirun --tag-output --report-bindings /usr/bin/time -f '%e %S %U %P %M' -o "timing.dat" --append ./out

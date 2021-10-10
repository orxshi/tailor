#!/bin/bash
#SBATCH -p short
##SBATCH -x barbun[2-20]
##SBATCH -x levrek[129-136,138-164]
##SBATCH -x sardalya[2-64,97-150]
#SBATCH -J helifine
#SBATCH --ntasks=64
##SBATCH --nodes=3-4
#SBATCH --hint=nomultithread
#SBATCH -t 00-04:00:00
##SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
##SBATCH --signal=TERM@300
mpirun --tag-output --report-bindings /usr/bin/time -f '%e %S %U %P %M' -o "timing.dat" --append ./out
#valgrind --undef-value-errors=no --leak-check=full \
    #--show-leak-kinds=all \
    #--verbose \
    #--log-file=valgrind-out.txt --suppressions=/truba/sw/centos7.3/lib/openmpi/4.0.1-gcc-7.0.1/share/openmpi/openmpi-valgrind.supp

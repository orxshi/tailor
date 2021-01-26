#!/bin/bash
#SBATCH -p short
##SBATCH -x barbun[2-100]
##SBATCH -x barbun4,barbun5,barbun10,barbun12,barbun13,barbun15,barbun16,barbun17,barbun24,barbun25,barbun34,barbun35,barbun36,barbun37,barbun38,barbun41,barbun53,barbun92
#SBATCH -J helifine
#SBATCH --ntasks=128
#SBATCH --hint=nomultithread
#SBATCH -t 00-01:00:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
##SBATCH --signal=TERM@300
mpirun --tag-output --report-bindings /usr/bin/time -f '%e %S %U %P %M' -o "timing.dat" --append ./out
#valgrind --undef-value-errors=no --leak-check=full \
    #--show-leak-kinds=all \
    #--verbose \
    #--log-file=valgrind-out.txt --suppressions=/truba/sw/centos7.3/lib/openmpi/4.0.1-gcc-7.0.1/share/openmpi/openmpi-valgrind.supp

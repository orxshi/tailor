#!/bin/bash
#SBATCH -p mid1
##SBATCH -x barbun[2-80]
##SBATCH -x barbun2,barbun3,barbun4,barbun5,barbun6,barbun9,barbun10,barbun12,barbun13,barbun15,barbun16,barbun17,barbun19,barbun20,barbun21,barbun24,barbun25,barbun28,barbun31,barbun33,barbun34,barbun35,barbun36,barbun37,barbun38,barbun39,barbun41,barbun42,barbun43,barbun48,barbun50,barbun51,barbun52,barbun53,barbun54,barbun58,barbun62,barbun63,barbun64,barbun71,barbun92
#SBATCH -J helifine
#SBATCH --ntasks=128
#SBATCH --hint=nomultithread
#SBATCH -t 04-00:00:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
##SBATCH --signal=TERM@300
mpirun --tag-output --report-bindings /usr/bin/time -f '%e %S %U %P %M' -o "timing.dat" --append ./out
#valgrind --undef-value-errors=no --leak-check=full \
    #--show-leak-kinds=all \
    #--verbose \
    #--log-file=valgrind-out.txt --suppressions=/truba/sw/centos7.3/lib/openmpi/4.0.1-gcc-7.0.1/share/openmpi/openmpi-valgrind.supp

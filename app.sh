###========================================
#!/bin/bash
#BSUB -n 2 
#BSUB -o lsf.out
#BSUB -e lsf.err
#BSUB -q "standard"
#BSUB -J mpi
#BSUB -x
#---------------------------------------------------------------------
#BSUB -R gpu
#BSUB -R "span[ptile=1]"
. /etc/profile.d/modules.sh
module load cuda/7.0.28
module load openmpi/1.10.2

mpirun -np 2 ./mpi_implementation output > output.txt
###end of script

#PBS -S /bin/bash
#PBS -q newest
#PBS -o column.out
#PBS -e column.err
#PBS -N column
#PBS -l nodes=1:ppn=2
#PBS -l walltime=00:00:10

cd $PBS_O_WORKDIR
mpirun -np 4 column.o

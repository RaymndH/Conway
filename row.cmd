#PBS -S /bin/bash
#PBS -q newest
#PBS -o row.out
#PBS -e row.err
#PBS -N row
#PBS -l nodes=1:ppn=2
#PBS -l walltime=00:00:10

cd $PBS_O_WORKDIR
mpirun -np 2 row.o

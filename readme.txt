To run column decomposition:

$ mpif90 -o column.o column.f90 
$ qsub column.cmd
$ vim column.out

To run row decomposition(not working):
$ mpif90 -o row.o row.f90
$ qsub row.cmd
$ vim row.out

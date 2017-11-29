# MPP_DataConversion for [COMPSCI590] Parallel Programming @ Duke University

* Parallel matrix conversion and redistribution

* Authors: Ryan Yoon, Victor Yu


Compile:

   In .src/Makefile, set properly the Fortran compiler,
   compilation flags, BLACS library, OpenMP settings.
   Then "make".

Test:

1) Matrix format conversion

   ./src/test_conv.x arg1 arg2 arg3

   arg1: number of rows in test matrix
   arg2: number of columns in test matrix
   arg3: sparsity of test matrix

2) Matrix redistribution

   mpirun -n np ./src/test_redist.x arg1 arg2

   mpirun: name of MPI executable
   np    : number of MPI tasks
   arg1  : number of rows (= columns) in test matrix
   arg2  : sparsity of test matrix

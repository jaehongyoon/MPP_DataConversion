program test_redist

   use mat_redist
   use mat_util

   implicit none

   include 'mpif.h'

   real(kind=8), dimension(:,:), allocatable :: a,a2

   real(kind=8) :: threshold
   real(kind=8) :: sparsity
   real(kind=8) :: t0, t1
   integer :: n, nnz, g_nnz
   integer :: i_row, i_col
   character(len=128) :: arg1, arg2

   integer :: n_proc, nprow, npcol, myid, myprow, mypcol
   integer :: mpi_comm_global, mpierr
   integer :: blk
   integer :: BLACS_CTXT
   integer :: local_row, local_col
   integer, external :: numroc

   integer :: status

   ! Initialize MPI
   call MPI_Init(mpierr)
   mpi_comm_global = MPI_COMM_WORLD
   call MPI_Comm_size(mpi_comm_global,n_proc,mpierr)
   call MPI_Comm_rank(mpi_comm_global,myid,mpierr)

   if(n_proc < 2) then
      write(*,"('Test with at least 2 MPI tasks!')")
      call MPI_Abort
   endif

   ! Read input
   if(COMMAND_ARGUMENT_COUNT() == 2) then
      call GET_COMMAND_ARGUMENT(1,arg1)
      call GET_COMMAND_ARGUMENT(2,arg2)
      read(arg1,*) n
      read(arg2,*) sparsity
   else
      if(myid == 0) then
         write(*,"('Wrong number of arguments!')")
         write(*,"('Expected:')")
         write(*,"(' 1) Number of rows')")
         write(*,"(' 2) Desired sparsity')")
      endif
      call MPI_Abort
   endif

   threshold = 1.0d-15

   call init_timer()

   if(myid == 0) then
      write(*,"('================================================')")
      write(*,"('== RUNNING MATRIX REDISTRIBUTION TEST PROGRAM ==')")
      write(*,"('================================================')")
      write(*,"('Generating dense random matrix..')")
      write(*,"('|  Number of rows      : ',I16)") n
      write(*,"('|  Number of columns   : ',I16)") n
      write(*,"('|  Zero threshold      : ',E16.1)") threshold
      write(*,"('|  Expected sparsity   : ',F16.3)") sparsity
   endif

   ! Set up square-like processor grid
   do npcol = nint(sqrt(real(n_proc))),2,-1
      if(mod(n_proc,npcol) == 0) exit
   enddo
   nprow = n_proc/npcol

   ! Set block size
   blk = 16

   ! Set up BLACS
   BLACS_CTXT = mpi_comm_global

   call BLACS_Gridinit(BLACS_CTXT,'r',nprow,npcol)
   call BLACS_Gridinfo(BLACS_CTXT,nprow,npcol,myprow,mypcol)

   local_row = numroc(n,blk,myprow,0,nprow)
   local_col = numroc(n,blk,mypcol,0,npcol)

   allocate(a(local_row,local_col))
   allocate(a2(local_row,local_col))

   ! Generate random matrix and number of non-zeros
   call get_wall_time(t0)
   call rand_mat(sparsity,local_row,local_col,a)
   call get_wall_time(t1)
   call get_nnz(threshold,local_row,local_col,a,nnz)

   ! Save a copy of a for correctness check
   a2 = a

   ! Compute sparsity
   sparsity = 1.d0-dble(nnz)/(n*n)

   if(myid == 0) then
      write(*,"('|  Number of non-zeros : ',I16)") nnz
      write(*,"('|  Actual sparsity     : ',F16.3)") sparsity
      write(*,"('|  Time                : ',F15.2,'s')") t1-t0
      write(*,"('================================================')")
   endif

   ! 2D block-cyclic DNS to 1D block CSC
   if(myid == 0) then
      write(*,"('Conversion: 2D block-cyclic ==> 1D block')")
   endif
   call get_wall_time(t0)
   call redist_2dbc_to_1db(threshold,n,local_row,local_col,a,&
                           blk,mpi_comm_global,BLACS_CTXT)
   call get_wall_time(t1)
   if(myid == 0) then
      write(*,"('|  Time : ',F8.3,'s')") t1-t0
   endif

   ! 1D block CSC to 2D block-cyclic DNS
   if(myid == 0) then
      write(*,"('Conversion: 1D block ==> 2D block-cyclic')")
   endif
   call get_wall_time(t0)
   call redist_1db_to_2dbc(n,local_row,local_col,a2,&
                           blk,mpi_comm_global,BLACS_CTXT)
   call get_wall_time(t1)
   if(myid == 0) then
      write(*,"('|  Time : ',F8.3,'s')") t1-t0
      write(*,"('================================================')")
   endif

   a = a-a2

   call get_nnz(threshold,local_row,local_col,a,nnz)

   call MPI_Reduce(nnz,g_nnz,1,mpi_integer,mpi_sum,0,mpi_comm_global,mpierr)

   if(myid == 0) then
      if(g_nnz > 0) then
         write(*,"('|  Failed!!')")
      else
         write(*,"('|  Passed!!')")
      endif
      write(*,"('================================================')")
   endif

   deallocate(a)
   deallocate(a2)

   call MPI_Finalize(mpierr)

end program

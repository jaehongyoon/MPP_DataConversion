program test_conv

   use mat_conv
   use mat_util

   implicit none

   real(kind=8), dimension(:,:), allocatable :: a,a2
   real(kind=8), dimension(:), allocatable :: val
   integer, dimension(:), allocatable :: row_idx, col_idx, row_ptr, col_ptr

   real(kind=8) :: threshold
   real(kind=8) :: sparsity
   real(kind=8) :: t0, t1
   integer :: m, n, nnz
   integer :: i_row, i_col
   character(len=128) :: arg1, arg2, arg3

   integer :: status

   ! Read input
   if(COMMAND_ARGUMENT_COUNT() == 3) then
      call GET_COMMAND_ARGUMENT(1,arg1)
      call GET_COMMAND_ARGUMENT(2,arg2)
      call GET_COMMAND_ARGUMENT(3,arg3)
      read(arg1,*) m
      read(arg2,*) n
      read(arg3,*) sparsity
   else
      write(*,"('Wrong number of arguments!')")
      write(*,"('Expected:')")
      write(*,"(' 1) Number of rows')")
      write(*,"(' 2) Number of columns')")
      write(*,"(' 3) Desired sparsity')")
      stop
   endif

   threshold = 1.0d-15

   call init_timer()

   write(*,"('============================================')")
   write(*,"('== RUNNING MATRIX CONVERSION TEST PROGRAM ==')")
   write(*,"('============================================')")
   write(*,"('Generating dense random matrix..')")
   write(*,"('|  Number of rows      : ',I16)") m
   write(*,"('|  Number of columns   : ',I16)") n
   write(*,"('|  Zero threshold      : ',E16.1)") threshold
   write(*,"('|  Expected sparsity   : ',F16.3)") sparsity

   allocate(a(m,n))
   allocate(a2(m,n))

   ! Generate random matrix and number of non-zeros
   call get_wall_time(t0)
   call rand_mat(sparsity,m,n,a)
   call get_wall_time(t1)
   call get_nnz(threshold,m,n,a,nnz)

   ! Save a copy of a for correctness check
   a2 = a

   ! Compute sparsity
   sparsity = 1.d0-dble(nnz)/(m*n)

   write(*,"('|  Number of non-zeros : ',I16)") nnz
   write(*,"('|  Actual sparsity     : ',F16.3)") sparsity
   write(*,"('|  Time                : ',F15.2,'s')") t1-t0
   write(*,"('============================================')")

   allocate(val(nnz))
   allocate(row_idx(nnz))
   allocate(col_idx(nnz))
   allocate(row_ptr(m+1))
   allocate(col_ptr(n+1))

   ! DNS to CCS
   write(*,"('Conversion: DNS ==> CCS')")
   call get_wall_time(t0)
   call dns_to_ccs(threshold,m,n,nnz,a,val,row_idx,col_ptr)
   call get_wall_time(t1)
   write(*,"('|  Time : ',F8.3,'s')") t1-t0

   ! CCS to DNS
   write(*,"('Conversion: CCS ==> DNS')")
   call get_wall_time(t0)
   call ccs_to_dns(m,n,nnz,val,row_idx,col_ptr,a2)
   call get_wall_time(t1)
   write(*,"('|  Time : ',F8.3,'s')") t1-t0
   write(*,"('============================================')")

   a = a-a2

   call get_nnz(threshold,m,n,a,nnz)

   if(nnz > 0) then
      write(*,"('|  Failed!!')")
   else
      write(*,"('|  Passed!!')")
   endif
   write(*,"('============================================')")

   ! DNS to CRS
   write(*,"('Conversion: DNS ==> CRS')")
   call get_wall_time(t0)
   call dns_to_crs(threshold,m,n,nnz,a,val,col_idx,row_ptr)
   call get_wall_time(t1)
   write(*,"('|  Time : ',F8.3,'s')") t1-t0

   ! CRS to DNS
   write(*,"('Conversion: CRS ==> DNS')")
   call get_wall_time(t0)
   call crs_to_dns(m,n,nnz,val,col_idx,row_ptr,a2)
   call get_wall_time(t1)
   write(*,"('|  Time : ',F8.3,'s')") t1-t0
   write(*,"('============================================')")

   a = a-a2

   call get_nnz(threshold,m,n,a,nnz)

   if(nnz > 0) then
      write(*,"('|  Failed!!')")
   else
      write(*,"('|  Passed!!')")
   endif
   write(*,"('============================================')")

   ! DNS to COO
   write(*,"('Conversion: DNS ==> COO')")
   call get_wall_time(t0)
   call dns_to_coo(threshold,m,n,nnz,a,val,row_idx,col_idx)
   call get_wall_time(t1)
   write(*,"('|  Time : ',F8.3,'s')") t1-t0

   ! COO to DNS
   write(*,"('Conversion: COO ==> DNS')")
   call get_wall_time(t0)
   call crs_to_dns(m,n,nnz,val,row_idx,col_idx,a2)
   call get_wall_time(t1)
   write(*,"('|  Time : ',F8.3,'s')") t1-t0
   write(*,"('============================================')")

   a = a-a2

   call get_nnz(threshold,m,n,a,nnz)

   if(nnz > 0) then
      write(*,"('|  Failed!!')")
   else
      write(*,"('|  Passed!!')")
   endif
   write(*,"('============================================')")

   deallocate(a)
   deallocate(a2)
   deallocate(val)
   deallocate(row_idx)
   deallocate(col_idx)
   deallocate(row_ptr)
   deallocate(col_ptr)

   call exit

end program

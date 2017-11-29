module mat_util

   implicit none

   private

   public :: rand_mat
   public :: get_nnz
   public :: get_global_idx
   public :: init_timer
   public :: get_wall_time

   integer :: clock_rate
   integer :: clock_max

contains

subroutine rand_mat(sparsity,n_rows,n_cols,matrix)

   implicit none

   real(kind=8), intent(in)  :: sparsity
   integer,      intent(in)  :: n_rows
   integer,      intent(in)  :: n_cols
   real(kind=8), intent(out) :: matrix(n_rows,n_cols) !< Dense matrix

   integer :: i_col
   integer :: i_row

   ! Gaussian random number generation
   call random_number(matrix)

   ! Make it sparse
   !$omp do
   do i_col = 1,n_cols
      do i_row = 1,n_rows
         if(matrix(i_row,i_col) < sparsity) then
            matrix(i_row,i_col) = 0.d0
         endif
      enddo
   enddo
   !$omp end do

end subroutine

! This routine computes the local number of non_zero elements.
subroutine get_nnz(zero_def,n_rows,n_cols,matrix,nnz)

   implicit none

   real(kind=8), intent(in)  :: zero_def
   integer,      intent(in)  :: n_rows
   integer,      intent(in)  :: n_cols
   real(kind=8), intent(in)  :: matrix(n_rows,n_cols) !< Dense matrix
   integer,      intent(out) :: nnz !< Number of non-zeros

   integer :: i_row
   integer :: i_col

   nnz = 0

   !$omp do reduction( + : nnz )
   do i_col = 1,n_cols
      do i_row = 1,n_rows
         if(abs(matrix(i_row,i_col)) > zero_def) then
            nnz = nnz+1
         endif
      enddo
   enddo
   !$omp end do

end subroutine

! This routine computes the global index from local index
subroutine get_global_idx(global_idx,local_idx,block_size,pid,np)

   implicit none

   integer, intent(out) :: global_idx !< Global index
   integer, intent(in)  :: local_idx  !< Local index
   integer, intent(in)  :: block_size
   integer, intent(in)  :: pid
   integer, intent(in)  :: np

   integer :: block !< Local block
   integer :: idx   !< Local index in block

   block = (local_idx-1)/block_size
   idx = local_idx-block*block_size

   global_idx = pid*block_size+block*block_size*np+idx

end subroutine

! This routine initializes system clock.
subroutine init_timer()

   implicit none

   integer :: initial_time

   call system_clock(initial_time,clock_rate,clock_max)

end subroutine

! This routine gets the wall time.
subroutine get_wall_time(wtime)

   implicit none
   real*8, intent(out) :: wtime
 
   integer :: tics

   call system_clock(tics)

   wtime = 1.d0*tics/clock_rate

end subroutine

end module

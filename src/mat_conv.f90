module mat_conv

   implicit none

   private

   ! fast: sparsity pattern is known
   public :: dns_to_ccs
   public :: dns_to_ccs_fast
   public :: dns_to_crs
   public :: dns_to_crs_fast
   public :: dns_to_coo
   public :: dns_to_coo_fast
   public :: ccs_to_dns
   public :: ccs_to_coo
   public :: ccs_to_crs
   public :: crs_to_dns
   public :: crs_to_ccs
   public :: crs_to_coo
   public :: coo_to_dns
   public :: coo_to_ccs
   public :: coo_to_crs

contains

subroutine dns_to_ccs(zero_def,n_rows,n_cols,nnz,matrix,val,row_ind,col_ptr)

   implicit none

   real(kind=8), intent(in)  :: zero_def
   integer,      intent(in)  :: n_rows
   integer,      intent(in)  :: n_cols
   integer,      intent(in)  :: nnz
   real(kind=8), intent(in)  :: matrix(n_rows,n_cols) !< Dense matrix
   real(kind=8), intent(out) :: val(nnz) !< Non-zero values
   integer,      intent(out) :: row_ind(nnz) !< Row index
   integer,      intent(out) :: col_ptr(n_cols+1) !< Column pointer

   integer :: i_row
   integer :: i_col
   integer :: i_val
   logical :: new_col

   i_val = 0
   col_ptr = 0

   do i_col = 1,n_cols
      new_col = .true.
      do i_row = 1,n_rows
         if(abs(matrix(i_row,i_col)) > zero_def) then
            i_val = i_val+1
            if(new_col) then
               col_ptr(i_col) = i_val
               new_col = .false.
            endif
            val(i_val) = matrix(i_row,i_col)
            row_ind(i_val) = i_row
         endif
      enddo
   enddo

   col_ptr(n_cols+1) = i_val+1

end subroutine

subroutine dns_to_crs(zero_def,n_rows,n_cols,nnz,matrix,val,col_ind,row_ptr)

   implicit none

   real(kind=8), intent(in)  :: zero_def
   integer,      intent(in)  :: n_rows
   integer,      intent(in)  :: n_cols
   integer,      intent(in)  :: nnz
   real(kind=8), intent(in)  :: matrix(n_rows,n_cols) !< Dense matrix
   real(kind=8), intent(out) :: val(nnz) !< Non-zero values
   integer,      intent(out) :: col_ind(nnz) !< Column index
   integer,      intent(out) :: row_ptr(n_rows+1) !< Row pointer

   integer :: i_row
   integer :: i_col
   integer :: i_val
   logical :: new_row

   i_val = 0
   row_ptr = 0

   do i_row = 1,n_rows
      new_row = .true.
      do i_col = 1,n_cols
         if(abs(matrix(i_row,i_col)) > zero_def) then
            i_val = i_val+1
            if(new_row) then
               row_ptr(i_row) = i_val
               new_row = .false.
            endif
            val(i_val) = matrix(i_row,i_col)
            col_ind(i_val) = i_col
         endif
      enddo
   enddo

   row_ptr(n_rows+1) = i_val+1

end subroutine

subroutine dns_to_coo(zero_def,n_rows,n_cols,nnz,matrix,val,row_ind,col_ind)

   implicit none

   real(kind=8), intent(in)  :: zero_def
   integer,      intent(in)  :: n_rows
   integer,      intent(in)  :: n_cols
   integer,      intent(in)  :: nnz
   real(kind=8), intent(in)  :: matrix(n_rows,n_cols) !< Dense matrix
   real(kind=8), intent(out) :: val(nnz) !< Non-zero values
   integer,      intent(out) :: row_ind(nnz) !< Row index
   integer,      intent(out) :: col_ind(nnz) !< Column index

   integer :: i_row
   integer :: i_col
   integer :: i_val

   i_val = 0

   do i_col = 1,n_cols
      do i_row = 1,n_rows
         if(abs(matrix(i_row,i_col)) > zero_def) then
            i_val = i_val+1
            val(i_val) = matrix(i_row,i_col)
            row_ind(i_val) = i_row
            col_ind(i_val) = i_col
         endif
      enddo
   enddo

end subroutine

subroutine dns_to_ccs_fast(n_rows,n_cols,nnz,matrix,row_ind,col_ptr,val)

   implicit none

   integer,      intent(in)  :: n_rows
   integer,      intent(in)  :: n_cols
   integer,      intent(in)  :: nnz
   real(kind=8), intent(in)  :: matrix(n_rows,n_cols) !< Dense matrix
   integer,      intent(in)  :: row_ind(nnz) !< Row index
   integer,      intent(in)  :: col_ptr(n_cols+1) !< Column pointer
   real(kind=8), intent(out) :: val(nnz) !< Non-zero values

   integer :: i_row
   integer :: i_col
   integer :: i_val

   i_col = 0

   do i_val = 1,nnz
      if(i_val == col_ptr(i_col+1) .and. i_col /= n_cols) then
         i_col = i_col+1
      endif
      i_row = row_ind(i_val)
      val(i_val) = matrix(i_row,i_col)
   enddo

end subroutine

subroutine dns_to_crs_fast(n_rows,n_cols,nnz,matrix,col_ind,row_ptr,val)

   implicit none

   integer,      intent(in)  :: n_rows
   integer,      intent(in)  :: n_cols
   integer,      intent(in)  :: nnz
   real(kind=8), intent(in)  :: matrix(n_rows,n_cols) !< Dense matrix
   integer,      intent(in)  :: col_ind(nnz) !< Column index
   integer,      intent(in)  :: row_ptr(n_rows+1) !< Row pointer
   real(kind=8), intent(out) :: val(nnz) !< Non-zero values

   integer :: i_row
   integer :: i_col
   integer :: i_val

   i_row = 0

   do i_val = 1,nnz
      if(i_val == row_ptr(i_row+1) .and. i_row /= n_rows) then
         i_row = i_row+1
      endif
      i_col = col_ind(i_val)
      val(i_val) = matrix(i_row,i_col)
   enddo

end subroutine

subroutine dns_to_coo_fast(n_rows,n_cols,nnz,matrix,row_ind,col_ind,val)

   implicit none

   integer,      intent(in)  :: n_rows
   integer,      intent(in)  :: n_cols
   integer,      intent(in)  :: nnz
   real(kind=8), intent(in)  :: matrix(n_rows,n_cols) !< Dense matrix
   integer,      intent(in)  :: row_ind(nnz) !< Row index
   integer,      intent(in)  :: col_ind(nnz) !< Column index
   real(kind=8), intent(out) :: val(nnz) !< Non-zero values

   integer :: i_row
   integer :: i_col
   integer :: i_val

   !$omp do
   do i_val = 1,nnz
      i_row = row_ind(i_val)
      i_col = col_ind(i_col)
      val(i_val) = matrix(i_row,i_col)
   enddo
   !$omp end do

end subroutine

subroutine ccs_to_dns(n_rows,n_cols,nnz,val,row_ind,col_ptr,matrix)

   implicit none

   integer,      intent(in)  :: n_rows
   integer,      intent(in)  :: n_cols
   integer,      intent(in)  :: nnz
   real(kind=8), intent(in)  :: val(nnz) !< Non-zero values
   integer,      intent(in)  :: row_ind(nnz) !< Row index
   integer,      intent(in)  :: col_ptr(n_cols+1) !< Column pointer
   real(kind=8), intent(out) :: matrix(n_rows,n_cols) !< Dense values

   integer :: i_row
   integer :: i_col
   integer :: i_val

   matrix = 0.d0
   i_col = 0

   do i_val = 1,nnz
      if(i_val == col_ptr(i_col+1) .and. i_col /= n_cols) then
         i_col = i_col+1
      endif
      i_row = row_ind(i_val)
      matrix(i_row,i_col) = val(i_val)
   enddo

end subroutine

subroutine crs_to_dns(n_rows,n_cols,nnz,val,col_ind,row_ptr,matrix)

   implicit none

   integer,      intent(in)  :: n_rows
   integer,      intent(in)  :: n_cols
   integer,      intent(in)  :: nnz
   real(kind=8), intent(in)  :: val(nnz) !< Non-zero values
   integer,      intent(in)  :: col_ind(nnz) !< Column index
   integer,      intent(in)  :: row_ptr(n_rows+1) !< Row pointer
   real(kind=8), intent(out) :: matrix(n_rows,n_cols) !< Dense values

   integer :: i_row
   integer :: i_col
   integer :: i_val

   matrix = 0.d0
   i_row = 0

   do i_val = 1,nnz
      if(i_val == row_ptr(i_row+1) .and. i_row /= n_rows) then
         i_row = i_row+1
      endif
      i_col = col_ind(i_val)
      matrix(i_row,i_col) = val(i_val)
   enddo

end subroutine

subroutine coo_to_dns(n_rows,n_cols,nnz,val,row_ind,col_ind,matrix)

   implicit none

   integer,      intent(in)  :: n_rows
   integer,      intent(in)  :: n_cols
   integer,      intent(in)  :: nnz
   real(kind=8), intent(in)  :: val(nnz) !< Non-zero values
   integer,      intent(in)  :: row_ind(nnz) !< Row index
   integer,      intent(in)  :: col_ind(nnz) !< Column index
   real(kind=8), intent(out) :: matrix(n_rows,n_cols) !< Dense matrix

   integer :: i_row
   integer :: i_col
   integer :: i_val

   matrix = 0.d0

   !$omp do
   do i_val = 1,nnz
      i_row = row_ind(i_val)
      i_col = col_ind(i_val)
      matrix(i_row,i_col) = val(i_val)
   enddo
   !$omp end do

end subroutine

subroutine ccs_to_coo(n_rows,n_cols,nnz,val,row_ind,col_ptr,col_ind)

   implicit none

   integer,      intent(in)  :: n_rows !< Number of rows
   integer,      intent(in)  :: n_cols !< Number of columns
   integer,      intent(in)  :: nnz !< Number of non-zeros
   real(kind=8), intent(in)  :: val(nnz) !< Non-zero values
   integer,      intent(in)  :: row_ind(nnz) !< Row index
   integer,      intent(in)  :: col_ptr(n_cols+1) !< Column pointer
   integer,      intent(out) :: col_ind(nnz) !< Column index

   integer :: i_col
   integer :: offset
   integer :: col_len

   offset = 1

   do i_col = 1,n_cols
      col_len = col_ptr(i_col+1)-col_ptr(i_col)
      col_ind(offset:offset+col_len-1) = i_col
      offset = offset+col_len
   enddo

end subroutine

subroutine crs_to_coo(n_rows,n_cols,nnz,val,col_ind,row_ptr,row_ind)

   implicit none

   integer,      intent(in)  :: n_rows !< Number of rows
   integer,      intent(in)  :: n_cols !< Number of columns
   integer,      intent(in)  :: nnz !< Number of non-zeros
   real(kind=8), intent(in)  :: val(nnz) !< Non-zero values
   integer,      intent(in)  :: col_ind(nnz) !< Column index
   integer,      intent(in)  :: row_ptr(n_rows+1) !< Row pointer
   integer,      intent(out) :: row_ind(nnz) !< Row index

   integer :: i_row
   integer :: offset
   integer :: row_len

   offset = 1

   do i_row = 1,n_rows
      row_len = row_ptr(i_row+1)-row_ptr(i_row)
      row_ind(offset:offset+row_len-1) = i_row
      offset = offset+row_len
   enddo

end subroutine

subroutine coo_to_ccs(n_rows,n_cols,nnz,val,row_ind,col_ind,col_ptr)

   implicit none

   integer,      intent(in)    :: n_rows !< Number of rows
   integer,      intent(in)    :: n_cols !< Number of columns
   integer,      intent(in)    :: nnz !< Number of non-zeros
   real(kind=8), intent(inout) :: val(nnz) !< Non-zero values
   integer,      intent(inout) :: row_ind(nnz) !< Row index
   integer,      intent(inout) :: col_ind(nnz) !< Column index
   integer,      intent(out)   :: col_ptr(n_cols+1) !< Column pointer

   real(kind=8) :: tmp_real
   integer      :: tmp_int
   integer      :: i_val
   integer      :: i_col
   integer      :: min_id

   integer, dimension(:), allocatable :: aux_ind

   allocate(aux_ind(nnz))

   do i_val = 1,nnz
      aux_ind(i_val) = row_ind(i_val)+(col_ind(i_val)-1)*n_rows
   enddo

   do i_val = 1,nnz
      min_id = minloc(aux_ind(i_val:nnz),1)+i_val-1

      tmp_int = aux_ind(i_val)
      aux_ind(i_val) = aux_ind(min_id)
      aux_ind(min_id) = tmp_int

      tmp_int = col_ind(i_val)
      col_ind(i_val) = col_ind(min_id)
      col_ind(min_id) = tmp_int

      tmp_int = row_ind(i_val)
      row_ind(i_val) = row_ind(min_id)
      row_ind(min_id) = tmp_int

      tmp_real = val(i_val)
      val(i_val) = val(min_id)
      val(min_id) = tmp_real
   enddo

   deallocate(aux_ind)

   allocate(aux_ind(n_cols))
   aux_ind = 0

   do i_val = 1,nnz
      aux_ind(col_ind(i_val)) = aux_ind(col_ind(i_val))+1
   enddo

   col_ptr(1) = 1

   do i_col = 1,n_cols
      col_ptr(i_col+1) = col_ptr(i_col)+aux_ind(i_col)
   enddo

end subroutine

subroutine coo_to_crs(n_rows,n_cols,nnz,val,row_ind,col_ind,row_ptr)

   implicit none

   integer,      intent(in)    :: n_rows !< Number of rows
   integer,      intent(in)    :: n_cols !< Number of columns
   integer,      intent(in)    :: nnz !< Number of non-zeros
   real(kind=8), intent(inout) :: val(nnz) !< Non-zero values
   integer,      intent(inout) :: row_ind(nnz) !< Row index
   integer,      intent(inout) :: col_ind(nnz) !< Column index
   integer,      intent(out)   :: row_ptr(n_rows+1) !< Row pointer

   real(kind=8) :: tmp_real
   integer      :: tmp_int
   integer      :: i_val
   integer      :: i_row
   integer      :: min_id

   integer, dimension(:), allocatable :: aux_ind

   allocate(aux_ind(nnz))

   do i_val = 1,nnz
      aux_ind(i_val) = col_ind(i_val)+(row_ind(i_val)-1)*n_cols
   enddo

   do i_val = 1,nnz
      min_id = minloc(aux_ind(i_val:nnz),1)+i_val-1

      tmp_int = aux_ind(i_val)
      aux_ind(i_val) = aux_ind(min_id)
      aux_ind(min_id) = tmp_int

      tmp_int = col_ind(i_val)
      col_ind(i_val) = col_ind(min_id)
      col_ind(min_id) = tmp_int

      tmp_int = row_ind(i_val)
      row_ind(i_val) = row_ind(min_id)
      row_ind(min_id) = tmp_int

      tmp_real = val(i_val)
      val(i_val) = val(min_id)
      val(min_id) = tmp_real
   enddo

   deallocate(aux_ind)

   allocate(aux_ind(n_rows))
   aux_ind = 0

   do i_val = 1,nnz
      aux_ind(row_ind(i_val)) = aux_ind(row_ind(i_val))+1
   enddo

   row_ptr(1) = 1

   do i_row = 1,n_rows
      row_ptr(i_row+1) = row_ptr(i_row)+aux_ind(i_row)
   enddo

end subroutine

subroutine ccs_to_crs(n_rows,n_cols,nnz,val,row_ind,col_ptr,col_ind,row_ptr)

   implicit none

   integer,      intent(in)    :: nnz !< Number of non-zeros
   integer,      intent(in)    :: n_rows !< Number of rows
   integer,      intent(in)    :: n_cols !< Number of columns
   real(kind=8), intent(inout) :: val(nnz) !< Non-zero values
   integer,      intent(in)    :: row_ind(nnz) !< Row index
   integer,      intent(in)    :: col_ptr(n_cols+1) !< Column pointer
   integer,      intent(out)   :: col_ind(nnz) !< Column index
   integer,      intent(out)   :: row_ptr(n_rows+1) !< Row pointer

   integer, dimension(:), allocatable :: aux_ind

   call ccs_to_coo(n_rows,n_cols,nnz,val,row_ind,col_ptr,col_ind)

   allocate(aux_ind(nnz))

   aux_ind = row_ind

   call coo_to_crs(n_rows,n_cols,nnz,val,aux_ind,col_ind,row_ptr)

   deallocate(aux_ind)

end subroutine

subroutine crs_to_ccs(n_rows,n_cols,nnz,val,col_ind,row_ptr,row_ind,col_ptr)

   implicit none

   integer,      intent(in)    :: nnz !< Number of non-zeros
   integer,      intent(in)    :: n_rows !< Number of rows
   integer,      intent(in)    :: n_cols !< Number of columns
   real(kind=8), intent(inout) :: val(nnz) !< Non-zero values
   integer,      intent(in)   :: col_ind(nnz) !< Column index
   integer,      intent(in)   :: row_ptr(n_rows+1) !< Row pointer
   integer,      intent(out)    :: row_ind(nnz) !< Row index
   integer,      intent(out)    :: col_ptr(n_cols+1) !< Column pointer

   integer, dimension(:), allocatable :: aux_ind

   call crs_to_coo(n_rows,n_cols,nnz,val,col_ind,row_ptr,row_ind)

   allocate(aux_ind(nnz))

   aux_ind = col_ind

   call coo_to_ccs(n_rows,n_cols,nnz,val,aux_ind,row_ind,col_ptr)

   deallocate(aux_ind)

end subroutine

end module

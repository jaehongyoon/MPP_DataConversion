module mat_redist

   use mat_util, only: get_nnz,get_global_idx

   implicit none

   integer :: l_nnz
   integer :: l_nnz2
   integer :: l_cols2

   real(kind=8), dimension(:), allocatable :: val
   integer,      dimension(:), allocatable :: row_idx
   integer,      dimension(:), allocatable :: col_ptr

contains

subroutine redist_2dbc_to_1db(threshold,global_size,l_rows,l_cols,matrix,&
                              block_size,mpi_comm,blacs_ctxt)

   implicit none

   include "mpif.h"

   real(kind=8), intent(in) :: threshold
   integer,      intent(in) :: global_size
   integer,      intent(in) :: l_rows
   integer,      intent(in) :: l_cols
   real(kind=8), intent(in) :: matrix(l_rows,l_cols)
   integer,      intent(in) :: block_size
   integer,      intent(in) :: mpi_comm
   integer,      intent(in) :: blacs_ctxt

   integer :: mpi_id         !< MPI task id
   integer :: mpi_size       !< Size of MPI communicator
   integer :: mpierr
   integer :: pid_row
   integer :: pid_col
   integer :: psize_row
   integer :: psize_col
   integer :: i_row          !< Row counter
   integer :: i_col          !< Col counter
   integer :: i_val          !< Value counter
   integer :: j_val          !< Value counter
   integer :: i_proc         !< Process counter
   integer :: global_col_idx !< Global column id
   integer :: global_row_idx !< Global row id
   integer :: local_col_idx  !< Local column id in 1D block distribution
   integer :: local_row_idx  !< Local row id in 1D block distribution
   integer :: tmp_int
   integer :: min_pos
   integer :: min_idx
   real(kind=8) :: tmp_real
   integer, dimension(:), allocatable :: dest !< Destination of each element

   ! See documentation of MPI_Alltoallv
   real(kind=8), dimension(:), allocatable :: val_send_buffer !< Send buffer for Hamiltonian
   integer,      dimension(:), allocatable :: idx_send_buffer !< Send buffer for global 1D id
   integer,      dimension(:), allocatable :: send_count      !< Number of elements to send to each processor
   integer,      dimension(:), allocatable :: send_displ      !< Displacement from which to take the outgoing data
   real(kind=8), dimension(:), allocatable :: val_recv_buffer !< Receive buffer for Hamiltonian
   integer,      dimension(:), allocatable :: idx_recv_buffer !< Receive buffer for global 1D id
   integer,      dimension(:), allocatable :: recv_count      !< Number of elements to receive from each processor
   integer,      dimension(:), allocatable :: recv_displ      !< Displacement at which to place the incoming data

   integer :: send_displ_aux !< Auxiliary variable used to set displacement 
   integer :: recv_displ_aux !< Auxiliary variable used to set displacement

   call MPI_Comm_rank(mpi_comm,mpi_id,mpierr)
   call MPI_Comm_size(mpi_comm,mpi_size,mpierr)

   call blacs_gridinfo(blacs_ctxt,psize_row,psize_col,pid_row,pid_col)

   call get_nnz(threshold,l_rows,l_cols,matrix,l_nnz)

   allocate(dest(l_nnz))
   allocate(idx_send_buffer(l_nnz))
   allocate(val_send_buffer(l_nnz))
   dest = 0
   idx_send_buffer = 0
   val_send_buffer = 0.d0

   l_cols2 = global_size/mpi_size
   if(mpi_id == 0) then
      l_cols2 = global_size-(mpi_size-1)*l_cols2
   endif

   ! Compute destination and global 1D id
   i_val = 0
   do i_col = 1,l_cols
      do i_row = 1,l_rows
         if(abs(matrix(i_row,i_col)) > threshold) then
            i_val = i_val+1
            call get_global_idx(global_row_idx,i_row,block_size,&
                                pid_row,psize_row)
            call get_global_idx(global_col_idx,i_col,block_size,&
                                pid_col,psize_col)

            ! Compute destination
            dest(i_val) = (global_col_idx-1)/(global_size/mpi_size)
            ! The last process may take more
            if(dest(i_val) > (mpi_size-1)) dest(i_val) = mpi_size-1

            ! Compute the global id
            ! Pack global id and data into buffers
            idx_send_buffer(i_val) = (global_col_idx-1)*global_size+global_row_idx
            val_send_buffer(i_val) = matrix(i_row,i_col)
        endif
     enddo
   enddo

   allocate(send_count(mpi_size))
   allocate(recv_count(mpi_size))
   send_count = 0
   recv_count = 0

   ! Set send_count
   do i_proc = 1,mpi_size
      do i_val = 1,l_nnz
         if(dest(i_val) == i_proc-1) then
            send_count(i_proc) = send_count(i_proc)+1
         endif
      enddo
   enddo

   deallocate(dest)

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer,recv_count,&
                     1,mpi_integer,mpi_comm,mpierr)

   ! Set local number of nonzero
   l_nnz2 = sum(recv_count,1)

   allocate(send_displ(mpi_size))
   allocate(recv_displ(mpi_size))
   send_displ = 0
   recv_displ = 0

   ! Set send and receive displacement
   send_displ_aux = 0
   recv_displ_aux = 0
   do i_proc = 1,mpi_size
      send_displ(i_proc) = send_displ_aux
      send_displ_aux = send_displ_aux+send_count(i_proc)
      recv_displ(i_proc) = recv_displ_aux
      recv_displ_aux = recv_displ_aux+recv_count(i_proc)
   enddo

   ! Send and receive the packed data
   allocate(idx_recv_buffer(l_nnz2))
   idx_recv_buffer = 0

   call MPI_Alltoallv(idx_send_buffer,send_count,send_displ,mpi_integer,&
                      idx_recv_buffer,recv_count,recv_displ,mpi_integer,&
                      mpi_comm,mpierr)

   deallocate(idx_send_buffer)

   allocate(val_recv_buffer(l_nnz2))
   val_recv_buffer = 0.d0

   call MPI_Alltoallv(val_send_buffer,send_count,send_displ,mpi_real8,&
                      val_recv_buffer,recv_count,recv_displ,mpi_real8,&
                      mpi_comm,mpierr)

   deallocate(send_count)
   deallocate(recv_count)
   deallocate(send_displ)
   deallocate(recv_displ)

   ! Unpack and reorder
   do i_val = 1,l_nnz2
      min_idx = minloc(idx_recv_buffer(i_val:l_nnz2),1)+i_val-1

      tmp_int = idx_recv_buffer(i_val)
      idx_recv_buffer(i_val) = idx_recv_buffer(min_idx)
      idx_recv_buffer(min_idx) = tmp_int

      tmp_real = val_recv_buffer(i_val)
      val_recv_buffer(i_val) = val_recv_buffer(min_idx)
      val_recv_buffer(min_idx) = tmp_real
   enddo

   allocate(val(l_nnz2))
   val = 0.d0

   val = val_recv_buffer

   deallocate(val_recv_buffer)

   allocate(row_idx(l_nnz2))
   allocate(col_ptr(l_cols2+1))
   row_idx = 0
   col_ptr = 0

   ! Compute row index and column pointer
   i_col = (idx_recv_buffer(1)-1)/global_size
   do i_val = 1,l_nnz2
      row_idx(i_val) = mod(idx_recv_buffer(i_val)-1,global_size)+1
      if((idx_recv_buffer(i_val)-1)/global_size+1 > i_col) then
         i_col = i_col+1
         col_ptr(i_col-(idx_recv_buffer(1)-1)/global_size) = i_val
      endif
   enddo

   col_ptr(l_cols2+1) = l_nnz2+1

   deallocate(idx_recv_buffer)

end subroutine

subroutine redist_1db_to_2dbc(global_size,l_rows,l_cols,matrix,&
                              block_size,mpi_comm,blacs_ctxt)

   implicit none

   include "mpif.h"

   integer,      intent(in)  :: global_size
   integer,      intent(in)  :: l_rows
   integer,      intent(in)  :: l_cols
   real(kind=8), intent(out) :: matrix(l_rows,l_cols)
   integer,      intent(in)  :: block_size
   integer,      intent(in)  :: mpi_comm
   integer,      intent(in)  :: blacs_ctxt

   integer :: mpi_id         !< MPI task id
   integer :: mpi_size       !< Size of MPI communicator
   integer :: mpierr
   integer :: pid_row
   integer :: pid_col
   integer :: psize_row
   integer :: psize_col
   integer :: i_row          !< Row counter
   integer :: i_col          !< Col counter
   integer :: i_val          !< Value counter
   integer :: j_val          !< Value counter
   integer :: k_val          !< Value counter
   integer :: i_proc         !< Process counter
   integer :: global_col_idx !< Global column id
   integer :: global_row_idx !< Global row id
   integer :: local_col_idx  !< Local column id in 1D block distribution
   integer :: local_row_idx  !< Local row id in 1D block distribution
   integer :: proc_col_idx   !< Column id in process grid
   integer :: proc_row_idx   !< Row id in process grid

   integer, dimension(:), allocatable :: dest       !< Destination of each element
   integer, dimension(:), allocatable :: global_idx !< Global 1d id

   ! See documentation of MPI_Alltoallv
   real(kind=8), dimension(:), allocatable :: val_send_buffer !< Send buffer for Hamiltonian
   integer,      dimension(:), allocatable :: idx_send_buffer !< Send buffer for global 1D id
   integer,      dimension(:), allocatable :: send_count      !< Number of elements to send to each processor
   integer,      dimension(:), allocatable :: send_displ      !< Displacement from which to take the outgoing data
   real(kind=8), dimension(:), allocatable :: val_recv_buffer !< Receive buffer for Hamiltonian
   integer,      dimension(:), allocatable :: idx_recv_buffer !< Receive buffer for global 1D id
   integer,      dimension(:), allocatable :: recv_count      !< Number of elements to receive from each processor
   integer,      dimension(:), allocatable :: recv_displ      !< Displacement at which to place the incoming data

   integer :: send_displ_aux !< Auxiliary variable used to set displacement 
   integer :: recv_displ_aux !< Auxiliary variable used to set displacement

   call MPI_Comm_rank(mpi_comm,mpi_id,mpierr)
   call MPI_Comm_size(mpi_comm,mpi_size,mpierr)

   call blacs_gridinfo(blacs_ctxt,psize_row,psize_col,pid_row,pid_col)

   allocate(val_send_buffer(l_nnz2))
   allocate(idx_send_buffer(l_nnz2))
   allocate(send_count(mpi_size))
   allocate(global_idx(l_nnz2))
   allocate(dest(l_nnz2))
   val_send_buffer = 0.d0
   idx_send_buffer = 0
   send_count = 0
   global_idx = 0
   dest = 0

   i_col = 0
   ! Compute destination and global 1D id
   do i_val = 1,l_nnz2
      if(i_val == col_ptr(i_col+1) .and. i_col /= l_cols2) then
         i_col = i_col+1
      endif
      i_row = row_idx(i_val)

      ! Compute global id
      global_row_idx = i_row
      global_col_idx = i_col+mpi_id*(global_size/mpi_size)
      global_idx(i_val) = (global_col_idx-1)*global_size+global_row_idx

      ! Compute destination
      proc_row_idx = mod((global_row_idx-1)/block_size,psize_row)
      proc_col_idx = mod((global_col_idx-1)/block_size,psize_col)
      dest(i_val)  = proc_col_idx+proc_row_idx*psize_col
   enddo

   j_val = 0
   k_val = l_nnz2+1

   ! Set send_count
   do i_proc = 1,mpi_size/2
      do i_val = 1,l_nnz2
         if(dest(i_val) == i_proc-1) then
            j_val = j_val+1
            val_send_buffer(j_val) = val(i_val)
            idx_send_buffer(j_val) = global_idx(i_val)
            send_count(i_proc) = send_count(i_proc)+1
         endif
         if(dest(i_val) == mpi_size-i_proc) then
            k_val = k_val-1
            val_send_buffer(k_val) = val(i_val)
            idx_send_buffer(k_val) = global_idx(i_val)
            send_count(mpi_size+1-i_proc) = send_count(mpi_size+1-i_proc)+1
         endif
      enddo
   enddo

   deallocate(global_idx)
   deallocate(dest)

   allocate(recv_count(mpi_size))
   recv_count = 0

   ! Set recv_count
   call MPI_Alltoall(send_count,1,mpi_integer,recv_count,&
                     1,mpi_integer,mpi_comm,mpierr)

   l_nnz = sum(recv_count,1)

   ! Set send and receive displacement
   allocate(send_displ(mpi_size))
   allocate(recv_displ(mpi_size))
   send_displ = 0
   recv_displ = 0

   send_displ_aux = 0
   recv_displ_aux = 0

   do i_proc = 1,mpi_size
      send_displ(i_proc) = send_displ_aux
      send_displ_aux = send_displ_aux+send_count(i_proc)

      recv_displ(i_proc) = recv_displ_aux
      recv_displ_aux = recv_displ_aux+recv_count(i_proc)
   enddo

   ! Send and receive the packed data
   allocate(val_recv_buffer(l_nnz))
   val_recv_buffer = 0.d0

   call MPI_Alltoallv(val_send_buffer,send_count,send_displ,mpi_real8,&
                      val_recv_buffer,recv_count,recv_displ,mpi_real8,&
                      mpi_comm,mpierr)

   deallocate(val_send_buffer)

   allocate(idx_recv_buffer(l_nnz))
   idx_recv_buffer = 0

   call MPI_Alltoallv(idx_send_buffer,send_count,send_displ,mpi_integer,&
                      idx_recv_buffer,recv_count,recv_displ,mpi_integer,&
                      mpi_comm,mpierr)

   deallocate(idx_send_buffer)
   deallocate(send_count)
   deallocate(recv_count)
   deallocate(send_displ)
   deallocate(recv_displ)

   matrix = 0.d0

   ! Unpack density matrix
   do i_val = 1,l_nnz
      ! Compute global 2d id
      global_col_idx = (idx_recv_buffer(i_val)-1)/global_size+1
      global_row_idx = mod(idx_recv_buffer(i_val)-1,global_size)+1

      ! Compute local 2d id
      local_row_idx = (global_row_idx-1)/(psize_row*block_size)*block_size&
                      +mod((global_row_idx-1),block_size)+1
      local_col_idx = (global_col_idx-1)/(psize_col*block_size)*block_size&
                      +mod((global_col_idx-1),block_size)+1

      ! Put value to correct position
      matrix(local_row_idx,local_col_idx) = val_recv_buffer(i_val)
   enddo

   deallocate(val_recv_buffer)
   deallocate(idx_recv_buffer)
   deallocate(val)
   deallocate(row_idx)
   deallocate(col_ptr)

end subroutine

end module

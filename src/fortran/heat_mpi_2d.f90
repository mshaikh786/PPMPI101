
#define N           27
#define MAX_ITER    4000
#define TOL         1e-4
#define MAX_TEMP    100.0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Creating a module for globalizing the allocatable arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module grid
    real, dimension(:,:), allocatable ::T_old,T_new
end module grid


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Creating a module for globalizing mpi related data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module my_mpi_vars
    use mpi

  ! variables to store virtual topology information
    ! Global rank and number of MPI processes
    integer ::  grank, gprocs
    integer ::  globalN, localN
    ! Local rank and number of MPI processes
    integer lrank, lprocs
    integer :: ndims = 2                                    ! Since we are decomposing grid in 2D
    integer, dimension(2) :: dims                           ! placeholder for dimensions
    integer, dimension(2) :: mycoords                       ! placeholder for Cartesian coordinates
    integer, dimension(2) :: period                         ! placeholder for periodicity information
    integer :: reorder = 0

    integer :: x_axis = 1, y_axis = 0                       ! convinience variables for referencing axes
    integer :: right, left, up, down                        ! placeholder for ranks of nearest neighbours

    integer :: cartcomm                                     ! Declaration of new communicator
    integer :: ierr                                         ! error handler for MPI calls

    integer, dimension(8)  :: req                           ! placeholder for MPI Isend and Irecv requests
    integer, dimension(MPI_STATUS_SIZE,8)  :: status        ! placeholder for maintaining status of the MPI requests. Size defined in mpif.h
end module my_mpi_vars



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program heat_mpi_2d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use grid
    use mpi
    use my_mpi_vars

    implicit none


    ! Local rows and columns along with ghost region
    integer rows, cols;


    integer :: i,j,iter,rem
    real :: dT,dT_local



    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD,grank,ierr)
    call MPI_Comm_size(MPI_COMM_WORLD,gprocs,ierr)

    ! Create a Cartesian grid communicator
    call MPI_Dims_create(gprocs, ndims, dims, ierr)

    ! Let's create a new communicator with these dimensions
    period(1) = 0
    period(2) = 0
    call MPI_Cart_create (MPI_COMM_WORLD, ndims, dims, period, reorder,cartcomm, ierr)

    ! Know who my neighbors are
    call MPI_Comm_rank(cartcomm, lrank,ierr)                     ! Who am I in new communicator
    call MPI_Cart_coords(cartcomm, lrank, ndims, mycoords,ierr)  !Get my Cartesian coordinates
    call MPI_Comm_size(cartcomm, lprocs,ierr)

    call MPI_Cart_shift(cartcomm, x_axis, 1, left, right,ierr);
    call MPI_Cart_shift(cartcomm, y_axis, 1, up, down,ierr);

    ! decide my local N
    globalN = N
    localN =  globalN / dims(1)
    rem = mod(globalN,gprocs)
    if (grank == gprocs - 1) then
        localN = localN + rem
    endif
    rows = localN + 2
    cols = localN + 2

    ! allocate memory
    allocate(T_old(rows,cols))
    allocate(T_new(rows,cols))

    call initialize(rows,cols)
    iter = 0
    dT = MAX_TEMP

    do while ( (dT > TOL) .AND. (iter <= MAX_ITER) )
        ! update boundaries
        call halo_update(rows,cols)
        do i = 2,rows-1
            do j = 2,cols-1
               T_new(i,j) = 0.25 * ( T_old(i-1,j) + T_old(i+1,j) + T_old(i,j-1) + T_old(i,j+1) )
            end do
        end do

        ! Find the maximum error on each rank and let every rank know what is the maximum on any rank
        dT = 0.0
        dT_local = 0.0
        do i = 2,rows-1
            do j = 2,cols-1
               dT_local = max(abs(T_new(i,j) - T_old(i,j)) , dT)
               T_old(i,j) = T_new(i,j)
            end do
        end do
        call MPI_Allreduce(dT_local, dT, 1, MPI_REAL, MPI_MAX, cartcomm,ierr);
        iter=iter+1
    end do

    !writing the output to the file is optional and is there for validation.
    call write_grid(rows,cols,iter)
    call MPI_Barrier(cartcomm,ierr)
    !call print_grid(rows,cols,iter)

    if (lrank == 0) then
        if ((iter - 1) == MAX_ITER) then
            print *, "Reached maximum iterations ",iter,". Error = ", dT
        else
            print *, lrank,": Converged in ",iter," iterations with and error of ",dT
        end if
    end if

call MPI_FINALIZE(ierr)
end program heat_mpi_2d


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialize the grids with appropirate boundary and initial conditions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine initialize(rows,cols)
    use grid
    use my_mpi_vars
    implicit none
    integer :: rows,cols
    integer :: i,j

    do i=1,rows
        do j=1,cols
                T_old(i,j) = 0.75 * MAX_TEMP
                T_new(i,j) = 0.75 * MAX_TEMP
        end do
    end do

    if (mycoords(2) == 0) then
        do i=1,rows
            T_old(i,1) = MAX_TEMP
            T_new(i,1) = MAX_TEMP
        end do
    end if
    if (mycoords(1) == dims(1)-1) then
        do j=1,cols
            T_old(rows,j) = MAX_TEMP
            T_new(rows,j) = MAX_TEMP
        end do
    end if
end subroutine initialize



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine for Halo exchange
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine halo_update(rows,cols)
    use grid
    use my_mpi_vars
    implicit none
    integer :: rows,cols
    integer :: i,j
    integer :: tag =100         ! not important for us but MPI expects a value.

    ! We will need to make send and receive buffers to exchange left and right columns
    real, dimension(rows) :: sendbuf_l
    real, dimension(rows) :: sendbuf_r
    real, dimension(cols) :: sendbuf_u
    real, dimension(cols) :: sendbuf_d

    real, dimension(rows) :: recvbuf_l
    real, dimension(rows) :: recvbuf_r
    real, dimension(cols) :: recvbuf_u
    real, dimension(cols) :: recvbuf_d


    do i =1,rows
        sendbuf_l(i) = T_old(i,2)
        sendbuf_r(i) = T_old(i,cols-1)
        ! The recvbuf value will change only for those messages with destination rank > -1
        recvbuf_l(i) = T_old(i,1)
        recvbuf_r(i) = T_old(i,cols)
    end do

    do j = 1,cols
        sendbuf_u(j) = T_old(2,j)
        sendbuf_d(j) = T_old(rows-1,j)
        ! The recvbuf value will change only for those messages with destination rank > -1
        recvbuf_u(j) = T_old(1,j)
        recvbuf_d(j) = T_old(rows,j)
    end do

    call MPI_Irecv(recvbuf_u, cols, MPI_REAL, up, tag, cartcomm, req(5),ierr);
    call MPI_Irecv(recvbuf_d, cols, MPI_REAL, down, tag, cartcomm, req(6),ierr);
    call MPI_Irecv(recvbuf_l, rows, MPI_REAL, left, tag, cartcomm, req(7),ierr);
    call MPI_Irecv(recvbuf_r, rows, MPI_REAL, right, tag, cartcomm, req(8),ierr);

    call MPI_Isend(sendbuf_u, cols, MPI_REAL, up, tag, cartcomm, req(1),ierr)
    call MPI_Isend(sendbuf_d, cols, MPI_REAL, down, tag, cartcomm,req(2),ierr);
    call MPI_Isend(sendbuf_l, rows, MPI_REAL, left, tag, cartcomm, req(3),ierr);
    call MPI_Isend(sendbuf_r, rows, MPI_REAL, right, tag, cartcomm,req(4),ierr);

    call MPI_Waitall(8, req, status,ierr);

    !Copy the received columns from left and right neighbors into the ghost region
    do i=1,rows
        T_old(i,1) = recvbuf_l(i)
        T_old(i,cols) = recvbuf_r(i)
    end do

    do j=1,cols
        T_old(1,j) = recvbuf_u(j)
        T_old(rows,j) = recvbuf_d(j)
    end do
 end subroutine halo_update



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! write grid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_grid(rows,cols,iter)
    use grid
    use my_mpi_vars
    integer :: i,j,iter,rows,cols,root
    integer :: x,y,indx                            ! some scratch variables
    integer :: numelems_local,recv_buf_size
    ! following allocations are significant on to root rank
    integer :: recv_nelems(lprocs)                 ! number of elements in a vector to receive from each rank
    integer :: disp(lprocs)                        ! array to hold offsets where data from each rank will go in the receive buffer

    real, dimension(:), allocatable :: sendbuf,recvbuf      ! send buffer to send the data to root and recvbuffer is only signifcant to root
    character(len=256)   ::  fname                         ! output filename

    root=0
    numelems_local = localN * localN
    ! tell Writer how many elements each rank is going to send
    call MPI_Gather(numelems_local, 1, MPI_INT, recv_nelems, 1, MPI_INT, root,cartcomm,ierr)

    ! disp is an array of offsets or displacements
    ! i.e. start index where to receive vector from each rank (useful only on master)
    if (lrank == root) then
        disp(1) = 0
        do i = 2,lprocs
            disp(i) = disp(i-1)+ recv_nelems(i)
        end do
    end if

    ! Now prepare for gathering the actual data
    allocate(sendbuf(numelems_local))
    indx=1
    do i = 2, rows-1
       do j = 2, cols-1
            sendbuf(indx) = T_old(i,j)
            indx=indx+1
       end do
    end do
    recv_buf_size=0
    do i = 1,lprocs
        recv_buf_size = recv_buf_size + recv_nelems(i)
    end do

    if (lrank == root) then
        allocate(recvbuf(recv_buf_size))
    end if

    call MPI_Gatherv(sendbuf, numelems_local, MPI_REAL, recvbuf, recv_nelems,disp, MPI_REAL, root, cartcomm,ierr);

    if (lrank == root) then
        write(fname,'(A,I0,A)') 'grid_t',iter,'.txt'
        open (unit = 10, file = fname, status='replace')
        indx=0
        do y=1,recv_buf_size,dims(2)*numelems_local
            do i=1,numelems_local,localN
                do x=1,dims(1)*numelems_local,numelems_local
                    do j=1,localN
                        indx = (y-1) + (i-1) + (x-1) + j
                        write(10,fmt='((F10.5))',advance='no') recvbuf(indx)
                    end do
                end do
                write(10,*)
            end do
        end do
        close(10)
    end if
end subroutine write_grid






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! debugging printer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine print_grid(rows,cols,iter)
    use grid
    use my_mpi_vars
    implicit none
    integer :: rows,cols,iter
    integer :: i,j
    character(len=32) :: fname
    write(fname,'(A,I0,A,I0,A)') 'grid_',lrank,'_t',iter,'.txt'

    open (unit = 10, file = fname, status='replace')
    do i=2,rows-1
       write(10,*) (T_old(i,j), j=2,cols-1)
    end do
    close(10)
 end subroutine print_grid


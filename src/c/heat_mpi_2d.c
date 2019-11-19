/* 
 * 1D decomposition of the domain in the column dimension
 * such that each MPI task will get an equal number of rows.
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N			27
#define MAX_ITER	4000
#define TOL			1e-4
#define MAX_TEMP	100.0

float** allocate(float**);
void init_grid(float**);
void halo_update(float**);
void print_grid(int, float**);
// Global rank and number of MPI processes
int grank, gprocs;
int globalN, localN;
// Local rank and number of MPI processes
int lrank, lprocs;
// Local rows and columns along with ghost region
int rows, cols;

// variables to store virtual topology information
int ndims = 2;			// Since we are decomposing grid in 2D
int dims[2] = { 0, 0 };  	// placeholder for dimensions
int mycoords[2] = { 0, 0 };	// placeholder for Cartesian coordinates
int period[2] = { 0, 0 };	// placeholder for periodicity information
int reorder = 0;

int x_axis = 1, y_axis = 0;
int right, left, up, down;

MPI_Comm cartcomm;		// Declaration of new communicator

int main(int argc, char *argv[]) {

	int i, j;
	float **T_old, **T_new;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &grank);
	MPI_Comm_size(MPI_COMM_WORLD, &gprocs);

	// Create a Cartesian grid communicator

	MPI_Dims_create(gprocs, ndims, dims);	// Get a good guess on dimensions

	// Let's create a new communicator with these dimensions
	MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, period, reorder, &cartcomm);
	// Know who my neighbors are
	MPI_Comm_rank(cartcomm, &lrank);			// Who am I in new communicator
	MPI_Cart_coords(cartcomm, lrank, ndims, mycoords);//Get my Cartesian coordinates
	MPI_Comm_size(cartcomm, &lprocs);

	MPI_Cart_shift(cartcomm, x_axis, 1, &left, &right);
	MPI_Cart_shift(cartcomm, y_axis, 1, &up, &down);
	// decide my local N
	globalN = N;
	localN = (int) (globalN / dims[0]);
	int rem = globalN % gprocs;
	if (grank == gprocs - 1) {
		localN += rem;
	}
	rows = localN + 2;
	cols = localN + 2;

	T_old = allocate(T_old);
	T_new = allocate(T_new);

	init_grid(T_old);

	print_grid(0,T_old);


	int iter = 0;

	float dT=MAX_TEMP,
		  dT_local = MAX_TEMP;
	while (dT > TOL && iter <= MAX_ITER) {
		halo_update(T_old);
		for (i = 1; i < rows - 1; i++) {
			for (j = 1; j < cols - 1; j++) {
				T_new[i][j] = 0.25
						* (T_old[i - 1][j] + T_old[i + 1][j] + T_old[i][j - 1]
								+ T_old[i][j + 1]);
			}
		}
		dT = 0.0;
		dT_local = 0.0;
		for (i = 1; i < rows - 1; i++) {
			for (j = 1; j < cols - 1; j++) {
				dT_local = fmaxf(fabsf(T_new[i][j] - T_old[i][j]), dT_local);
				T_old[i][j] = T_new[i][j];
			}
		}

		MPI_Allreduce(&dT_local, &dT, 1, MPI_FLOAT, MPI_MAX, cartcomm);

		iter++;
	}

	// Printing results
	if (lrank == 0) {
		if ((iter - 1) == MAX_ITER)
			printf("Reached maximum iterations %d. Error = %2.4f\n",
					iter - 1, dT);
		else
			printf("Converged in %d iterations with and error of %2.4f\n",
					iter - 1, dT);
	}

	MPI_Finalize();
	return 0;
}

float** allocate(float **T) {
	int i = 0, j = 0;

	// Allocate including some extra for ghost regions

	T = (float**) malloc(rows * sizeof(float*));
	if (T != NULL) {
		T[0] = (float*) malloc(rows * cols * sizeof(float));
	}
	if (T[0] != NULL) {
		for (i = 0; i < rows; i++)
			T[i] = (*T + i * cols);
	}
	return T;
}

void init_grid(float **T) {
	int i, j;

	for (i = 0; i < rows; i++)
		for (j = 0; j < cols; j++)
			T[i][j] = 0.75 * MAX_TEMP;

	if (mycoords[0] == dims[0] - 1) {
		for (j = 0; j < cols; j++)
			T[rows - 1][j] = 0.75 * MAX_TEMP;
	}
	if (mycoords[1] == 0) {
		for (i = 0; i < rows; i++)
			T[i][0] = MAX_TEMP;
	}

}

void halo_update(float **T) {

	int i, j;
	// We will need to make send and receive buffers to exchange left and right columns
	float *sendbuf_l = (float*) malloc(rows * sizeof(float));
	float *recvbuf_l = (float*) malloc(rows * sizeof(float));
	float *sendbuf_r = (float*) malloc(rows * sizeof(float));
	float *recvbuf_r = (float*) malloc(rows * sizeof(float));

	for (i = 0; i < rows; i++) {
		sendbuf_l[i] = T[i][1];
		sendbuf_r[i] = T[i][cols - 2];
		recvbuf_l[i] = 0.0;
		recvbuf_r[i] = 0.0;
	}
	MPI_Request req[8];
	MPI_Status status[8];
	int tag = 100; // not important for us but MPI expects a value.
	MPI_Isend(&T[1][0], cols, MPI_FLOAT, up, tag, cartcomm, &req[0]);
	MPI_Isend(&T[rows - 2][0], cols, MPI_FLOAT, down, tag, cartcomm, &req[1]);
	MPI_Isend(sendbuf_l, rows, MPI_FLOAT, left, tag, cartcomm, &req[2]);
	MPI_Isend(sendbuf_r, rows, MPI_FLOAT, right, tag, cartcomm, &req[3]);

	MPI_Irecv(&T[0][0], cols, MPI_FLOAT, up, tag, cartcomm, &req[4]);
	MPI_Irecv(&T[rows - 1][0], cols, MPI_FLOAT, down, tag, cartcomm, &req[5]);
	MPI_Irecv(recvbuf_l, rows, MPI_FLOAT, left, tag, cartcomm, &req[6]);
	MPI_Irecv(recvbuf_r, rows, MPI_FLOAT, right, tag, cartcomm, &req[7]);

	MPI_Waitall(8, req, status);

	// Copy the received columns from left and right neighbors into the ghost region
	for (i = 0; i < rows; i++) {
		T[i][0] = recvbuf_l[i];
		T[i][cols - 1] = recvbuf_r[i];
	}
	free(sendbuf_l);
	free(sendbuf_r);
	free(recvbuf_l);
	free(recvbuf_r);
}

void print_grid(int iter, float **T) {
	int i, j;
	char fname[256];
	sprintf(fname, "grid%d_%d.txt", lrank, iter);
	FILE *fw = fopen(fname, "w+");
	fprintf(fw, "grank: %d, lrank: %d\n", grank, lrank);
	fprintf(fw, "r,l,u,d: %d,%d,%d,%d\n", right, left, up, down);
	fprintf(fw, "dims: %d,%d\n", dims[0], dims[1]);
	fprintf(fw, "Coords: %d,%d\n", mycoords[0], mycoords[1]);
	fprintf(fw, "globalN: %d, localN: %d, rows: %d, cols: %d\n", globalN,
			localN, rows, cols);
	for (i = 0; i < rows; i++) {
		for (j = 0; j < cols; j++) {
			fprintf(fw, "%2.2lf ", T[i][j]);
		}
		fprintf(fw, "\n");
	}
	fclose(fw);
}

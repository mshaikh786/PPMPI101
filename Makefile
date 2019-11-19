CC=cc

all:
	$(CC)  -g ./src/c/heat_mpi_2d.c -o heat_mpi_2d -lm

clean:
	rm heat_mpi_2d

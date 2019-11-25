FC=ftn
FFLAGS=-g -eZ

all: heat_mpi_2d
heat_mpi_2d: 
	$(FC) $(FFLAGS) -c heat_mpi_2d.f90 
	$(FC) heat_mpi_2d.o -o heat_mpi_2d
clean:
	rm heat_mpi_2d.o heat_mpi_2d *.mod *.i

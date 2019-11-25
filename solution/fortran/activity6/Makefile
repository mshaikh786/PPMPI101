FC=ftn
FFLAGS=-g -eZ

all: heat_serial heat_mpi_2d

heat_serial:
	$(FC) $(FFLAGS) -c heat_serial.f90 
	$(FC) heat_serial.o -o heat_serial

heat_mpi_2d: 
	$(FC) $(FFLAGS) -c heat_mpi_2d.f90 
	$(FC) heat_mpi_2d.o -o heat_mpi_2d
clean:
	rm *.o heat_serial heat_mpi_2d *.mod *.i

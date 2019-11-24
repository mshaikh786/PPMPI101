CC=cc
FC=ftn
all: heat_serial_c heat_mpi_2d_c heat_serial_f90

heat_serial_c:
	$(CC) -g ./src/c/heat_serial.c -c 
	$(CC) -g heat_serial.o -o heat_serial_c -lm

heat_mpi_2d_c: 
	$(CC)  -g ./src/c/heat_mpi_2d.c -c
	$(CC)  -g heat_mpi_2d.o -o heat_mpi_2d_c -lm


heat_serial_f90:
	$(FC) -g -eZ ./src/fortran/heat_serial.f90 -c 
	$(FC) -g heat_serial.o -o heat_serial_f90 -lm


clean:
	rm heat_serial_c heat_mpi_2d_c heat_serial_f90  *.o *.mod *.i

CC=cc
CFLAGS=-g

all: heat_serial heat_mpi_2d

heat_serial: 
	$(CC) $(CFLAGS) -c heat_serial.c 
	$(CC) heat_serial.o -o heat_serial

heat_mpi_2d: 
	$(CC) $(CFLAGS) -c heat_mpi_2d.c 
	$(CC) heat_mpi_2d.o -o heat_mpi_2d
clean:
	rm *.o heat_mpi_2d heat_serial

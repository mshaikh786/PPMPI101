CC=cc
CFLAGS=-g
LIBS=-lm
all: heat_serial heat_mpi

heat_serial: 
	$(CC) $(CFLAGS) -c src/heat_serial.c 
	$(CC) $(CFLAGS) heat_serial.o -o heat_serial $(LIBS)

heat_mpi: 
	$(CC) $(CFLAGS) -c src/heat_mpi_2d.c 
	$(CC) $(CFLAGS) heat_mpi_2d.o -o heat_mpi_2d $(LIBS)
clean:
	rm *.o heat_mpi_2d heat_serial





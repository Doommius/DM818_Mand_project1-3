#
# Edison - NERSC 
#
# Intel Compilers are loaded by default; for other compilers please check the module list
#
CC = g++
MPCC = mpic++
OPENMP = -fopenmp
CFLAGS = -O3
LIBS = -lm


TARGETS = mpi

all:	$(TARGETS)


mpi: mpi.o common.o
	$(MPCC) -o $@ $(MPILIBS) mpi.o common.o $(LIBS)

mpi.o: mpi.cpp common.h
	$(MPCC) -c $(CFLAGS) mpi.cpp
common.o: common.cpp common.h
	$(CC) -c $(CFLAGS) common.cpp

clean:
	rm -f *.o $(TARGETS) *.stdout *.txt

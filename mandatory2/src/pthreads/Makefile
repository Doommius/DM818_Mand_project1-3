#
# Franklin and Hopper
#
CC = CC
MPCC = CC
OPENMP = -mp
CFLAGS = -O3
LIBS =

#
# Bassi  - Being retired soon/ Don't have access
#
#CC = cc -+ -qsuppress=1500-036
#MPCC = mpcc
#OPENMP = -qsmp=omp
#CFLAGS = -O3
#LIBS = -lm

#
# Jacquard
#
#CC = pathCC
#MPCC = mpicxx 
#OPENMP = -mp
#LIBS = -lm
#CFLAGS = -O3

#
# DaVinci - Being retired soon/ Don't have access
#
#CC = g++
#MPCC = g++
#OPENMP = -fopenmp
#LIBS = -lm
#CFLAGS = -O3
#MPILIBS = -lmpi++ -lmpi

TARGETS = pthreads

all:	$(TARGETS)

pthreads: pthreads.o common.o
	$(CC) -o $@ $(LIBS) -lpthread pthreads.o common.o

pthreads.o: pthreads.cpp common.h
	$(CC) -c $(CFLAGS) pthreads.cpp

clean:
	rm -f *.o $(TARGETS)

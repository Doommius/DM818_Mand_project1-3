#
# on Franklin:
# type "module load acml" to use ACML (required)
# type "module swap PrgEnv-pgi PrgEnv-gnu" to use GNU compiler
# type "module swap PrgEnv-pgi PrgEnv-pathscale" to use PathScale compiler
#
# On other machines you might need to use other library than ACML
# this could be MKL on an Intel platform or ESSL on an IBM platform.
#
#define __AVX2__ 1
#CC = CC
CC = /opt/intel/bin/icc
#LIBS = -lacml
#LIBS = $(shell pkg-config cblas --libs)

#multi-thread cblas for IMADA i5 3.2 GHz machine
#use the following line only with a hand compiled atlas
#LIBS = -L/home/daniel/atlas/lib -lptcblas -latlas -lpthread

#otherwise use
#LIBS = -lcblas -latlas -lpthread

#single-thread cblas for IMADA i5 3.2 GHz machine
#LIBS = -L/home/daniel/atlas/lib -lcblas -latlas

#not optimized Ubuntu version
LIBS= -lcblas

CFLAGS = -O3

all:	benchmark-naive benchmark-blocked benchmark-blas benchmark

benchmark-naive: benchmark.o dgemm-naive.o
	$(CC) -o $@ benchmark.o dgemm-naive.o $(LIBS)
benchmark-blocked: benchmark.o dgemm-blocked.o
	$(CC) -o $@ benchmark.o dgemm-blocked.o $(LIBS)
benchmark-blas: benchmark.o dgemm-blas.o
	$(CC) -o $@ benchmark.o dgemm-blas.o $(LIBS)


benchmark-dgemm : benchmark.o dgemm.o
	$(CC) -o $@ benchmark.o dgemm.o $(LIBS)

test-naive : test.o dgemm-naive.o
	$(CC) -o $@ test.o dgemm-naive.o

test-dgemm : test.o dgemm.o
	$(CC) -o $@ test.o dgemm.o

test-blocked : test.o dgemm-blocked.o
	$(CC) -o $@ test.o dgemm-blocked.o


%.o: %.cpp
	$(CC) -c $(CFLAGS) $<

%.o: %.c
	$(CC) -c $(CFLAGS) $<

clean:
	rm -f benchmark-naive benchmark-blocked benchmark-blas benchmark-dgemm test-blocked test-dgemm *.o

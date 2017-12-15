# Makefile for DNS-algorithm skeleton
# Made for Franklin
# For other machines, replace CC = cc with
# mpicc or other appropriate compiler

CC = cc
CFLAGS = -O3

CFLAGS = -Wall $(OPT) -std=gnu99 -funroll-loops -funit-at-a-time -mfpmath=sse -fstrict-aliasing -ftree-vectorize -ffast-math
LDFLAGS = -Wall

LIBS = -lrt

dnsmat: main.o matrix.o comm.o
	$(CC) $(CFLAGS) -o dnsmat main.o matrix.o $(LIBS)

%.o: %.c
	$(CC) -c $(CFLAGS) $<

clean:
	rm -f dnsmat *.o

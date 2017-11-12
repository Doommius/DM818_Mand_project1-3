




#not optimized Ubuntu version
LIBS= -lcblas


benchmark-cuda: benchmark.o dgemm-cuda.o
	$(CC) -o -I/usr/local/cuda/include $@ benchmark.o dgemm-cuda.o $(LIBS)

cuda-alone: dgemm-cuda.o

	$(CC) -o -I/usr/local/cuda/include $@ dgemm-cuda.o $(LIBS)

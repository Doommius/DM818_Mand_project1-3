rm mpi.txt
./mpi -n 500 -no -s mpi.txt
./mpi -n 1000 -no -s mpi.txt
./mpi -n 2000 -no -s mpi.txt
./mpi -n 4000 -no -s openmp.txt
./mpi -n 8000 -no -s mpi.txt
./mpi -n 16000 -no -s mpi.txt
./mpi -n 32000 -no -s mpi.txt
./mpi -n 64000 -no -s mpi.txt
./mpi -n 128000 -no -s mpi.txt
./mpi -n 256000 -no -s mpi.txt
./autograder -v mpi -s mpi.txt
rm mpi.txt

mpirun -n 10 -hostfile imada-hostfile.alive mand2_mpi_opti -n 500 -no -s mpi.txt
mpirun -n 10 -hostfile imada-hostfile.alive mand2_mpi_opti -n 1000 -no -s mpi.txt
mpirun -n 10 -hostfile imada-hostfile.alive mand2_mpi_opti -n 2000 -no -s mpi.txt
mpirun -n 10 -hostfile imada-hostfile.alive mand2_mpi_opti -n 4000 -no -s openmp.txt
mpirun -n 10 -hostfile imada-hostfile.alive mand2_mpi_opti -n 8000 -no -s mpi.txt
mpirun -n 10 -hostfile imada-hostfile.alive mand2_mpi_opti -n 16000 -no -s mpi.txt
mpirun -n 10 -hostfile imada-hostfile.alive mand2_mpi_opti -n 32000 -no -s mpi.txt
mpirun -n 10 -hostfile imada-hostfile.alive mand2_mpi_opti -n 64000 -no -s mpi.txt
mpirun -n 10 -hostfile imada-hostfile.alive mand2_mpi_opti -n 128000 -no -s mpi.txt
mpirun -n 10 -hostfile imada-hostfile.alive mand2_mpi_opti -n 256000 -no -s mpi.txt
./autograder -v mpi -s mpi.txt
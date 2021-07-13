. ap.setenv
cmake .
make -j4
#mpirun -n 2 --oversubscribe ./main --nx 64 --ny 64 --nz 32 --bs 32 --layers 1
mpirun -n 4 --oversubscribe ./main --nx 32 --ny 32 --nz 32 --bs 8 --layers 1
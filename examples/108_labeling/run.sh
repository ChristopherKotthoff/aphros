. ap.setenv
cmake .
make -j4
mpirun -n 2 --oversubscribe ./main --nx 64 --ny 64 --nz 32 --bs 32 --layers 1
#mpirun -n 8 --oversubscribe ./main --nx 64 --ny 64 --nz 64 --bs 32 --layers 8
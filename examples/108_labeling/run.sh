. ap.setenv
cmake .
make -j4
#mpirun -n 2 --oversubscribe ./main --nx 64 --ny 64 --nz 32 --bs 32 --layers 1
mpirun -n 4 ./main --nx 50 --ny 50 --nz 40 --bs 10 --layers 16 --new_recolor 1

#. ap.setenv
#cmake .
#make -j4
bsub -n 2 mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 1 snake.conf
bsub -n 4 mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 2 snake.conf
bsub -n 6 mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 3 snake.conf
bsub -n 12 mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 4 snake.conf
bsub -n 24 -R fullnode mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 5 snake.conf
bsub -n 48 -R fullnode mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 6 snake.conf
bsub -n 96 -R fullnode mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 7 snake.conf
bsub -n 192 -R fullnode mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 8 snake.conf
bsub -n 384 -R fullnode mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 9 snake.conf
bsub -n 768 -R fullnode mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 10 snake.conf










bsub -n 2 mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 11 snake.conf
bsub -n 4 mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 12 snake.conf
bsub -n 6 mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 13 snake.conf
bsub -n 12 mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 14 snake.conf
bsub -n 24 -R fullnode mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 15 snake.conf
bsub -n 48 -R fullnode mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 16 snake.conf
bsub -n 96 -R fullnode mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 17 snake.conf
bsub -n 192 -R fullnode mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 18 snake.conf
bsub -n 384 -R fullnode mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 19 snake.conf
bsub -n 768 -R fullnode mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 20 snake.conf










bsub -n 2 mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 21 snake.conf
bsub -n 4 mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 22 snake.conf
bsub -n 6 mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 23 snake.conf
bsub -n 12 mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 24 snake.conf
bsub -n 24 -R fullnode mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 25 snake.conf
bsub -n 48 -R fullnode mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 26 snake.conf
bsub -n 96 -R fullnode mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 27 snake.conf
bsub -n 192 -R fullnode mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 28 snake.conf
bsub -n 384 -R fullnode mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 29 snake.conf
bsub -n 768 -R fullnode mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 30 snake.conf










bsub -n 2 mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 31 snake.conf
bsub -n 4 mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 32 snake.conf
bsub -n 6 mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 33 snake.conf
bsub -n 12 mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 34 snake.conf
bsub -n 24 -R fullnode mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 35 snake.conf
bsub -n 48 -R fullnode mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 36 snake.conf
bsub -n 96 -R fullnode mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 37 snake.conf
bsub -n 192 -R fullnode mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 38 snake.conf
bsub -n 384 -R fullnode mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 39 snake.conf
bsub -n 768 -R fullnode mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 40 snake.conf










bsub -n 2 mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 41 snake.conf
bsub -n 4 mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 42 snake.conf
bsub -n 6 mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 43 snake.conf
bsub -n 12 mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 44 snake.conf
bsub -n 24 -R fullnode mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 45 snake.conf
bsub -n 48 -R fullnode mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 46 snake.conf
bsub -n 96 -R fullnode mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 47 snake.conf
bsub -n 192 -R fullnode mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 48 snake.conf
bsub -n 384 -R fullnode mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 49 snake.conf
bsub -n 768 -R fullnode mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 50 snake.conf










bsub -n 2 mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 51 snake.conf
bsub -n 4 mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 52 snake.conf
bsub -n 6 mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 53 snake.conf
bsub -n 12 mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 54 snake.conf
bsub -n 24 -R fullnode mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 55 snake.conf
bsub -n 48 -R fullnode mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 56 snake.conf
bsub -n 96 -R fullnode mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 57 snake.conf
bsub -n 192 -R fullnode mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 58 snake.conf
bsub -n 384 -R fullnode mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 59 snake.conf
bsub -n 768 -R fullnode mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 60 snake.conf










bsub -n 2 mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 61 snake.conf
bsub -n 4 mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 62 snake.conf
bsub -n 6 mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 63 snake.conf
bsub -n 12 mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 64 snake.conf
bsub -n 24 -R fullnode mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 65 snake.conf
bsub -n 48 -R fullnode mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 66 snake.conf
bsub -n 96 -R fullnode mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 67 snake.conf
bsub -n 192 -R fullnode mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 68 snake.conf
bsub -n 384 -R fullnode mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 69 snake.conf
bsub -n 768 -R fullnode mpirun ./main --nx 128 --ny 128 --nz 192 --bs 16 --layers 1 --new_recolor 0 --file_prefix 70 snake.conf











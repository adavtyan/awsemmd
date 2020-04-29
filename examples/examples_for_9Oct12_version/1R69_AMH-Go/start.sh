mpirun -np 4 ./lmp_linux < h4_traj9.in

#bsub -n 4 -o h4.out -e err mpirun  /ifs1/home/adavtyan/h4/lmp_linux -in h4_traj9.in -log input.log &

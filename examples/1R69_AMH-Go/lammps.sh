#!/bin/bash

lammps_dir="/home/aram/lammps-18Feb11/src"
project_dir="/home/aram/Project_Files_2009.05.18/amnesiacmd"
simulation_dir=`pwd`

pr_file1="fix_backbone.cpp"
pr_file2="fix_backbone.h"
pr_file3="fix_go-model.cpp"
pr_file4="fix_go-model.h"
pr_file5="pair_excluded_volume.cpp"
pr_file6="pair_excluded_volume.h"
#pr_file7="style.h"
pr_file8="smart_matrix_lib.h"
pr_file9="fix_qbias.cpp"
pr_file10="fix_qbias.h"
pr_file11="fragment_memory.cpp"
pr_file12="fragment_memory.h"

lmp_name="lmp_serial"
in_file="1r69.in"

cp -v $project_dir/$pr_file1 $project_dir/$pr_file2 $project_dir/$pr_file3 $project_dir/$pr_file4 $project_dir/$pr_file5 $project_dir/$pr_file6 $project_dir/$pr_file8 $project_dir/$pr_file9 $project_dir/$pr_file10 $project_dir/$pr_file11 $project_dir/$pr_file12  $lammps_dir
pwd

cd  $lammps_dir
#make clean-all
make serial

cp -v $lammps_dir/$lmp_name $simulation_dir
pwd

cd $simulation_dir

for i in $* ; do
        case $i in
          run) ./$lmp_name < $in_file ;;
	   mpi_run) ./start.sh ;;
        esac
done

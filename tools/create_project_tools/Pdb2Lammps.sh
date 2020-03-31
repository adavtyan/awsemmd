#!/bin/bash

export s=$HOME/opt/script

if [[ ! $# -eq 2 ]]
then
	echo
	echo "> $0 pdb_id project_name"
	echo
	exit
fi

pdb_file=$1
output_file=$2

echo $pdb_file
echo $output_file

python $s/PDBToSequanceFile.py $pdb_file $output_file".se"
python $s/SequanceToZ-Matrix.py $output_file".se" $output_file".zm" -d
python $s/Z-MatrixToCoordinates.py $output_file".zm" $output_file".coord"
python $s/CoordinatesToWorkLammpsDataFile.py $output_file".coord" "data."$output_file -b

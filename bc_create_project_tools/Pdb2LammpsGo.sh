#!/bin/bash

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

python PDBToSequanceFile.py $pdb_file $output_file".se"
python SequanceToZ-Matrix.py $output_file".se" $output_file".zm" -d
python Z-MatrixToCoordinates.py $output_file".zm" $output_file".coord"
python CoordinatesToWorkLammpsDataFile.py $output_file".coord" "data."$output_file -b -go

python GetContactMapFromPDB.py $pdb_file fix_gomodel_coeff.data

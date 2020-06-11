#!/bin/bash

#########################################################################
# Author: Haiqing Zhao
# Created Time: 05.10.2017
# Description: this script generates needed input files for protein_DNA simulations by AWSEM_3SPN model.
# Including: three DNA_list files, a combined data file and a seperately given pair_coeff file.
# The default sets for DNA are the curved B-dna and with no ions.

# Modified by Hao Wu
# Updated: Oct 13 2017
#########################################################################

export awsem_dir=/homes/haowu/software/awsemmd-master-May2017/create_project_tools
export dna_dir=/homes/haowu/software/lammps-31Mar17_3spn2c/src/USER-3SPN2

if [ "$#" -ne 3 ]; then
    echo
    echo "$0 protein_pdb_id dna_pdb dna_seq"
    echo
    exit
fi

protein_pdb_id=$1
#protein_res_num=$2
dna_pdb=$2
dna_seq=$3

echo $protein_pdb_id
echo $dna_pdb

# Create protein input files from AWSEM
python $awsem_dir/PDBToCoordinates.py $protein_pdb_id protein.coord
python $awsem_dir/CoordinatesToWorkLammpsDataFile.py protein.coord data.protein -b

# Create ideal atomistic model from x3dna for sequence-based parameters in list files 
python $dna_dir/utils/make_bp_params.py $dna_seq
x3dna_utils cp_std BDNA
rebuild -atomic bp_step.par atomistic.pdb
python $dna_dir/utils/pdb2cg_dna.py atomistic.pdb

# If no explicit ions are needed, use DSIM_ICNF_noIon: delete atom type 15-18
$dna_dir/DSIM_ICNF_noIon/icnf.exe $dna_seq 1 1 . 0

# Replace atom coordinates with ideal atomistic model
$dna_dir/utils/replace_atoms.sh conf_lammps.in dna_conf.in bdna_curv_conf.in

# Make list files, DNA atom indices are after the number of protein atoms
#protein_atom_num=$(expr $protein_res_num \* 3)
# read atom number from datafile
atom_data_line=$(sed -n '3p' data.protein) 
protein_atom_num=$(echo $atom_data_line | awk '{print $1;}')
$dna_dir/utils/make_list_files_protDNA.py bdna_curv_conf.in $protein_atom_num

# Replace DNA coordinates with that in dna_pdb provided
python $dna_dir/utils/pdb2cg_dna.py $dna_pdb
$dna_dir/utils/replace_atoms.sh conf_lammps.in dna_conf.in data.dna

# Combine protein and dna files into one LAMMPS data file
python proteinDna_combine.py data.protein data.dna protein.seq

# Create protein beads group definition file (CA, CB and O)
awk '/group/{print}' protein.in > protein_group_def.lmp

# Create DNA phosphate beads group definition file (for protein-DNA energy output)
# python write_phosphate_group.py $protein_atom_num

# Create a DNA bond definition file based on PSF file
awk '/!NBONDS/{flag=1}/!NTHETA/{flag=0}flag' in00_cvmd.psf > dnaBondFromPsf.txt

# Create a CG DNA PDB file
dna_bp_num="$(head -1 $dna_seq)"
python write_cg_dna_pdb.py $dna_bp_num in00_conf.crd cg_dna.pdb

####TO CLEAN UP UNNECESSARY FILES
#rm -rf dna_conf.in dna_bonds.in dna_angles.in dna_dihedrals.in in00_cvmd.psf in00_conf.crd in00_conf.xyz
#rm -rf conf_lammps.in dna_conf.in tmp.start tmp.end bdna_curv_conf.in protein.in
#rm -rf ref_frames.dat protein.coord data.protein data.dna bp_step.par atomistic.pdb Atomic*.pdb

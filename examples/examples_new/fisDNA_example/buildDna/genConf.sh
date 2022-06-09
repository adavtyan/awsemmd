#########################################################################
# Author: Charlse.Zhang
# Created Time: Thu 03 Apr 2014 05:47:21 PM CDT
# File Name: genConf.sh
# Description: Build 3SPN model for the DNA sequence of interest and generate
#              input files for running simulation using LAMMPS
# Last revision: 18 May 2017 by Min-Yeh (Victor), Tsai
#########################################################################
#!/bin/bash

spn2cdir="[PATH TO LOCAL FOLDER]/buildDna/USER-3SPN2" # setup 3SPN path

##### Generate All-Atom topology #####
python ${spn2cdir}/utils/make_bp_params.py dnaSeq.txt # output: bp_step.par
x3dna_utils cp_std BDNA # output: Atomic*.pdb
rebuild -atomic bp_step.par atomistic.pdb # output: atomistic.pdb, ref_frames.dat 

##### Map into a coarse-grained representation #####
python ${spn2cdir}/utils/pdb2cg_dna.py atomistic.pdb # output: dna_conf.in, dna_bonds.in, dna_angles.in, dna_dihedrals.in
                                                     #         in00_conf.xyz, in00_conf.crd, in00_cvmd.psf 
#python ${spn2cdir}/utils/pdb2cg_dna.py align_dna/atomistic_dna_align.pdb # Use when you need your DNA with customized coordinate

##### Generate LAMMPS input files #####
${spn2cdir}/DSIM_ICNF/icnf.exe dnaSeq.txt 1 1 . 0 # conf_lammps.in, in00_conf.xyz, in00_cvmd.psf 
${spn2cdir}/utils/replace_atoms.sh conf_lammps.in dna_conf.in bdna_curv_conf.in # bdna_curv_conf.in tmp.start, tmp.end
python ${spn2cdir}/utils/make_list_files.py bdna_curv_conf.in # in00_bond.list, in00_angl.list, in00_dihe.list 

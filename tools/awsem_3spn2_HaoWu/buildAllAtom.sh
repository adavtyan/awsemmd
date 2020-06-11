#!/bin/bash

python BuildAllAtomsFromLammps_multiChain_wDNA.py dump.lammpstrj dump.pdb -seq protein.seq  \
        -dnaPdb ./cg_dna.pdb -dnaBond ./dnaBondFromPsf.txt

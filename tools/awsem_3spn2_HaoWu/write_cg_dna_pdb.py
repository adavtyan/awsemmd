#!/usr/local/bin/python

#######################################################################

# Copyright (2017) Hao Wu
# Papoian Research Group
# University of Maryland, College Park
# http://papoian.chem.umd.edu/
# Last update: Jun 08 2017

#----------------------------------------------------------------------

# Write coarse grained DNA PDB file based on psf and coordinates
# Follow PDB format:
# http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM

#######################################################################

import sys

if len(sys.argv) != 4:
    print ("\n" + sys.argv[0] + " res_num crd_file pdb_file\n")
    exit()

res_num = int(sys.argv[1])
crdFile = sys.argv[2]
pdbFile = sys.argv[3]

# chainID: subject to change if more DNA chain is needed!
chain1 = 'Y'
chain2 = 'Z'

atomName2Ele = {'S': 'C', 'P': 'P', 'A': 'C', 'C': 'C', 'G': 'C', 'T': 'C'}

with open(pdbFile, 'w') as fw:
    fw.write('CRYST1  105.950  181.170  109.490  90.00  90.00  90.00 P 1           1\n')
    with open(crdFile, 'r') as fr:
        lines = fr.readlines()[2:] # skip first two lines and the last line
        for line in lines:
            line = line.split()
            atom_id = int(line[0])
            res_id = int(line[1])
            res_name = line[2]
            atom_name = line[3]
            x = float(line[4])
            y = float(line[5])
            z = float(line[6])
            fw.write('ATOM  ')
            fw.write('{:5d}'.format(atom_id))
            fw.write(' ')
            fw.write('{:^4}'.format(atom_name))
            fw.write(' ')
            fw.write(res_name)
            fw.write(' ')
            if res_id <= res_num:
                fw.write(chain1)
                fw.write('{:4d}'.format(res_id))
            else:
                fw.write(chain2)
                fw.write('{:4d}'.format(res_id - res_num))
            fw.write('    ')
            fw.write('{:8.3f}'.format(x))
            fw.write('{:8.3f}'.format(y))
            fw.write('{:8.3f}'.format(z))
            fw.write('  0.00')
            fw.write('  0.00')
            fw.write('          ')
            fw.write('{:>2}'.format(atomName2Ele[atom_name]))
            fw.write('\n')
    fw.write('END\n')

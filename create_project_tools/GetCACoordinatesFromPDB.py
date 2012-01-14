#!/shared/local/bin/python

import sys

def three2one(prot):
    """ translate a protein sequence from 3 to 1 letter code"""
    
    code = {"GLY" : "G", "ALA" : "A", "LEU" : "L", "ILE" : "I",
            "ARG" : "R", "LYS" : "K", "MET" : "M", "CYS" : "C",
            "TYR" : "Y", "THR" : "T", "PRO" : "P", "SER" : "S",
            "TRP" : "W", "ASP" : "D", "GLU" : "E", "ASN" : "N",
	    "GLN" : "Q", "PHE" : "F", "HIS" : "H", "VAL" : "V",
	    "M3L" : "K", "MSE" : "M", "CAS" : "C" }
    
    newprot = ""
    for a in prot:
        newprot += code.get(a, "?")
    
    return newprot

if len(sys.argv)!=3:
    print "\nGetCACoordinatesFromPDB.py PDB_Id Output_file\n"
    exit()

from Bio.PDB.PDBParser import PDBParser

p = PDBParser(PERMISSIVE=1)

struct_id = sys.argv[1]
filename = struct_id + ".pdb"

output_fn = ""
output_fn = sys.argv[2]

out = open( output_fn, 'w' )

xyz_CA = []

s = p.get_structure(struct_id, filename)
chains = s[0].get_list()
for ch in chains:
    sequance = []
    for res in ch:
	is_regular_res = res.has_id('CA') and res.has_id('O')
	res_id = res.get_id()[0]
        if (res_id==' ' or res_id=='H_MSE' or res_id=='H_M3L' or res_id=='H_CAS') and is_regular_res:
            sequance.append(res.get_resname())
            xyz_CA.append(res['CA'].get_coord())

out.write( str(len(xyz_CA)) )
out.write('\n')
for ixyz in xyz_CA:
	out.write( str(round(ixyz[0], 5)) )
	out.write('\t')
	out.write( str(round(ixyz[1], 5)) )
	out.write('\t')
	out.write( str(round(ixyz[2], 5)) )
	out.write('\n')

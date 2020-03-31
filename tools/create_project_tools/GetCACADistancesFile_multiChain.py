#!/usr/bin/python
#Modified from GetCACADistancesFile.py

import sys
from VectorAlgebra import *

def calc_dis(p1, p2):
    v = vector(p1, p2)
    return vabs(v)

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

def checkIfNative(ires, jres):
    xyz_CAi = ires['CA'].get_coord()
    xyz_CAj = jres['CA'].get_coord()
    v = vector(xyz_CAi, xyz_CAj)
    r = vabs(v)
    if r<12.0: return True
    else: return False

if len(sys.argv)==1:
    print "\nGetCACADistancesFile_multi.py PDB_Id Output_file \n"
    print "-s\tSplit into files for each chain"
#    sys.argv.append("1BG8")
    exit()

from Bio.PDB.PDBParser import PDBParser

p = PDBParser(PERMISSIVE=1)

struct_id = sys.argv[1]
filename = struct_id + ".pdb"

splite = False

output_fn = ""
if len(sys.argv)>2: output_fn = sys.argv[2]
if output_fn[-4:]==".dat": output_fn = output_fn[:-4]

if output_fn!="" and not splite:
    out = open( (output_fn+".dat"), 'w' )

k_bond = 100
k_angle = 20
k_dihedral = [1, 0.5]
epsilon = 1
epsilon2 = 1
sigma0 = 4

xyz_CA1 = []
xyz_CA2 = []
xyz_CA3 = []
xyz_CA4 = []

s = p.get_structure(struct_id, filename)
chains = s[0].get_list()
sequence = []
dis = []
all_res = []

for chain in chains:
    for resid in chain:
        is_regular_res = resid.has_id('CA') and resid.has_id('O')
        res_id = resid.get_id()[0]
        if (res_id==' ' or res_id=='H_MSE' or res_id=='H_M3L' or res_id=='H_CAS') and is_regular_res:
            all_res.append(resid)
            sequence.append(resid.get_resname())

for i in range( 0, len(all_res) ):
    dis.append([]);
    ires = all_res[i]
    xyz_CAi = ires['CA'].get_coord()
    for j in range( 0, len(all_res) ):
	jres = all_res[j]
	xyz_CAj = jres['CA'].get_coord()
	r = calc_dis(xyz_CAi, xyz_CAj)
	dis[i].append(r);

for ri in dis:
    for rij in ri:
        out.write( str(round(rij, 6)) )
        out.write( ' ' )
    out.write( '\n' )	
if not splite:
    out.write('\n')
    
out.close()

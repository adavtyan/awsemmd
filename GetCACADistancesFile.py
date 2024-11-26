#!/shared/local/bin/python

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
    print "\nReadingPDBFile.py PDB_Id [Output_file [-s]]\n"
    print "-s\tSplit into files for each chain"
#    sys.argv.append("1BG8")
    exit()

from Bio.PDB.PDBParser import PDBParser

p = PDBParser(PERMISSIVE=1)

struct_id = sys.argv[1]
filename = struct_id + ".pdb"

splite = False
for av in sys.argv:
    if av=="-s":
        splite = True
        sys.argv.remove(av)
        break

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
for ch in chains:
    sequance = []
    dis = []
    if output_fn!="":
	pass
    else:
        print "Chain:", ch.get_id()
    four_res = [None, None, None, None]
    all_res = []
    for res in ch:
	is_regular_res = res.has_id('CA') and res.has_id('O')
	res_id = res.get_id()[0]
        if (res_id==' ' or res_id=='H_MSE' or res_id=='H_M3L' or res_id=='H_CAS') and is_regular_res:
            all_res.append(res)
            sequance.append(res.get_resname())

    for i in range( 0, len(all_res) ):
	dis.append([]);
	ires = all_res[i]
	xyz_CAi = ires['CA'].get_coord()
	for j in range( 0, len(all_res) ):
		jres = all_res[j]
		xyz_CAj = jres['CA'].get_coord()
		r = calc_dis(xyz_CAi, xyz_CAj)
		dis[i].append(r);

    if output_fn!="":
        if splite and len(sequance)==0: continue
        if splite:
	    if len(chains)==1:
		file_name = output_fn+".dat"
	    else:
		file_name = output_fn+"_"+ch.get_id()+".dat"
            out = open( file_name, 'w' )
	
	for ri in dis:
		for rij in ri:
			out.write( str(round(rij, 6)) )
			out.write( ' ' )
		out.write( '\n' )	
        if not splite:
            out.write('\n')
        if splite:
            out.close()
    else:
        print three2one(sequance)
	for ri in dis:
		for rij in ri:
			print str(round(rij, 2)), ' ',
		print
        print '\n'

if output_fn!="" and not splite:
    out.close()

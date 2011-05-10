#!/shared/local/bin/python

# This is a script written by Aram Davtyan that was modified by Nick Schafer on 4/27/2011
# It differs from the GetCACADistances.py in that it returns all of the applicable
# CACA CACB CBCA and CBCB distances

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
	    "GLN" : "Q", "PHE" : "F", "HIS" : "H", "VAL" : "V"}
    
    newprot = ""
    for a in prot:
        newprot += code.get(a, "?")
    
    return newprot

if len(sys.argv)==1:
    print "\nReadingPDBFile.py PDB_Id [Output_file [-s]]\n"
    print "-s\tSplit into files for each chain"
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
if output_fn[-5:]==".data": output_fn = output_fn[:-5]

if output_fn!="" and not splite:
    out = open( (output_fn+".data"), 'w' )

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
        if res.get_id()[0]==' ' and is_regular_res:
            all_res.append(res)
            sequance.append(res.get_resname())

    for i in range( 0, len(all_res) ):
	dis.append([]);
	ires = all_res[i]
        if ires.get_resname()=='GLY':
            xyz_i_ca = ires['CA'].get_coord()
            xyz_i_ca[0] = round(xyz_i_ca[0]*100)/100.0
            xyz_i_ca[1] = round(xyz_i_ca[1]*100)/100.0
            xyz_i_ca[2] = round(xyz_i_ca[2]*100)/100.0
        else:
            xyz_i_ca = ires['CA'].get_coord()
            xyz_i_ca[0] = round(xyz_i_ca[0]*100)/100.0
            xyz_i_ca[1] = round(xyz_i_ca[1]*100)/100.0
            xyz_i_ca[2] = round(xyz_i_ca[2]*100)/100.0
            xyz_i_cb = ires['CB'].get_coord()
            xyz_i_cb[0] = round(xyz_i_cb[0]*100)/100.0
            xyz_i_cb[1] = round(xyz_i_cb[1]*100)/100.0
            xyz_i_cb[2] = round(xyz_i_cb[2]*100)/100.0
	for j in range( 0, len(all_res) ):
		jres = all_res[j]
                if jres.get_resname()=='GLY':
                    xyz_j_ca = jres['CA'].get_coord()
                    xyz_j_ca[0] = round(xyz_j_ca[0]*100)/100.0
                    xyz_j_ca[1] = round(xyz_j_ca[1]*100)/100.0
                    xyz_j_ca[2] = round(xyz_j_ca[2]*100)/100.0
                else:
                    xyz_j_ca = jres['CA'].get_coord()
                    xyz_j_ca[0] = round(xyz_j_ca[0]*100)/100.0
                    xyz_j_ca[1] = round(xyz_j_ca[1]*100)/100.0
                    xyz_j_ca[2] = round(xyz_j_ca[2]*100)/100.0
                    xyz_j_cb = jres['CB'].get_coord()
                    xyz_j_cb[0] = round(xyz_j_cb[0]*100)/100.0
                    xyz_j_cb[1] = round(xyz_j_cb[1]*100)/100.0
                    xyz_j_cb[2] = round(xyz_j_cb[2]*100)/100.0

                    
                rcaca = calc_dis(xyz_i_ca, xyz_j_ca)
		dis[i].append(rcaca);
                
                if jres.get_resname()!='GLY':
                    rcacb = calc_dis(xyz_i_ca, xyz_j_cb)
                    dis[i].append(rcacb);
                else:
                    dis[i].append(0.0);
                
                if ires.get_resname()!='GLY':
                    rcbca = calc_dis(xyz_i_cb, xyz_j_ca)
                    dis[i].append(rcbca);
                else:
                    dis[i].append(0.0);
                
                if ires.get_resname()!='GLY':
                    rcbcb = calc_dis(xyz_i_cb, xyz_j_cb)
                    dis[i].append(rcbcb);
                else:
                    dis[i].append(0.0);

    if output_fn!="":
        if splite and len(sequance)==0: continue
        if splite:
	    if len(chains)==1:
		file_name = output_fn+".data"
	    else:
		file_name = output_fn+"_"+ch.get_id()+".data"
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

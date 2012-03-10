#!/usr/bin/python

# ----------------------------------------------------------------------
# Copyright (2010) Aram Davtyan and Garegin Papoian

# Papoian's Group, University of Maryland at Collage Park
# http://papoian.chem.umd.edu/

# Last Update: 03/04/2011
# ----------------------------------------------------------------------

import sys
from VectorAlgebra import *

def calc_dihedral_angle(p1, p2, p3, p4):
    v1 = vector(p1, p2)
    v2 = vector(p2, p3)
    v3 = vector(p3, p4)
    return 180*dihedral_angle(v1, v2, v3)/pi

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
if output_fn[-3:]==".se": output_fn = output_fn[:-3]

if output_fn!="" and not splite:
    out = open( (output_fn+".se"), 'w' )
    se_out = open( (output_fn+".seq"), 'w' )

standartPhi = 60.0
standartPsi = 180.0

s = p.get_structure(struct_id, filename)
chains = s[0].get_list()
for ch in chains:
    sequance = []
    angles = []
    phi = standartPhi
    psi = standartPsi
    if output_fn!="":
	pass
    else:
        print "Chain:", ch.get_id()
    two_res = [None, None]
    for res in ch:
        is_regular_res = res.has_id('N') and res.has_id('CA') and res.has_id('C')
        #if res.get_id()[0]==' ' and is_regular_res:
        res_id = res.get_id()[0]
        if (res_id==' ' or res_id=='H_MSE' or res_id=='H_M3L' or res_id=='H_CAS') and is_regular_res:
            two_res.append(res)
            p_res = two_res.pop(0)
            if p_res:
                sequance.append(p_res.get_resname())
		print p_res.get_resname()
            if two_res[0] and two_res[1]:
                if two_res[0].get_id()[1]+1!=two_res[1].get_id()[1]:
                    print "Error: Wrong residue order"
                xyz_N1 = two_res[0]['N'].get_coord()
                xyz_CA1 = two_res[0]['CA'].get_coord()
                xyz_C1 = two_res[0]['C'].get_coord()
                xyz_N2 = two_res[1]['N'].get_coord()
                xyz_CA2 = two_res[1]['CA'].get_coord()
                xyz_C2 = two_res[1]['C'].get_coord()
                psi = calc_dihedral_angle(xyz_N1, xyz_CA1, xyz_C1, xyz_N2)
                angles.append([phi, psi])
                phi = calc_dihedral_angle(xyz_C1, xyz_N2, xyz_CA2, xyz_C2)

    psi = standartPsi
    if two_res[1] and two_res[1].has_id('O'):
        xyz_N = two_res[1]['N'].get_coord()
        xyz_CA = two_res[1]['CA'].get_coord()
        xyz_C = two_res[1]['C'].get_coord()
        xyz_O = two_res[1]['O'].get_coord()
        psi = calc_dihedral_angle(xyz_N, xyz_CA, xyz_C, xyz_O)
        if psi>0:  psi-= 180.0
        else: psi+= 180.0
        angles.append([phi, psi])
    if two_res[0]: sequance.append(two_res[0].get_resname())
    if two_res[1]: sequance.append(two_res[1].get_resname())
    if output_fn!="":
        if splite and len(sequance)==0: continue
        if splite:
	    if len(chains)==1:
		file_name = output_fn+".se"
	    else:
		file_name = output_fn+"_"+ch.get_id()+".se"
            out = open( file_name, 'w' )
            se_out = open( (file_name+"q"), 'w' )
        out.write(three2one(sequance))
	se_out.write(three2one(sequance))
        out.write('\n')
        for ia in angles:
            out.write( str(round(ia[0],3)) )
            out.write(' ')
            out.write( str(round(ia[1],3)) )
            out.write('\n')
        if not splite:
            out.write('\n')
        if splite:
            out.close()
            se_out.close()
    else:
        print three2one(sequance)
        for ia in angles:
            print str(round(ia[0],3)), str(round(ia[1],3))
        print '\n'

if output_fn!="" and not splite:
    out.close()
    se_out.close()

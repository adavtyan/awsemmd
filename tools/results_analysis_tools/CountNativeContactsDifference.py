#!/usr/bin/python

# ----------------------------------------------------------------------
# Copyright (2010) Aram Davtyan and Garegin Papoian

# Papoian's Group, University of Maryland at Collage Park
# http://papoian.chem.umd.edu/

# Last Update: 03/04/2011
# ----------------------------------------------------------------------

import sys
from VectorAlgebra import *
from Bio.PDB.PDBParser import PDBParser

def calc_dihedral_angle(p1, p2, p3, p4):
    v1 = vector(p1, p2)
    v2 = vector(p2, p3)
    v3 = vector(p3, p4)
    return pi*dihedral_angle(v1, v2, v3)/180

def calc_angle(p1, p2, p3):
    v1 = vector(p1, p2)
    v2 = vector(p2, p3)
    return vangle(v1, v2)

def calc_bond(p1, p2):
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

def checkIfNative(xyz_CAi, xyz_CAj):
    v = vector(xyz_CAi, xyz_CAj)
    r = vabs(v)
    if r<12.0: return True
    else: return False

def isNative(r):
	if r<12.0: return True
    	else: return False

def CalcGoModelCoeffs(atoms):
	for i in range( 0, len(atoms) ):
		sigma.append([])
		for j in range( i+4, len(atoms) ):
			if abs(j-i)<4:
				sigma[i].append(0.0)
			elif j>i:
				xyz_CAi = atoms[i]
				xyz_CAj = atoms[j]
				r = sqrt((xyz_CAj[0]- xyz_CAi[0])**2 + (xyz_CAj[1]- xyz_CAi[1])**2 + (xyz_CAj[2]- xyz_CAi[2])**2)
				sigma[i].append(r)
			else:
				sigma[i].append(sigma[j][i])

def Compare():
	dContactTotal = 0.0
	for i in range(0, len(sigma) ):
		for j in range(0, len(sigma[i]) ):
			if isNative(sigmaN[i][j]):
				if sigma[i][j]>1.2*sigmaN[i][j]:
					dContactTotal = dContactTotal + 1
	if frac:
		dContactTotal = (nNative - dContactTotal)/nNative
	diff.append(dContactTotal)

#Variables

ca_atoms_pdb = []
sigmaN = []
sigma = []
sigma0 = 4
diff = []
frac = False # Output fraction of native contacts
nNative = 0

if len(sys.argv)<=2:
    print ("\nExtractGoModelCGCoeffs.py Input_file PDB_id [Output_file [-s]] [-f]\n")
    print ("-s\tSplit into files for each chain\n")
    print ("-f\tOutput fraction of native contacts\n")
    exit()

filename = sys.argv[1]
pdb_id = sys.argv[2]

if pdb_id[-4:].lower()==".pdb":
	pdb_file = pdb_id
else:
	pdb_file = pdb_id + ".pdb"

splite = False
for av in sys.argv:
    if av=="-s":
        splite = True
#        sys.argv.reimove(av)
    if av=="-f":
	frac = True

output_fn = ""
if len(sys.argv)>3: output_fn = sys.argv[3]
if output_fn[-5:]==".data": output_fn = output_fn[:-5]


p = PDBParser(PERMISSIVE=1)

s = p.get_structure(pdb_id, pdb_file)

chains = s[0].get_list()
chain = chains[0]
for res in chain:
	is_regular_res = res.has_id('CA') and res.has_id('O')
	res_id = res.get_id()[0]
        if (res_id==' ' or res_id=='H_MSE' or res_id=='H_M3L') and is_regular_res:
		ca_atoms_pdb.append(res['CA'].get_coord())

for i in range( 0, len(ca_atoms_pdb) ):
	sigmaN.append([])
	for j in range( i+4, len(ca_atoms_pdb) ):
		if abs(j-i)<4: 
			sigmaN[i].append(0.0)
		elif j>i:
			xyz_CAi = ca_atoms_pdb[i]
			xyz_CAj = ca_atoms_pdb[j]
			v = vector(xyz_CAi, xyz_CAj)
			sigmaN[i].append(vabs(v))
			if isNative(vabs(v)): nNative = nNative + 1
		else:
			sigmaN[i].append(sigmaN[j][i])

if output_fn!="" and not splite:
    out = open( (output_fn+".data"), 'w' )


infile = open(filename, 'r')

box = []
A = []
nFrame = 0
nAtoms = 0
ca_atoms = []
for l in infile:
	l = l.strip()
	if l[:5]=="ITEM:":
		item = l[6:]
	else:
		if item == "TIMESTEP":
			if len(ca_atoms)>0:
				sigma = []
				CalcGoModelCoeffs(ca_atoms)
				Compare()
			nFrame = nFrame + 1
			ca_atoms = []
			step = int(l)
			box = []
			A = []
		elif item == "NUMBER OF ATOMS":
			nAtoms = int(l)
		elif item[:10] == "BOX BOUNDS":
			box.append(l)
			l = l.split()
			A.append([float(l[0]), float(l[1])])
		elif item[:5] == "ATOMS":
			l = l.split() 
			if l[1]=='1':
				x = float(l[2])
				y = float(l[3])
				z = float(l[4])
				x = (A[0][1] - A[0][0])*x + A[0][0]
				y = (A[1][1] - A[1][0])*y + A[1][0]
				z = (A[2][1] - A[2][0])*z + A[2][0]
				ca_atoms.append( [x, y, z] )
infile.close()

if len(ca_atoms)>0:
	sigma = []
	CalcGoModelCoeffs(ca_atoms)
	Compare()

ldiff = len(diff)

print ("Number of native contacts: ", nNative)

if output_fn!="":
#        if splite and len(sequance)==0: continue
        if splite:
	    if len(chains)==1:
		file_name = output_fn+".data"
	    else:
		file_name = output_fn+"_"+ch.get_id()+".data"
            out = open( file_name, 'w' )
		
	idf = 0
	for df in diff:
		idf = idf + 1
		if frac: 
			out.write(str(round(df,4)))
		else:
			out.write(str(int(df)))
		if idf != ldiff:
			out.write(" ")
	out.write("\n")
	
        if not splite:
            out.write('\n')
        if splite:
            out.close()
else:
	idf = 0
	for df in diff:
		idf = idf + 1
		if frac:
			print (str(round(df,4)), end="")
		else:
			print (str(int(df)), end="")
		if idf!=ldiff:
			print (" ", end="")
	print ("\n")

if output_fn!="" and not splite:
    out.close()

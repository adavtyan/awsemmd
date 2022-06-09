#!/usr/bin/python

# ----------------------------------------------------------------------
# Copyright (2010) Aram Davtyan and Garegin Papoian

# Papoian's Group, University of Maryland at Collage Park
# http://papoian.chem.umd.edu/

# Last Update: 03/04/2011
# ----------------------------------------------------------------------

from pylab import *

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
			if abs(j-i)<3:
				sigma[i].append(0.0)
			elif j>i:
				xyz_CAi = atoms[i]
				xyz_CAj = atoms[j]
				r = sqrt((xyz_CAj[0]- xyz_CAi[0])**2 + (xyz_CAj[1]- xyz_CAi[1])**2 + (xyz_CAj[2]- xyz_CAi[2])**2)
				sigma[i].append(r)
			else:
				sigma[i].append(sigma[j][i])

#Variables

ca_atoms_pdb = []
sigmaN = []
sigma = []
sigma0 = 4

if len(sys.argv)<=3:
    print "\nExtractGoModelCGCoeffs.py Input_file PDB_id snapshot\n"
    print "-s\tSplit into files for each chain"
    exit()

filename = sys.argv[1]
pdb_id = sys.argv[2]

pdb_file = pdb_id + ".pdb"

frame = int(sys.argv[3])

p = PDBParser(PERMISSIVE=1)

s = p.get_structure(pdb_id, pdb_file)

chains = s[0].get_list()
chain = chains[0]
for res in chain:
	is_regular_res = res.has_id('CA') and res.has_id('O')
        if res.get_id()[0]==' ' and is_regular_res:
		ca_atoms_pdb.append(res['CA'].get_coord())

for i in range( 0, len(ca_atoms_pdb) ):
	sigmaN.append([])
	for j in range( i+4, len(ca_atoms_pdb) ):
		if abs(j-i)<3: 
			sigmaN[i].append(0.0)
		elif j>i:
			xyz_CAi = ca_atoms_pdb[i]
			xyz_CAj = ca_atoms_pdb[j]
			v = vector(xyz_CAi, xyz_CAj)
			sigmaN[i].append(vabs(v))
		else:
			sigmaN[i].append(sigmaN[j][i])

figure()

xlabel('Native')
ylabel('Prediction')
title('Protein contact map')
grid(True)

ln = len(ca_atoms_pdb)

plot([0, ln],[0, ln], color='black')
axis([0, ln+1, 0, ln+1])

for i in range( 0, len(ca_atoms_pdb) ):
        for j in range( i+4, len(ca_atoms_pdb) ):
		if isNative(sigmaN[i][j-i-4]):
#			plot([i+1],[j],'bs')
			plot([j],[i+1],'rs')

infile = open(filename, 'r')

box = []
A = []
nAtoms = 0
nFrame = 0
ca_atoms = []
found = False
for l in infile:
	l = l.strip()
	if l[:5]=="ITEM:":
		item = l[6:]
	else:
		if item == "TIMESTEP":
			if found: break
			step = int(l)
			if nFrame==frame: found = True
			nFrame = nFrame + 1
		if found:
			if item == "NUMBER OF ATOMS":
				nAtoms = int(l)
			elif item == "BOX BOUNDS":
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
	for i in range( 0, len(ca_atoms) ):
	        for j in range( i+4, len(ca_atoms) ):
        	        if isNative(sigma[i][j-i-4]):
                	        plot([i+1],[j],'bs')
#                        	plot([j],[i+1],'rs')

show()

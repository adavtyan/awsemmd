#!/usr/bin/python

# ----------------------------------------------------------------------
# Copyright (2010) Aram Davtyan and Garegin Papoian

# Papoian's Group, University of Maryland at Collage Park
# http://papoian.chem.umd.edu/

# Modified by Nick Schafer on 6/3/11
# ----------------------------------------------------------------------

import sys
from VectorAlgebra import *
from Bio.PDB.PDBParser import PDBParser

numFoldons = 0

# Read the foldon definition file
def readFoldonFile(fileName):
    foldonFile = open(fileName, 'r')
    iLine = 0
    for line in foldonFile:
        line = line.split()
        for residue in line:
            residue = residue.strip()
            foldonmap[int(residue)-1] = iLine

        iLine = iLine + 1

    numFoldons = iLine
    return numFoldons

# Given a distance r, determine if it is within the contact threshold
def isContact(r):
	if r<8.0: return True
    	else: return False

# Given a set of coordinates (atoms), determine all pairwise distances
def CalcCADistances(atoms):
	for i in range( 0, len(atoms) ):
		dist.append([])
		for j in range( 0, len(atoms) ):
                    xyz_CAi = atoms[i]
                    xyz_CAj = atoms[j]
                    r = sqrt((xyz_CAj[0]- xyz_CAi[0])**2 + (xyz_CAj[1]- xyz_CAi[1])**2 + (xyz_CAj[2]- xyz_CAi[2])**2)
                    dist[i].append(r)

# Given a set of pairwise distances and the native pairwise distances
def CountNativeContacts():
    # Zero contacts array
    for i in xrange(numFoldons):
        for j in xrange(numFoldons):
            numNativeContacts[i][j] = 0

    # Count Native Contacts
    for i in range(0, len(dist) ):
        for j in range(0, len(dist) ):
            if isContact(distN[i][j]) and isContact(dist[i][j]) and abs(i-j)>3:
                iFoldon = foldonmap[i]
                jFoldon = foldonmap[j]
                numNativeContacts[iFoldon][jFoldon] = numNativeContacts[iFoldon][jFoldon] + 1
                if iFoldon!=jFoldon:
                    numNativeContacts[jFoldon][iFoldon] = numNativeContacts[jFoldon][iFoldon] + 1

#Variables
ca_atoms_pdb = []
distN = []
dist = []
nNative = 0

if len(sys.argv)<=3:
    print "\nCountFoldonContacts.py Input_file PDB_id Foldon_file Output_file_name\n"
    exit()

filename = sys.argv[1]
pdb_id = sys.argv[2]
foldon_file = sys.argv[3]

if pdb_id[-4:].lower()==".pdb":
	pdb_file = pdb_id
else:
	pdb_file = pdb_id + ".pdb"

output_fn = ""
if len(sys.argv)>4: output_fn = sys.argv[4]
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

# Calculate Native Distances
for i in range( 0, len(ca_atoms_pdb) ):
	distN.append([])
	for j in range( 0, len(ca_atoms_pdb) ):
            xyz_CAi = ca_atoms_pdb[i]
            xyz_CAj = ca_atoms_pdb[j]
            v = vector(xyz_CAi, xyz_CAj)
            distN[i].append(vabs(v))
            if isContact(vabs(v)) and abs(i-j)>3:
                nNative = nNative + 1
                        
# Initialize foldon map
foldonmap = []
for i in xrange(len(ca_atoms_pdb)):
    foldonmap.append(0)
# Read in foldon file
numFoldons = readFoldonFile(foldon_file)
# Declare numNativeContacts
numNativeContacts = []
# Initialize contacts array
for i in xrange(numFoldons):
    numNativeContacts.append([])
    for j in xrange(numFoldons):
        numNativeContacts[i].append(0)

if output_fn!="":
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
				dist = []
				CalcCADistances(ca_atoms)
				CountNativeContacts()
                                print numNativeContacts
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
    dist = []
    CalcCADistances(ca_atoms)
    CountNativeContacts()
    print numNativeContacts



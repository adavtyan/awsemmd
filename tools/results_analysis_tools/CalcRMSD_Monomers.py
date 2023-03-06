#!/usr/bin/python

# Created by Aram Davtyan, Aug 2018
# Use BioPython to calculate monomer RMSD from the ref. structure.
# Only C-alpha atoms are used.

# ----------------------------------------------------------------------
# Copyright (2010) Aram Davtyan and Garegin Papoian

# Papoian's Group, University of Maryland at Collage Park
# http://papoian.chem.umd.edu/

# Last Update: 08/13/2018
# ----------------------------------------------------------------------

import sys

##
import numpy
##

from VectorAlgebra import *

atom_type = {'1' : 'C', '2' : 'N', '3' : 'O', '4' : 'C', '5' : 'H', '6' : 'C'}
atom_desc = {'1' : 'C-Alpha', '2' : 'N', '3' : 'O', '4' : 'C-Beta', '5' : 'H-Beta', '6' : 'C-Prime'}
PDB_type = {'1' : 'CA', '2' : 'N', '3' : 'O', '4' : 'CB', '5' : 'HB', '6' : 'C' }

class PDB_Atom:
	no = 0
	ty = ''
	res = 'UNK'
	res_no = 0
	x = 0.0
	y = 0.0
	z = 0.0
	atm = 'C'
	
	def __init__(self, no, ty, res, res_no, x, y, z, atm):
		self.no = no
		self.ty = ty
		self.res = res
		self.res_no = res_no
		self.x = x
		self.y = y
		self.z = z
		self.atm = atm
		
	def write_(self, f):
		f.write('ATOM')
		f.write(('       '+str(self.no))[-7:])
		f.write('  ')
		f.write((self.ty+'    ')[:4])
		f.write(self.res)
		f.write(' ')
		f.write('T')
		f.write(('    '+str(self.res_no))[-4:])
		f.write(('            '+str(round(self.x,3)))[-12:])
		f.write(('        '+str(round(self.y,3)))[-8:])
		f.write(('        '+str(round(self.z,3)))[-8:])
		f.write('  1.00')
		f.write('  0.00')
		f.write(('            '+self.atm)[-12:]+'  ')
		f.write('\n')

class Atom:
	No = 0
	ty = ''
	x = 0.0
	y = 0.0
	z = 0.0
	desc = ''
	
	def __init__(self, No, ty, No_m, x, y, z, desc=''):
		self.No = No
		self.ty = ty
		self.No_m = No_m
		self.x = x
		self.y = y
		self.z = z
		self.desc = desc
	
	def write_(self, f):
		f.write(str(self.No))
		f.write(' ')
		f.write(PDB_type[self.No_m])
		f.write(' ')
		f.write(str(round(self.x,8)))
		f.write(' ')
		f.write(str(round(self.y,8)))
		f.write(' ')
		f.write(str(round(self.z,8)))
		f.write(' ')
		f.write(self.desc)
		f.write('\n')

if len(sys.argv)!=4:
	print ("\CalcRMSD.py PDB_Id Input_file(lammpstrj) Output_file(rmsd)\n")
	exit()

struct_id = sys.argv[1]
if struct_id[-4:].lower()==".pdb":
        pdb_file = struct_id
        struct_id = struct_id[:-4]
else:
        pdb_file = struct_id + ".pdb"

lammps_file = sys.argv[2]
output_file = sys.argv[3]

n_atoms = 0
i_atom = 0
item = ''
step = 0
chain_id_list = []
chain_id = []
ca_atoms_pdb = []
ca_atoms = []
box = []
A = []

out = open(output_file, 'w')

from Bio.PDB.PDBParser import PDBParser

##
from Bio.SVDSuperimposer import SVDSuperimposer
##

p = PDBParser(PERMISSIVE=1)

##
def computeRMSD():
	if len(ca_atoms)!=len(ca_atoms_pdb):
		print ("Error. Length mismatch!")
		exit()
	l = len(ca_atoms)

	res = {}
	for ch in chain_id_list:
		fixed_coord = []
		moving_coord = []
		
		for i in range(l):
			if chain_id[i]==ch:
				fixed_coord.append([ca_atoms_pdb[i][0], ca_atoms_pdb[i][1], ca_atoms_pdb[i][2]])
				moving_coord.append([ca_atoms[i][0], ca_atoms[i][1], ca_atoms[i][2]])

		if len(fixed_coord)>0:
			fixed_coord = numpy.array(fixed_coord)
			moving_coord = numpy.array(moving_coord)

			sup = SVDSuperimposer()
			sup.set(fixed_coord, moving_coord)
			sup.run()
			rms = sup.get_rms()

			res[ch] = rms

	return res

s = p.get_structure(struct_id, pdb_file)
chains = s[0].get_list()
for chain in chains:
        chain_id_list.append(chain.get_id())
	for res in chain:
		is_regular_res = res.has_id('CA') and res.has_id('O')
		res_id = res.get_id()[0]
        	if (res_id==' ' or res_id=='H_MSE' or res_id=='H_M3L') and is_regular_res:
			ca_atoms_pdb.append(res['CA'].get_coord())
			chain_id.append(chain.get_id())

first = True
lfile = open(lammps_file)
for l in lfile:
	l = l.strip()
	if l[:5]=="ITEM:":
		item = l[6:]
	else:
		if item == "TIMESTEP":
			if len(ca_atoms)>0:
				rmsd = computeRMSD()
				if first:
					out.write("#")
					for ch in chain_id_list:
						if rmsd.has_key(ch): out.write(" %s" % ch)
					out.write("\n")
					first = False
				for ch in chain_id_list:
					if rmsd.has_key(ch):
						out.write(str(round(rmsd[ch],3)))
						out.write(' ')
				out.write("\n")
				n_atoms = len(ca_atoms)
			step = int(l)
			ca_atoms = []
			box = []
			A = []
		elif item == "NUMBER OF ATOMS":
			n_atoms = int(l)
		elif item[:10] == "BOX BOUNDS":
			box.append(l)
			l = l.split()
			A.append([float(l[0]), float(l[1])])
		elif item[:5] == "ATOMS":
			l = l.split()
			i_atom = l[0]
			x = float(l[2])
			y = float(l[3])
			z = float(l[4])
			x = (A[0][1] - A[0][0])*x + A[0][0]
			y = (A[1][1] - A[1][0])*y + A[1][0]
			z = (A[2][1] - A[2][0])*z + A[2][0]
			desc = atom_desc[l[1]]
			if desc=='C-Alpha':
#				atom = Atom(i_atom, atom_type[l[1]], l[1], x, y, z, desc)
				atom = [x,y,z]
				ca_atoms.append(atom)
lfile.close()

if len(ca_atoms)>0:
	rmsd = computeRMSD()
	for ch in chain_id_list:
		if rmsd.has_key(ch):
			out.write(str(round(rmsd[ch],3)))
			out.write(' ')
	out.write("\n")
	n_atoms = len(ca_atoms)

out.close()

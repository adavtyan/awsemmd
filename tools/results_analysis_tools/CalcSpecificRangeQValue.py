#!/usr/bin/python

# ----------------------------------------------------------------------
# Copyright (2010) Aram Davtyan and Garegin Papoian

# Papoian's Group, University of Maryland at Collage Park
# http://papoian.chem.umd.edu/

# Last Update: 03/04/2011
# ----------------------------------------------------------------------

import sys
from VectorAlgebra import *

#from Bio.PDB.PDBParser import PDBParser

atom_type = {'1' : 'C', '2' : 'N', '3' : 'O', '4' : 'C', '5' : 'H', '6' : 'C'}
atom_desc = {'1' : 'C-Alpha', '2' : 'N', '3' : 'O', '4' : 'C-Beta', '5' : 'H-Beta', '6' : 'C-Prime'}
PDB_type = {'1' : 'CA', '2' : 'N', '3' : 'O', '4' : 'CB', '5' : 'HB', '6' : 'C' }

class PDB_Atom:
	no = 0
	ty = ''
	mol = 0
	res = 'UNK'
	res_no = 0
	x = 0.0
	y = 0.0
	z = 0.0
	atm = 'C'
	
	def __init__(self, no, ty, mol, res, res_no, x, y, z, atm):
		self.no = no
		self.ty = ty
		self.mol = mol
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
		f.write(self.mol)
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

if len(sys.argv)!=7:
	print "\nCalcQValue.py PDB_Id Input_file Output_file sigma_exp minsep maxsep [-i]\n"
	print
	print "\t\t-i\tcalculate individual q values for each chain"
	print
	exit()

splitq = False
for iarg in range(0, len(sys.argv)):
	if sys.argv[iarg]=="-i":
		splitq = True
		sys.argv.pop(iarg)

struct_id = sys.argv[1]
if struct_id[-4:].lower()==".pdb":
	pdb_file = struct_id
else:
	pdb_file = struct_id + ".pdb"
lammps_file = sys.argv[2]

output_file = ""
if len(sys.argv)>3: output_file = sys.argv[3]

sigma_exp = float(sys.argv[4])

minsep = int(sys.argv[5])
maxsep = int(sys.argv[6])

n_atoms = 0
i_atom = 0
item = ''
step = 0
ca_atoms_pdb = []
pdb_chain_id = []
ca_atoms = []
box = []
A = []
sigma = []
sigma_sq = []


out = open(output_file, 'w')

from Bio.PDB.PDBParser import PDBParser

p = PDBParser(PERMISSIVE=1)

def computeQ():
	if len(ca_atoms)!=len(ca_atoms_pdb):
		print "Error. Length mismatch!"
		print "Pdb: ", len(ca_atoms_pdb), "trj: ", len(ca_atoms)
		exit()
	Q = {}
	norm = {}
	N = len(ca_atoms)
	for ia in range(0, N):
		for ja in range(ia+minsep, ia+maxsep):
			if (ja >= N):
				continue
			if (splitq and pdb_chain_id[ia]==pdb_chain_id[ja]) or not splitq:
				r = vabs(vector(ca_atoms[ia], ca_atoms[ja]))
				rn = vabs(vector(ca_atoms_pdb[ia], ca_atoms_pdb[ja]))
				dr = r - rn
				if splitq: index = pdb_chain_id[ia]
				else: index = 1
				if not Q.has_key(index):
					Q[index] = 0.0
					norm[index] = 0
				Q[index] = Q[index] + exp(-dr*dr/(2*sigma_sq[ja-ia]))
				norm[index] = norm[index] + 1
	for key in Q:
		Q[key] = Q[key]/norm[key]
	return Q

s = p.get_structure(struct_id, pdb_file)
chains = s[0].get_list()
#chain = chains[0]
ichain = 0
for chain in chains:
	ichain = ichain + 1
	for res in chain:
		is_regular_res = res.has_id('CA') and res.has_id('O')
		res_id = res.get_id()[0]
	        if (res_id==' ' or res_id=='H_MSE' or res_id=='H_M3L' or res_id=='H_CAS' ) and is_regular_res:
			ca_atoms_pdb.append(res['CA'].get_coord())
			pdb_chain_id.append(ichain)

for i in range(0, len(ca_atoms_pdb)+1):
	sigma.append( (1+i)**sigma_exp )
	sigma_sq.append(sigma[-1]*sigma[-1])

lfile = open(lammps_file)
for l in lfile:
	l = l.strip()
	if l[:5]=="ITEM:":
		item = l[6:]
	else:
		if item == "TIMESTEP":
			if len(ca_atoms)>0:
				q = computeQ()
				for key in q:
					out.write(str(round(q[key],3)))
					out.write(' ')
				out.write('\n')
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
	q = computeQ()
	for key in q:
		out.write(str(round(q[key],3)))
		out.write(' ')
	out.write('\n')
	n_atoms = len(ca_atoms)

out.close()

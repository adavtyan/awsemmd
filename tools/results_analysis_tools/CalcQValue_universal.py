#!/usr/bin/python

# ----------------------------------------------------------------------
# Copyright (2010) Aram Davtyan and Garegin Papoian

# Papoian's Group, University of Maryland at Collage Park
# http://papoian.chem.umd.edu/

# Last Update: 03/04/2011
# ----------------------------------------------------------------------

import os
from VectorAlgebra import *
import argparse
import gzip

from Bio.PDB.PDBParser import PDBParser

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

parser = argparse.ArgumentParser(description='Parameters for Q calculations.')
parser.add_argument('--reference', '-R', action='store', required=True,
										help='reference PDB file name')
parser.add_argument('--input', '-I', action='store', required=True,
										help='input lammps trajectory file name')
parser.add_argument('--output', '-O', action='store', required=True,
										help='output file name')
parser.add_argument('--cutoff', '-C', action='store', type=float,
										help='set cutoff for monomer and complex Q calculations. By default no cutoff is used')
parser.add_argument('--sigma_exp', '-E', action='store', type=float, default=0.15,
		help='set sigma exponent for monomer and complex calculations (default: 0.15)')
parser.add_argument('--min_sep', '-M', action='store', type=int, default=3,
		help='minimal sequnce seperation for monomer and complex Q calculations (default: 3)')
parser.add_argument('--interface_sigma', '-S', action='store', type=float, default=3.0,
		help='set interchain sigma value (default: 3.0)')
parser.add_argument('--interface_cutoff', '-IC', action='store', type=float, default=10.0,
		help='set cutoff for interface Q calculations, with values of zero or less interpreted as no cutoff required (default: 10.0)')
parser.add_argument('--chains', '-CH', action='store', type=str, nargs="+", default="ALL",
		help='specify chains to be used in the order needed (default:)')

args = parser.parse_args()

pdb_file = args.reference
struct_id = os.path.splitext(os.path.basename(pdb_file))[0]

lammps_file = args.input
output_file = args.output

sigma_exp = args.sigma_exp
sigma_interchain = args.interface_sigma
sigma_sq_interchain = sigma_interchain*sigma_interchain

min_sep = args.min_sep

b_cutoff = False
cutoff = 0.0
if args.cutoff is not None and args.cutoff>0.0:
  b_cutoff = True
  cutoff = args.cutoff
  
b_interface_cutoff = False
interface_cutoff = 0.0
if args.interface_cutoff is not None and args.interface_cutoff>0.0:
  b_interface_cutoff = True
  interface_cutoff = args.interface_cutoff

b_max_cutoff = b_cutoff and b_interface_cutoff
max_cutoff = max(cutoff, interface_cutoff)

use_chains = args.chains

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

p = PDBParser(PERMISSIVE=1)

def computeQ():
	if len(ca_atoms)!=len(ca_atoms_pdb):
		print ("Error. Length mismatch!")
		print ("Pdb: ", len(ca_atoms_pdb), "trj: ", len(ca_atoms))
		exit()
	Q = {'Complex' : 0.0}
	norm = {'Complex' : 0.0}
	N = len(ca_atoms)
	for ia in range(0, N):
		for ja in range(ia+1, N):
			if pdb_chain_id[ia]==pdb_chain_id[ja] and ja-ia<min_sep:
				continue

			rn = vabs(vector(ca_atoms_pdb[ia], ca_atoms_pdb[ja]))
			if b_max_cutoff and rn>max_cutoff:
				continue

			r = vabs(vector(ca_atoms[ia], ca_atoms[ja]))
			dr = r - rn

			if pdb_chain_id[ia]==pdb_chain_id[ja]:
				sigma_sqc = sigma_sq[ja-ia]
				index = f"Chain_{pdb_chain_id[ia]}"
			else:
				sigma_sqc = sigma_sq_interchain
				index = f"Chains_{pdb_chain_id[ia]}:{pdb_chain_id[ja]}"
			Q_one = exp(-0.5*dr*dr/sigma_sqc)

			if index not in Q:
				Q[index] = 0.0
				norm[index] = 0.0

			if not b_cutoff or rn<=cutoff:
				Q['Complex'] += Q_one
				norm['Complex'] += 1.0
				if pdb_chain_id[ia]==pdb_chain_id[ja]:
					Q[index] += Q_one
					norm[index] += 1.0
			if (not b_interface_cutoff or rn<=interface_cutoff) and pdb_chain_id[ia]!=pdb_chain_id[ja]:
				Q[index] += Q_one
				norm[index] += 1.0

	for key in Q:
		Q[key] = Q[key]/norm[key]
	return Q

s = p.get_structure(struct_id, pdb_file)
model = s[0]
if use_chains=='' or use_chains=='ALL' or use_chains=='all':
	chains = model.get_list()
else:
	chains = []
	for chain in use_chains:
		if chain in model: 
			chains.append(model[chain])
		else:
			print ("Error. the specified chain is not found")
			exit ()

for chain in chains:
	ichain = chain.get_id()
	for res in chain:
		is_regular_res = res.has_id('CA') and res.has_id('O')
		res_id = res.get_id()[0]
		if (res_id==' ' or res_id=='H_MSE' or res_id=='H_M3L' or res_id=='H_CAS' ) and is_regular_res:
			ca_atoms_pdb.append(res['CA'].get_coord())
			pdb_chain_id.append(ichain)

for i in range(0, len(ca_atoms_pdb)+1):
	sigma.append( (1+i)**sigma_exp )
	sigma_sq.append(sigma[-1]*sigma[-1])

b_first = True
binary_file = lammps_file.endswith('.gz')
if binary_file:
	lfile = gzip.open(lammps_file, 'rb')
else:
	lfile = open(lammps_file)
for l in lfile:
	if binary_file: l = l.decode()
	l = l.strip()
	if l[:5]=="ITEM:":
		item = l[6:]
	else:
		if item == "TIMESTEP":
			if len(ca_atoms)>0:
				q = computeQ()
				if b_first:
					out.write('#')
					for key in q:
						if len(q)!=2 or key!='Complex':
							out.write(" "+key)
					out.write('\n')
					b_first = False
				for key in q:
					if len(q)!=2 or key!='Complex':
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
		if len(q)!=2 or key!='Complex':
			out.write(str(round(q[key],3)))
			out.write(' ')
	out.write('\n')
	n_atoms = len(ca_atoms)

out.close()

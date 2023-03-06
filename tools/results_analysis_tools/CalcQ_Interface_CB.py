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

#if len(sys.argv)!=5 and len(sys.argv)!=4:
if len(sys.argv)!=4 :
	print ("\n.py native_PDB_Id dump_file Output_Qi_file\n")
	print()
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
#len_chainA = int(sys.argv[4])

sigma_exp = 0.15
n_atoms = 0
i_atom = 0
item = ''
step = 0
cb_atoms_pdb = []
cb_res_pdb = []
pdb_chain_id = []
cb_atoms = []
ca_atoms = []
box = []
A = []
sigma = []
sigma_sq = []


out = open(output_file, 'w')

from Bio.PDB.PDBParser import PDBParser

p = PDBParser(PERMISSIVE=1)

def computeQ_inter():
	if len(cb_atoms)!=len(cb_atoms_pdb):
		print ("Error. Length mismatch!")
		print ("Pdb: ", len(cb_atoms_pdb), "trj: ", len(cb_atoms))
		exit()
	Q=0
	N = len(cb_atoms)
	count = 0
	for ia in range(0, len_chainA):
		for ja in range(len_chainA, N):
			rn = vabs(vector(cb_atoms_pdb[ia], cb_atoms_pdb[ja]))
			if rn > 9.5:
				continue
			count += 1
			r = vabs(vector(cb_atoms[ia], cb_atoms[ja]))
			dr = r - rn
			Q = Q + exp(-dr*dr/(2*sigma_sq[N/2]));
	Q = Q/count
	return Q

#push in all the CA atoms of pdb file
s = p.get_structure(struct_id, pdb_file)
chains = s[0].get_list()
print ("Number of chains:",  len(chains))
ichain = 0
for chain in chains:
	ichain = ichain + 1
	for res in chain:
		is_regular_res = res.has_id('CA') and res.has_id('O')
		res_id = res.get_id()[0]
	        if (res_id==' ' or res_id=='H_MSE' or res_id=='H_M3L' or res_id=='H_CAS' ) and is_regular_res:

			if res.get_resname() == 'GLY':
				cb_atoms_pdb.append(res['CA'].get_coord())
			else:
				cb_atoms_pdb.append(res['CB'].get_coord())

			cb_res_pdb.append(res.get_resname())
			pdb_chain_id.append(ichain)
print ("Number of residues in the PDB: ", len(cb_atoms_pdb))
len_chainA = pdb_chain_id.count(1)
#print (len_chainA)

for i in range(0, len(cb_atoms_pdb)+1):
	sigma.append( (1+i)**sigma_exp )
	sigma_sq.append(sigma[-1]*sigma[-1])

#Construct direct contact list among the chain: < 6.5 A
#out_contact=open("contact_list", 'w')
#for ia in range(0, len(cb_atoms_pdb)/2):
#	for ja in range(len(cb_atoms_pdb)/2, len(cb_atoms_pdb)):
#		#print (pdb_chain_id[ia], pdb_chain_id[ja])
#		#if  pdb_chain_id[ia] == pdb_chain_id[ja]: #only count the interface
#		#	continue
#		r_N=vabs(vector(cb_atoms_pdb[ia], cb_atoms_pdb[ja]));
#		#print (r_N)
#		if r_N < 6.5:
#			out_contact.write('1 '+str(ia+1)+' '+str(ja+1)+' '+cb_res_pdb[ia]+' '+cb_res_pdb[ja]+' '+str(round(r_N, 3))+'\n')
#		else:
#			if r_N < 9.5:
#				out_contact.write('2 '+str(ia+1)+' '+str(ja+1)+' '+cb_res_pdb[ia]+' '+cb_res_pdb[ja]+' '+str(round(r_N, 3))+'\n')
################################

#push in all the CA atoms of dump file
lfile = open(lammps_file)
for l in lfile:
	l = l.strip()
	if l[:5]=="ITEM:":
		item = l[6:]
	else:
		if item == "TIMESTEP":
			if len(cb_atoms)>0:
				q = computeQ_inter()
				#for key in q:
				out.write(str(round(q,3)))
				#out.write(' ')
				out.write('\n')
				n_atoms = len(cb_atoms)
			step = int(l)
			cb_atoms = []
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
				atom_CA = [x,y,z]
				#ca_atoms.append(atom)
			if desc=='C-Beta':
				atom = [x,y,z]
				cb_atoms.append(atom) 
			if desc=='H-Beta' :
				cb_atoms.append(atom_CA)
lfile.close()

if len(cb_atoms)>0:
	q = computeQ_inter()
	#for key in q:
	out.write(str(round(q,3)))
	#out.write(' ')
	out.write('\n')
	n_atoms = len(cb_atoms)

out.close()

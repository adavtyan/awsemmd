#!/usr/bin/python

# ----------------------------------------------------------------------
# Copyright (2010) Aram Davtyan and Garegin Papoian

# Papoian's Group, University of Maryland at Collage Park
# http://papoian.chem.umd.edu/

# Last Update: 03/04/2011
# ----------------------------------------------------------------------

import sys

#from Bio.PDB.PDBParser import PDBParser

atom_type = {'1' : 'C', '2' : 'N', '3' : 'O', '4' : 'C', '5' : 'H', '6' : 'C'}
atom_desc = {'1' : 'C-Alpha', '2' : 'N', '3' : 'O', '4' : 'C-Beta', '5' : 'H-Beta', '6' : 'C-Prime'}
PDB_type = {'1' : 'CA', '2' : 'N', '3' : 'O', '4' : 'CB', '5' : 'HB', '6' : 'C' }

def one2three(one_letter_code): 
    """ translate a protein sequence from 3 to 1 letter code"""
    
    code = {
	    "R": "ARG", "K": "LYS", "N": "ASN", "Q": "GLN", "E": "GLU",
	    "D": "ASP", "H": "HIS", "Y": "TYR", "W": "TRP", "S":"SER",
	    "T":"THR", "G":"GLY", "P":"PRO", "A":"ALA", "M":"MET",
	    "C":"CYS", "F":"PHE", "L":"LEU", "V":"VAL", "I":"ILE" }    
    index = code[one_letter_code]
    return index

class PDB_Atom:
	no = 0
	ty = ''
	res = 'UNK'
	chain = 'Z'
	res_no = 0
	x = 0.0
	y = 0.0
	z = 0.0
	atm = 'C'
	
	def __init__(self, no, ty, res, chain, res_no, x, y, z, atm):
		self.no = no
		self.ty = ty
		self.res = res
		self.chain = chain
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
		f.write(self.chain)
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

if len(sys.argv)!=3 and len(sys.argv)!=4:
	print "\n" + sys.argv[0] + " lammps_Input pdb_Output pdbID.seq [snapshot]\n"
	exit()

lammps_file = sys.argv[1]

output_file = ""
if len(sys.argv)>2: output_file = sys.argv[2]
seq_file = sys.argv[3]
fh = open(seq_file, 'r')

seqs_all = ''
ch_lens = [0]
nres_tot = 0
nch = 0 
for line in fh.readlines():
  seq = line.strip()
  seqs_all += seq
  ch_lens.append(ch_lens[-1]+len(seq))
  nres_tot += len(seq)
  nch += 1
fh.close()
print "Number of sequences:", nch
print "Total length:", nres_tot

ich = 1
ch_map = {}
for i in range(nres_tot):
  if ich<len(ch_lens) and i>=ch_lens[ich]: ich += 1
  ch_map[i+1] = chr(64+ich)

psf_file = output_file
if output_file[-4:]!=".pdb": output_file = output_file + ".pdb"
if psf_file[-4:]==".pdb": psf_file = psf_file[:-3] + "psf"
if psf_file[-4:]!=".psf": psf_file = psf_file + ".psf"

snapshot = -1
if len(sys.argv)>4: snapshot = int(sys.argv[4])

an = 0.4831806
bn = 0.7032820
cn = -0.1864262
ap = 0.4436538
bp = 0.2352006
cp = 0.3211455

n_atoms = 0
i_atom = 0
item = ''
step = 0
atoms = []
atoms2 = []
atoms3 = []
bonds = []
box = []
A = []

out = open(output_file, 'w')


def convertToPDB():
  ires = 1
  for ia in atoms2:
    if ia.desc == 'N': ires = ires + 1
    resname = one2three(seqs_all[ires-1])
    if not ch_map.has_key(ires):
      print "Error! atom list and sequance file size mismatch!\n"
      sys.exit()
    ch = ch_map[ires]
    atom = PDB_Atom(ia.No, PDB_type[ia.No_m], resname, ch, ires, ia.x, ia.y, ia.z, ia.ty)
    atoms3.append(atom)

def buildAllAtoms():
	index = 0
	last_Ca_index = -1
	last_O_index = -1
	Cp_index = -1
	NullVal = Atom(0, '', '6', 0.0, 0.0, 0.0, '')
#	atoms2 = []
	for i in range(0, len(atoms)):
		ia = atoms[i]
		index = index + 1
		if ia.desc == 'O': last_O_index = i
		if ia.desc == 'C-Alpha':
			if last_Ca_index != -1:
				Cai = atoms[last_Ca_index]
				Cai1 = ia
				Oi = atoms[last_O_index]
				
				nx = an*Cai.x + bn*Cai1.x + cn*Oi.x
				ny = an*Cai.y + bn*Cai1.y + cn*Oi.y
				nz = an*Cai.z + bn*Cai1.z + cn*Oi.z
				
				px = ap*Cai.x + bp*Cai1.x + cp*Oi.x
                                py = ap*Cai.y + bp*Cai1.y + cp*Oi.y
                                pz = ap*Cai.z + bp*Cai1.z + cp*Oi.z
				
				N = Atom(index, 'N', '2', nx, ny, nz, 'N')
				index = index + 1
				Cp = Atom(int(Cai.No) + 1, 'C', '6', px, py, pz, 'C-Prime')
#				Cp = Atom(index, 'C', '6', px, py, pz, 'C-Prime')
#				index = index + 1
				
				atoms2.append(N)
				atoms2.pop(Cp_index)
				atoms2.insert(Cp_index, Cp)
#				atoms2.append(Cp)
			last_Ca_index = i
		ia.No = index
                atoms2.append(ia)
		if ia.desc == 'C-Alpha':
			atoms2.append(NullVal)
			Cp_index = index
			index = index + 1
	if atoms2[Cp_index].No==0: atoms2.pop(Cp_index)
	for i in range(Cp_index, len(atoms2)):
		atoms2[i].No = atoms2[i].No - 1 

def buildBonds():
	N_index = -1
	Ca_index = -1
	Cp_index = -1
	O_index = -1
	Cb_index = -1
	Hb_index = -1
	ires = 1
	for i in range(0, len(atoms2)):
		ia = atoms2[i]
		if ia.desc == 'N':
			ires += 1
			if N_index!=-1 and Ca_index!=-1:
				bonds.append([N_index, Ca_index])
			if Ca_index!=-1 and Cp_index!=-1:
				bonds.append([Ca_index, Cp_index])
			if Cp_index!=-1 and O_index!=-1:
				bonds.append([Cp_index, O_index])
			if Ca_index!=-1 and Cb_index!=-1:
				bonds.append([Ca_index, Cb_index])
			if Ca_index!=-1 and Hb_index!=-1:
				bonds.append([Ca_index, Hb_index])
			N_index = i+1
			if Cp_index!=-1 and (ires<=1 or ch_map[ires]==ch_map[ires-1]):
				bonds.append([Cp_index, N_index])
			Ca_index = -1
			Cp_index = -1
			O_index = -1
			Cb_index = -1
			Hb_index = -1
		if ia.desc == 'C-Alpha': Ca_index = i+1
		if ia.desc == 'C-Beta': Cb_index = i+1
		if ia.desc == 'H-Beta': Hb_index = i+1
		if ia.desc == 'C-Prime': Cp_index = i+1
		if ia.desc == 'O': O_index = i+1
	if N_index!=-1 and Ca_index!=-1:
		bonds.append([N_index, Ca_index])
	if Ca_index!=-1 and Cb_index!=-1:
		bonds.append([Ca_index, Cb_index])
	if Ca_index!=-1 and Hb_index!=-1:
		bonds.append([Ca_index, Hb_index])

def print_atom_array():
	out.write("ITEM: TIMESTEP\n")
	out.write(str(step))
	out.write("\n")
	out.write("ITEM: NUMBER OF ATOMS\n")
	out.write(str(n_atoms))
	out.write("\n")
	out.write("ITEM: BOX BOUNDS\n")
	for ib in box:
		out.write(ib)
		out.write("\n")
	out.write("ITEM: ATOMS\n")
	for ia in atoms2:
		ia.write_(out)

def print_pdb():
	for ia in atoms3:
		ia.write_(out)
	out.write("END\n");

def print_psf():
	space8 = "        "
	psfout = open(psf_file,'w')
	psfout.write("PDF\n\n\t2 !NTITLE\n\n")
	psfout.write((space8+str(len(atoms3)))[-8:]+" !NATOM\n")
	for ia in atoms2:
		psfout.write((space8+str(ia.No))[-8:]+" PROT 1")
		psfout.write("    R00")
		psfout.write("  "+ia.ty)
		psfout.write("       1")
		psfout.write("          0            1           0\n")
	psfout.write("\n")
	psfout.write((space8+str(len(bonds)))[-8:]+" !NBOND")
	for i in range(0, len(bonds)):
		ib = bonds[i]
		if i%4==0: psfout.write("\n") 
		psfout.write((space8+str(ib[0]))[-8:])
		psfout.write((space8+str(ib[1]))[-8:])
	psfout.close()

nFrame = 0
found = False
lfile = open(lammps_file)
if snapshot<0:
	for l in lfile:
		l = l.strip()
		if l[:5]=="ITEM:":
			item = l[6:]
		else:
			if item == "TIMESTEP":
				if len(atoms)>0:
					buildAllAtoms()
					convertToPDB()
					n_atoms = len(atoms2)
					print_pdb()
				step = int(l)
				atoms = []
				atoms2 = []
				atoms3 = []
				box = []
				A = []
				nFrame = nFrame + 1
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
				atom = Atom(i_atom, atom_type[l[1]], l[1], x, y, z, desc)
				atoms.append(atom)
	
	if len(atoms)>0:
		buildAllAtoms()
		convertToPDB()
		n_atoms = len(atoms2)
		print_pdb()
		buildBonds()
		print_psf()
else:
	for l in lfile:
		l = l.strip()
		if l[:5]=="ITEM:":
			item = l[6:]
			if item == "TIMESTEP":
				if found: break
				elif nFrame==snapshot: found = True
				nFrame = nFrame + 1
		elif found:
			if item == "TIMESTEP":
				step = int(l)
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
				atom = Atom(i_atom, atom_type[l[1]], l[1], x, y, z, desc)
				atoms.append(atom)
	if len(atoms)>0:
		buildAllAtoms()
		convertToPDB()
		n_atoms = len(atoms2)
		print_pdb()
		buildBonds()
		print_psf()

lfile.close()
out.close()

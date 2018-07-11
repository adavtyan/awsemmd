#!/usr/bin/python

# ----------------------------------------------------------------------
# Copyright (2010) Aram Davtyan and Garegin Papoian

# Papoian's Group, University of Maryland at Collage Park
# http://papoian.chem.umd.edu/

# Last Update: 07/11/2018
# ----------------------------------------------------------------------

import sys
from math import sqrt, sin, cos

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
	No_m = 0
	res_no = 0
	x = 0.0
	y = 0.0
	z = 0.0
	desc = ''
	
	def __init__(self, No, ty, No_m, res_no, x, y, z, desc=''):
		self.No = No
		self.ty = ty
		self.No_m = No_m
		self.res_no = res_no
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


def vnorm(v):
	return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])

def vscale(scale, v):
  return [scale*v[0], scale*v[1], scale*v[2]]

def vdot(v1, v2):
  return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]

def vcross(v1, v2):
	return [v1[1]*v2[2] - v1[2]*v2[1], v1[2]*v2[0] - v1[0]*v2[2], v1[0]*v2[1] - v1[1]*v2[0]]

if len(sys.argv)!=4 and len(sys.argv)!=5:
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

b_terminal = False

an = 0.4831806
bn = 0.7032820
cn = -0.1864262
ap = 0.4436538
bp = 0.2352006
cp = 0.3211455

rNCa = 1.45808
rCaCp = 1.52469
rCpO = 1.23156
psi_NCaCp = 1.94437215835
psi_CaCpO = 2.10317155324
theta_NCaCpO = 2.4

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
  ires = 0
  for ia in atoms2:
    #if ia.desc == 'N': ires = ires + 1
    ires = ia.res_no
    resname = one2three(seqs_all[ires-1])
    if not ch_map.has_key(ires):
      print "Error! atom list and sequance file size mismatch!\n"
      sys.exit()
    ch = ch_map[ires]
    atom = PDB_Atom(ia.No, PDB_type[ia.No_m], resname, ch, ires, ia.x, ia.y, ia.z, ia.ty)
    atoms3.append(atom)

def buildAllAtoms(build_terminal_atoms=True, build_bonds=False):
	N_atom_map = {}
	Ca_atom_map = {}
	Cp_atom_map = {}
	O_atom_map = {}
	Cb_atom_map = {}

	# Map CA, CB, and O atoms
	res_no = 0
	for i in range(0, len(atoms)):
		ia = atoms[i]
		if ia.desc == 'C-Alpha':
			res_no += 1
			ia.res_no = res_no
			Ca_atom_map[res_no] = ia
		elif ia.desc == 'O':
			ia.res_no = res_no
			O_atom_map[res_no] = ia
		elif ia.desc == 'C-Beta' or ia.desc == 'H-Beta':
			ia.res_no = res_no
			Cb_atom_map[res_no] = ia

	# Recovering N and Cp atoms except for termianl residues
	for i in range(1, nres_tot):
		if not Ca_atom_map.has_key(i):
			print "Error! missing Ca atom in residue %d!\n" % i
			sys.exit()
		if not Ca_atom_map.has_key(i+1):
			print "Error! missing Ca atom in residue %d!\n" % (i+1)
			sys.exit()
		if not O_atom_map.has_key(i):
			print "Error! missing O atom in residue %d!\n" % i
			sys.exit()

		if ch_map[i]==ch_map[i+1]:
			Cai = Ca_atom_map[i]
			Cai1 = Ca_atom_map[i+1]
			Oi = O_atom_map[i]

			nx = an*Cai.x + bn*Cai1.x + cn*Oi.x
			ny = an*Cai.y + bn*Cai1.y + cn*Oi.y
			nz = an*Cai.z + bn*Cai1.z + cn*Oi.z
				
			px = ap*Cai.x + bp*Cai1.x + cp*Oi.x
			py = ap*Cai.y + bp*Cai1.y + cp*Oi.y
			pz = ap*Cai.z + bp*Cai1.z + cp*Oi.z

			N = Atom(0, 'N', '2', i+1, nx, ny, nz, 'N')
			Cp = Atom(0, 'C', '6', i, px, py, pz, 'C-Prime')

			N_atom_map[i+1] = N
			Cp_atom_map[i] = Cp

	# Recovering N and Cp atoms for terminal residues
	if build_terminal_atoms:
		for i in range(len(ch_lens)):
			j = ch_lens[i]
			if i!=len(ch_lens)-1:
				Ca = Ca_atom_map[j+1]
				Cp = Cp_atom_map[j+1]
				O = O_atom_map[j+1]

				r = rNCa
				psi = psi_NCaCp
				theta = -theta_NCaCpO

				v1 = [Cp.x - O.x, Cp.y - O.y, Cp.z - O.z]
				v2 = [Cp.x - Ca.x, Cp.y - Ca.y, Cp.z - Ca.z]

				mz = v2
				my = vcross(v1, mz)
				mx = vcross(my, mz)

				mx = vscale(r*sin(psi)*cos(theta)/vnorm(mx), mx)
				my = vscale(r*sin(psi)*sin(theta)/vnorm(my), my)
				mz = vscale(r*cos(psi)/vnorm(mz), mz)

				nx = Ca.x + mx[0] + my[0] + mz[0]
				ny = Ca.y + mx[1] + my[1] + mz[1]
				nz = Ca.z + mx[2] + my[2] + mz[2]

				N = Atom(0, 'N', '2', j+1, nx, ny, nz, 'N')

				if N_atom_map.has_key(j+1): print "Warrning: N atom already exists"

				N_atom_map[j+1] = N
			if i!=0:
				Ca = Ca_atom_map[j]
				N = N_atom_map[j]
				O = O_atom_map[j]

				sign = 1 # 1 or -1
				r1 = rNCa
				r2 = rCaCp
				r3 = rCpO
				psi1 = psi_NCaCp
				psi2 = psi_CaCpO

				xn = [N.x - Ca.x, N.y - Ca.y, N.z - Ca.z]
				xo = [O.x - Ca.x, O.y - Ca.y, O.z - Ca.z]

				ro = vnorm(xo)
				ro_sq = ro*ro
				rn_sq = vdot(xn, xn)
				r1o = vdot(xn, xo)
				A = r1*r2*cos(psi1)
				B = ro_sq + r2*r2 - r3*r3
				T1 = xn[1]*xn[1] + xn[2]*xn[2]
				T2 = xo[1]*xo[1] + xo[2]*xo[2]
				T3 = xn[1]*xo[1] + xn[2]*xo[2]
				T4 = xn[2]*xo[1] - xn[1]*xo[2]
				T5 = xn[0]*xo[2] - xn[2]*xo[0]
				T6 = xn[0]*xo[1] - xn[1]*xo[0]
 
				cprod = vcross(xn, xo)
				cprod_sq = vdot(cprod, cprod)

				D = cprod_sq*r2*r2 - A*A*ro_sq - 0.25*B*B*rn_sq + A*B*r1o

				# If D<0 reduce the angle psi1 by changing A and B
				if D<0:
					D = 0
					if abs(B)>2.0*ro*r2:
						B = 2.0*ro*r2
						A = B*r1o/ro_sq
					else:
						A = 0.5*( B*r1o + sqrt(cprod_sq*(4.0*ro_sq*r2*r2 - B*B)) ) / ro_sq
						if A>r1*r2: print "Warrning: A value too large"

				px = ( A*( -xo[0]*T3 + xn[0]*T2 ) + 0.5*B*(xo[0]*T1 - xn[0]*T3) + sign*T4*sqrt(D) ) / cprod_sq
				py = ( -A*xo[2] + 0.5*B*xn[2] + px*T5 ) / T4
				pz = ( A*xo[1] - 0.5*B*xn[1] - px*T6 ) / T4

				px += Ca.x
				py += Ca.y
				pz += Ca.z

				Cp = Atom(0, 'C', '6', j, px, py, pz, 'C-Prime')

				if Cp_atom_map.has_key(j): print "Warrning: Cp atom already exists"

				Cp_atom_map[j] = Cp

	# Assigning correct atom indexes and saving to atoms2 array
	index = 1
	for i in range(1, nres_tot+1):
		if not Ca_atom_map.has_key(i):
			print "Error! missing Ca atom in residue!\n" % i
			sys.exit()
		if not O_atom_map.has_key(i):
			print "Error! missing O atom in residue!\n" % i
			sys.exit()
		if not Cb_atom_map.has_key(i):
			print "Error! missing Cb or Hb atom in residue!\n" % i
			sys.exit()
		if not N_atom_map.has_key(i):
			if build_terminal_atoms or (i>1 and ch_map[i]==ch_map[i-1]):
				print "Error! missing N atom in residue!\n" % i
				sys.exit()
		if not Cp_atom_map.has_key(i):
			if build_terminal_atoms or (i<nres_tot and ch_map[i]==ch_map[i+1]):
				print "Error! missing Cp atom in residue!\n" % i
				sys.exit()

		add_atoms = []
		if N_atom_map.has_key(i): add_atoms.append(N_atom_map[i])
		add_atoms.append(Ca_atom_map[i])
		if Cp_atom_map.has_key(i): add_atoms.append(Cp_atom_map[i])
		add_atoms.append(O_atom_map[i])
		add_atoms.append(Cb_atom_map[i])

		for j in range(len(add_atoms)):
			ja = add_atoms[j]
			ja.No = index
			atoms2.append(ja)
			index += 1

	# Building bonds if necessary
	if build_bonds:
		for i in range(1, nres_tot+1):
			N_index = -1
			Ca_index = -1
			Cp_index = -1
			O_index = -1
			Cb_index = -1
			N1_index = -1

			Ca_index = Ca_atom_map[i].No
			Cb_index = Cb_atom_map[i].No
			O_index = O_atom_map[i].No
			if N_atom_map.has_key(i): N_index = N_atom_map[i].No
			if Cp_atom_map.has_key(i): Cp_index = Cp_atom_map[i].No
			if i<nres_tot and ch_map[i]==ch_map[i+1]: N1_index = N_atom_map[i+1].No

			if N_index!=-1:
				bonds.append([N_index, Ca_index])
			if Cp_index!=-1:
				bonds.append([Ca_index, Cp_index])
				bonds.append([Cp_index, O_index])
			bonds.append([Ca_index, Cb_index])
			if Cp_index!=-1 and N1_index!=-1:
				bonds.append([Ca_index, N1_index])

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
	psfout.write("PSF\n\n\t2 !NTITLE\n\n")
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
					buildAllAtoms(b_terminal)
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
				atom = Atom(i_atom, atom_type[l[1]], l[1], 0, x, y, z, desc)
				atoms.append(atom)
	
	if len(atoms)>0:
		buildAllAtoms(b_terminal, build_bonds=True)
		convertToPDB()
		n_atoms = len(atoms2)
		print_pdb()
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
				atom = Atom(i_atom, atom_type[l[1]], l[1], 0, x, y, z, desc)
				atoms.append(atom)
	if len(atoms)>0:
		buildAllAtoms(b_terminal, build_bonds=True)
		convertToPDB()
		n_atoms = len(atoms2)
		print_pdb()
		print_psf()

lfile.close()
out.close()

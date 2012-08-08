#Last updated 6/29/12, 1:10 pm

import random
import sys
import string
from math import *

o=open("data.dna",'w')

class Atom:
	index = 0
	mol_no = 0
	res_no = 0
	type = 0
	charge = 0.0
	x = 0.0
	y = 0.0
	z = 0.0
	def __init__(self, id, imol, ires, ty, q, xx, yy, zz):
		self.index = id
		self.mol_no = imol
		self.res_no = ires
		self.type = ty
		self.charge = q
		self.x = xx
		self.y = yy
		self.z = zz

class Bond:
	index = 0
	type = 0
	atomid1 = 0
	atomid2 = 0
	def __init__(self, id, ty, a1, a2):
		self.index = id
		self.type = ty
		self.atomid1 = a1
		self.atomid2 = a2

class Angle:
	index = 0
	type = 0
	atomid1 = 0
	atomid2 = 0
	atomid3 = 0
	def __init__(self, id, ty, a1, a2, a3):
		self.index = id
		self.type = ty
		self.atomid1 = a1
		self.atomid2 = a2
		self.atomid3 = a3

#box sizes
xstart, xend = -400.0, 400.0
ystart, yend = -400.0, 400.0
zstart, zend = -400.0, 400.0
#ion parameters
r0 = 10
mM = 5
cation_charge = 1
anion_charge = -1

residue_number=0
type=1
charge=-1.0
lch = 0
chain_number = 0

n_atoms = 0
n_bonds = 0
n_angles = 0
n_atom_types = 0
n_bond_types = 0
n_angle_types = 0

atoms = []
bonds = []
angles = []
ions = []
chain_lengths = []
starting_points = [0]

excluded_pairs = []

fene_ty = {1 : 10, 2 : 8, 3 : 6, 4 : 4, 5 : 2, 6 : 12, 7 : 3, 8 : 5, 9 : 7, 10 :9 , 11 : 11}
ion_type = {cation_charge : 2, anion_charge: 3}

f=open("sequences",'r')

sequences=[]
for line in f.readlines():
	line=line.strip()
	sequences.append(line)

chain_lengths = []
for i in range(len(sequences)):
	chain_lengths.append(len(sequences[i]))

# DNA coordinates
x = 0.0
y = 0.0
z = 0.0
bead_type = 1
for ich in range(0,len(sequences)/2):	#molecule 0 and 1
	lch1 = len(sequences[ich*2])
	lch2 = len(sequences[ich*2+1])
	if lch1 != lch2:
			print "dna mismatch"
	molecule = ich+1		#molecule 1 and 2
	xx = (11.065*ich*2)
	for c in range(2):		#chain 0 and 1
		chain_number = chain_number+1	#chain 1 2 3 4
		s = c*2-1
		for bp in range(len(sequences[chain_number-1])):	
			charge = -1.0
			basepair = bp+1		#residue 1-4, 1-4, 1-8, 1-8
			residue_number = residue_number+1	#"atom" 1-24
			x = xx+(s*(11.065/2)*cos(bp*35*pi/180)); x = "%.6f" % x 
			y = s*(11.065/2)*sin(bp*35*pi/180); y = "%.6f" % y
			z = (11.065/2)*(bp*35*pi/180) - 360; z = "%.6f" % z
			if bp == 0:
				charge = 0.0
			atom=Atom(residue_number, molecule, basepair, bead_type, charge, x, y, z)
			atoms.append(atom)
			bead_type += 1
		starting_points.append(starting_points[chain_number-1]+chain_lengths[chain_number-1])

#ion parameters
atom_number = residue_number
basepair = 0
number_anions =int( mM*(10**-3)*1000*(10**-30)*(xend-xstart)*(yend-ystart)*(zend-zstart)*(6.02*(10**23)))
number_cations = number_anions+len(atoms)-2
for i in range(number_cations):
	ions.append(cation_charge)
for i in range(number_anions):
	ions.append(anion_charge)
#print ions
#print len(ions)
print "Salt concentration", mM, "mM"
print "Number of anions:", number_anions
print "Number of cations:", number_cations

for ion in ions:
	g = True
	for i in range(20):
		x = random.uniform(xstart,xend); x = float("%.6f" % x)
		y = random.uniform(ystart,yend); y = float("%.6f" % y)
		z = random.uniform(zstart,zend); z = float("%.6f" % z)
		g = True
		#check if overlaps with anything
		for ia in atoms:
			xd = x-float(ia.x)
			yd = y-float(ia.y)
			zd = z-float(ia.z)
			r = sqrt(xd*xd+yd*yd+zd*zd)
			if r < r0: g = False
		if g: break
	if not g:
		print "Warning: Cannot create distant ion"
	#add ion in atoms array	
	atom_number = atom_number+1
	molecule = molecule+1
	atom = Atom(atom_number, molecule, basepair, ion_type[ions[ion]]+bead_type-2, float(ions[ion]), x, y, z)
	atoms.append(atom)

ibond = 1
for ich in range(0,len(chain_lengths)/2):
	lch1 = chain_lengths[ich*2]
	lch2 = chain_lengths[ich*2+1]
	if lch1 != lch2:
		print "dna mismatch"
		sys.exit()
	istp1 = starting_points[ich*2]
	istp2 = starting_points[ich*2+1]
	for i in range(lch1):
		atom1=istp1+i+1
		atom2=istp2+i+1
		if atom1 < istp1+lch1:
			type = 1
			bonds.append(Bond(ibond,type,atom1,atom1+1))
			ibond=ibond+1
			excluded_pairs.append([atom1,atom1+1])
		if atom2 < istp2+lch2:
			type = 1
			bonds.append(Bond(ibond,type,atom2,atom2+1))
			ibond=ibond+1
			excluded_pairs.append([atom2,atom2+1])
		excluded_pairs.append([atom1,atom2])
		for j in range(-5,6):
			atom3=atom2+j
			if atom3 <= istp2+lch2 and atom3 > istp2:
				type=fene_ty[j+6]
				bonds.append(Bond(ibond,type,atom1,atom3))
				ibond=ibond+1
						
iangle=1
for ich in range(0,len(chain_lengths)/2):
	lch1= (chain_lengths[ich*2])
	lch2= (chain_lengths[ich*2+1])
	if lch1 != lch2:
		print "dna mismatch"
		sys.exit()
	istp1 = starting_points[ich*2]
	istp2 = starting_points[ich*2+1]
	for i in range(lch1):
		atom1=istp1+i+1
		atom2=istp2+i+1
		if atom1 <= istp2-2:
			type=1
			angles.append(Angle(iangle,type,atom1,atom1+1,atom1+2))
			iangle = iangle+1
		if atom2 <= istp2+lch2-2:
			type=1
			angles.append(Angle(iangle,type,atom2,atom2+1,atom2+2))
			iangle = iangle+1

n_atoms=len(atoms)
n_bonds=len(bonds)
n_angles=len(angles)
n_atom_types = bead_type+1
n_bond_types = 12
n_angle_types = 1

o.write("LAMMPS protein data file"+"\n"+"\n")
o.write("\t"+str(n_atoms)+"   atoms"+"\n")
o.write("\t"+str(n_bonds)+"   bonds"+"\n")
o.write("\t"+str(n_angles)+"   angles"+"\n")
o.write("\t"+"0"+"   dihedrals"+"\n")
o.write("\t"+"0"+"   impropers"+"\n"+"\n")
o.write("\t"+str(n_atom_types)+"   atom types"+"\n")
o.write("\t"+str(n_bond_types)+"   bond types"+"\n")
o.write("\t"+str(n_angle_types)+"   angle types"+"\n")
o.write("\t"+"0"+"   dihedral types"+"\n")
o.write("\t"+"0"+"   improper types"+"\n"+"\n")
o.write(str(xstart)+"  "+str(xend)+"   xlo"+" xhi"+"\n")
o.write(str(ystart)+"  "+str(yend)+"   ylo"+" yhi"+"\n")
o.write(str(zstart)+"  "+str(zend)+"   zlo"+" zhi"+"\n\n")
o.write("BondBond Coeffs"+"\n\n")
o.write("1   0 0 0"+"\n\n")
o.write("BondAngle Coeffs"+"\n\n")
o.write("1   0 0 0 0"+"\n\n")
o.write("Masses"+"\n\n")
o.write("\t"+"1*"+str(bead_type-1)+"   330"+"\n")
o.write("\t"+str(bead_type)+"   22.989769"+"\n")
o.write("\t"+str(bead_type+1)+"    35.453"+"\n")
o.write("\n")
o.write("Atoms\n\n")

for ia in atoms:
	o.write("\t"+str(ia.index)+"\t")	
	o.write(str(ia.mol_no)+"\t")
	o.write(str(ia.res_no)+"\t")
	o.write(str(ia.type)+"\t")
	o.write(str(ia.charge)+"\t")
	o.write(str(ia.x)+"\t")
	o.write(str(ia.y)+"\t")
	o.write(str(ia.z)+"\n")

o.write("\nBonds\n\n")
for ib in bonds:
	o.write("\t"+str(ib.index)+"\t")
	o.write(str(ib.type)+"\t")
	o.write(str(ib.atomid1)+"\t")
	o.write(str(ib.atomid2)+"\n")

o.write("\nAngles\n\n")
for ic in angles:
	o.write("\t"+str(ic.index)+"\t")
	o.write(str(ic.type)+"\t")
	o.write(str(ic.atomid1)+"\t")
	o.write(str(ic.atomid2)+"\t")
	o.write(str(ic.atomid3)+"\n")

o.close
f.close

excluded_pairs_string = ""
for ip in excluded_pairs:
	excluded_pairs_string += "pair_coeff "+str(ip[0])+" "+str(ip[1])+"\tnone\n"

replace_rules = [ ["``excluded_pairs",  excluded_pairs_string] ]
inp = open("pattern.dna.in")
inFile = inp.read()
inp.close()

for ir in replace_rules:
    inFile = inFile.replace(ir[0], ir[1])

out = open("dna.in",'w')
out.write(inFile)
out.close()

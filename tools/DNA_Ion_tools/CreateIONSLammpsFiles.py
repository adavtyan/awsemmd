#Last updated 6/29/12, 1:10 pm

import random
import sys
import string
from math import *

o=open("data.ions",'w')

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

#box sizes
xstart, xend = -50.0, 50.0
ystart, yend = -50.0, 50.0
zstart, zend = -50.0, 50.0
#ion parameters
r0 = 5
mM = 100
cation_charge = 1
anion_charge = -1

residue_number=0
type=1
charge=-1.0
lch = 0
chain_number = 0
molecule = 0 

n_atoms = 0
n_atom_types = 0

atoms = []
ions = []
chain_lengths = []

excluded_pairs = []

fene_ty = {1 : 10, 2 : 8, 3 : 6, 4 : 4, 5 : 2, 6 : 12, 7 : 3, 8 : 5, 9 : 7, 10 :9 , 11 : 11}
ion_type = {cation_charge : 2, anion_charge: 3}

bead_type = 1

#ion parameters
atom_number = residue_number
basepair = 0
number_anions =int( round(mM*(10**-3)*1000*(10**-30)*(xend-xstart)*(yend-ystart)*(zend-zstart)*(6.02*(10**23))))
number_cations = number_anions
for i in range(number_cations):
	ions.append(cation_charge)
for i in range(number_anions):
	ions.append(anion_charge)
#print (ions)
#print (len(ions))
print ("Salt concentration", mM, "mM")
print ("Number of anions:", number_anions)
print ("Number of cations:", number_cations)

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
		print ("Warning: Cannot create distant ion")
	#add ion in atoms array	
	atom_number = atom_number+1
	molecule = molecule+1
	atom = Atom(atom_number, molecule, basepair, ion_type[ions[ion]]+bead_type-2, float(ions[ion]), x, y, z)
	atoms.append(atom)

n_atoms=len(atoms)
n_atom_types = bead_type+1

o.write("LAMMPS protein data file"+"\n"+"\n")
o.write("\t"+str(n_atoms)+"   atoms"+"\n")
o.write("\t"+"0"+"   bonds"+"\n")
o.write("\t"+"0"+"   angles"+"\n")
o.write("\t"+"0"+"   dihedrals"+"\n")
o.write("\t"+"0"+"   impropers"+"\n"+"\n")
o.write("\t"+str(n_atom_types)+"   atom types"+"\n")
o.write("\t"+"0"+"   bond types"+"\n")
o.write("\t"+"0"+"   angle types"+"\n")
o.write("\t"+"0"+"   dihedral types"+"\n")
o.write("\t"+"0"+"   improper types"+"\n"+"\n")
o.write(str(xstart)+"  "+str(xend)+"   xlo"+" xhi"+"\n")
o.write(str(ystart)+"  "+str(yend)+"   ylo"+" yhi"+"\n")
o.write(str(zstart)+"  "+str(zend)+"   zlo"+" zhi"+"\n\n")
o.write("Masses"+"\n\n")
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

o.close

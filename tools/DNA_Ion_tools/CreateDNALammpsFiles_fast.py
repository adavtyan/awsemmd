#Last updated 6/29/12, 1:10 pm

import random
import sys
import string
from math import *
from progressbar import ProgressBar

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

def box_number(x):
	global Lx, Ly, Lz
	global xstart, ystart, zstart
	global ndiv

	p=[0,0,0]
	xx = [x[0]-xstart, x[1]-ystart, x[2]-zstart]
	a = [Lx/2, Ly/2, Lz/2]
	b = [0.0, 0.0, 0.0]
	for i in range(ndiv):
		for j in range(3):
			p[j]*=2
			if xx[j]>=a[j]+b[j]:
				p[j]+=1
				b[j]+=a[j]
			a[j]/=2
	return p

def add_ghost(ibox, atom):
        global b_ghosts
        global ndaxis
	global xlo, ylo, zlo
	global xhi, yhi, zhi

	g=False
	i0=0; i1=0; j0=0; j1=0; k0=0; k1=0
        if x-xlo[ibox]<r0:
		i0=-1
		g=True
        if xhi[ibox]-x<r0:
		i1=2
		g=True
	if y-ylo[ibox]<r0:
		j0=-1
		g=True
	if ylo[ibox]-y<r0:
		j1=2
		g=True
	if z-zlo[ibox]<r0:
		k0=-1
		g=True
	if zlo[ibox]-z<r0:
		k1=2
		g=True

	if g==True:
		pbox=[ibox/(ndaxis*ndaxis), ibox%(ndaxis*ndaxis)/ndaxis, ibox%(ndaxis*ndaxis)%ndaxis]
		for i in range(i0,i1):
			for j in range(j0,j1):
				for k in range(k0,k1):
					if i==0 and j==0 and k==0: continue
					ib=pbox[0]+i
					jb=pbox[1]+j
					kb=pbox[2]+k

					if ib<0: ib+=ndaxis
					if ib>=ndaxis: ib-=ndaxis
					if jb<0: jb+=ndaxis
					if jb>=ndaxis: jb-=ndaxis
					if kb<0: kb+=ndaxis
					if kb>=ndaxis: kb-=ndaxis
					
					b = (ib*ndaxis+jb)*ndaxis+kb
					b_ghosts[b].append(atom)

#box sizes
xstart, xend = -200.0, 200.0
ystart, yend = -200.0, 200.0
zstart, zend = -200.0, 200.0
Lx = xend-xstart
Ly = yend-ystart
Lz = zend-zstart
#ion parameters
r0 = 5
mM = 1000
cation_charge = 1
anion_charge = -1
#Algorithm parameters
ndiv=3 # number of divisions
nboxes=8**ndiv # total number of boxes
ndaxis=2**ndiv # number of divitions along each axes
b_atoms=[] # per box atom arrays
b_ghosts=[] # lists of ghosts atoms for each box
xlo=[]; ylo=[]; zlo=[]
xhi=[]; yhi=[]; zhi=[]
for i in range(nboxes):
	b_atoms.append([])
	b_ghosts.append([])
	xlo.append(0.0)
	ylo.append(0.0)
	zlo.append(0.0)
	xhi.append(0.0)
	yhi.append(0.0)
	zhi.append(0.0)
for i in range(ndaxis):
	for j in range(ndaxis):
		for k in range(ndaxis):
			ibox=(i*ndaxis+j)*ndaxis+k
			xlo[ibox] = xstart + i*Lx/ndaxis
			xhi[ibox] = xstart + (i+1)*Lx/ndaxis
			ylo[ibox] = ystart + j*Ly/ndaxis
			yhi[ibox] = ystart + (j+1)*Ly/ndaxis
			zlo[ibox] = zstart + k*Lz/ndaxis
			zhi[ibox] = zstart + (k+1)*Lz/ndaxis


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
n_dna_chains = 0

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
bead_type = 0
n_dna_chains = len(sequences)/2
for ich in range(n_dna_chains):	#molecule 0 and 1
	lch1 = len(sequences[ich*2])
	lch2 = len(sequences[ich*2+1])
	if lch1 != lch2:
			print ("dna mismatch")
	molecule = ich+1		#molecule 1 and 2
	xx = (11.065*ich*2)
	for c in range(2):		#chain 0 and 1
		chain_number = chain_number+1	#chain 1 2 3 4
		s = c*2-1
		for bp in range(len(sequences[chain_number-1])):	
			charge = -1.0
			basepair = bp+1		#residue 1-4, 1-4, 1-8, 1-8
			residue_number = residue_number+1	#"atom" 1-24
			x = xx+(s*(11.065/2)*cos(bp*35*pi/180))
			y = s*(11.065/2)*sin(bp*35*pi/180)
			z = (11.065/2)*(bp*35*pi/180) - 160
			if bp == 0:
				charge = 0.0
			bead_type += 1
			atom=Atom(residue_number, molecule, basepair, bead_type, charge, x, y, z)
			pbox=box_number([x,y,z])
			ibox=(pbox[0]*ndaxis+pbox[1])*ndaxis+pbox[2]
			b_atoms[ibox].append(atom)

			#add ion as a ghost if needed
			add_ghost(ibox, atom)
			
		starting_points.append(starting_points[chain_number-1]+chain_lengths[chain_number-1])

#ion parameters
atom_number = residue_number
basepair = 0
number_anions =int( round(mM*(10**-3)*1000*(10**-30)*(xend-xstart)*(yend-ystart)*(zend-zstart)*(6.02*(10**23))))
number_cations = number_anions+atom_number-2*n_dna_chains
for i in range(number_cations):
	ions.append(cation_charge)
for i in range(number_anions):
	ions.append(anion_charge)
#print (ions)
#print (len(ions))
print ("Salt concentration", mM, "mM")
print ("Number of anions:", number_anions)
print ("Number of cations:", number_cations)

j=0
fr=len(ions)/100
if fr<=0: fr=1
prb= ProgressBar('green', width=50, block='*', empty='o')
for ion in ions:
	ibox=random.randint(0,nboxes-1) # randomly choose a box
	g = True
	for i in range(20):
		x = random.uniform(xlo[ibox],xhi[ibox])
		y = random.uniform(ylo[ibox],yhi[ibox])
		z = random.uniform(zlo[ibox],zhi[ibox])
		g = True
		#check if overlaps with anything
		for ia in b_atoms[ibox]:
			xd = x-float(ia.x)
			yd = y-float(ia.y)
			zd = z-float(ia.z)
			r = sqrt(xd*xd+yd*yd+zd*zd)
			if r < r0: g = False
		if g: break
	if not g:
		print ("Warning: Cannot create distant ion\n\n")
	#add ion in atoms array	
	atom_number = atom_number+1
	molecule = molecule+1
	atom = Atom(atom_number, molecule, basepair, ion_type[ions[ion]]+bead_type-1, float(ions[ion]), x, y, z)
	b_atoms[ibox].append(atom)

	#add ion as a ghost if needed
	add_ghost(ibox, atom)

	j += 1
	if j%fr==0:
		p = round(100*float(j)/len(ions),2)
		prb.render(int(p), 'Adding ions\nNumber of tries %d' % i)
prb.render(100, 'Adding ions\nNumber of tries %d' % i)

#for ib in range(nboxes):
#	print ("N:", len(b_atoms[ib]))

# Sort atoms
natoms = atom_number
for i in range(natoms):
	atoms.append(None)
for ib in range(nboxes):
	for ia in b_atoms[ib]:
		i = ia.index - 1
		atoms[i] = ia
for i in range(natoms):
	if atoms[i] is None:
		print ("Atom", i+1, "missing")
		sys.exit()

# Add bonds
ibond = 1
for ich in range(0,len(chain_lengths)/2):
	lch1 = chain_lengths[ich*2]
	lch2 = chain_lengths[ich*2+1]
	if lch1 != lch2:
		print ("dna strand length mismatch")
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

# Add angles
iangle=1
for ich in range(0,len(chain_lengths)/2):
	lch1= (chain_lengths[ich*2])
	lch2= (chain_lengths[ich*2+1])
	if lch1 != lch2:
		print ("dna mismatch")
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
n_atom_types = bead_type+2
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
#o.write("Masses"+"\n\n")
#o.write("\t"+"1*"+str(bead_type-1)+"   330"+"\n")
#o.write("\t"+str(bead_type)+"   22.989769"+"\n")
#o.write("\t"+str(bead_type+1)+"    35.453"+"\n")
#o.write("\n")
o.write("Atoms\n\n")

for ia in atoms:
	o.write("\t"+str(ia.index)+"\t")	
	o.write(str(ia.mol_no)+"\t")
	o.write(str(ia.res_no)+"\t")
	o.write(str(ia.type)+"\t")
	o.write(str(ia.charge)+"\t")
	o.write(str(round(ia.x,6))+"\t")
	o.write(str(round(ia.y,6))+"\t")
	o.write(str(round(ia.z,6))+"\n")

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

masses_string = ""
masses_string += "mass		1*" + str(bead_type)  + "	330.0\n"
masses_string += "mass		" + str(bead_type+1)  + "	23.0\n"
masses_string += "mass		" + str(bead_type+2)  + "	35.4"

excluded_pairs_string = ""
for ip in excluded_pairs:
	excluded_pairs_string += "pair_coeff "+str(ip[0])+" "+str(ip[1])+"\tnone\n"

replace_rules = [ ["``excluded_pairs",  excluded_pairs_string], ["``masses", masses_string] ]
inp = open("pattern.dna.in")
inFile = inp.read()
inp.close()

for ir in replace_rules:
    inFile = inFile.replace(ir[0], ir[1])

out = open("dna.in",'w')
out.write(inFile)
out.close()

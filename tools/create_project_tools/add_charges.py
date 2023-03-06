#!/storage-home/k/kp24/Enthought/Canopy_64bit/User/bin/python

import sys
from pylab import *
from itertools import islice

import math
def roundup(x):
    return int(math.ceil(x / 3.0))


class Bead:
        def __init__(self,id):

                self.id = id
                self.mol =0 
                self.res = 0
		self.atype=0
                self.q = 0
                self.x = 0
                self.y = 0
                self.z = 0
	
	def read_coords(self,id,mol,res,atype,q,x,y,z):
		self.id=id
		self.mol=mol
		self.res=res
		self.atype=atype
		self.q=q
		self.x=x
		self.y=y
		self.z=z

class Bond:
        def __init__(self,id):

                self.id =id
                self.type =0
                self.a1 = 0
                self.a2 = 0

        def create_bond(self,id,type,a1,a2):
		self.id=id
		self.type=type
                self.a1 = a1
                self.a2 = a2
class Angle:
        def __init__(self,id):

                self.id =id
                self.type =0
                self.a1 = 0
                self.a2 = 0		
		self.a3 = 0

	def create_angle(self,id,type,a1,a2,a3):

                self.id=id
                self.type=type
                self.a1 = a1
                self.a2 = a2
		self.a3 = a3

class Dihedral:
        def __init__(self,id):

                self.id =id
                self.type =0
                self.a1 = 0
                self.a2 = 0
		self.a3 = 0
		self.a4 = 0

        def create_dihedral(self,id,type,a1,a2,a3,a4):
                self.id=id
                self.type=type
                self.a1 = a1
                self.a2 = a2
		self.a3 = a3
		self.a4 = a4


class Masses:
	def __init__(self,id):
                self.id = id
                self.mass =0
	def read_coords(self,id,mass):
                self.id = id
		self.mass=mass

class Coeff:
        def __init__(self,id):
                self.id = id
                self.a =0
		self.b =0
        def read_coords(self,id,a,b):
                self.id=id
		self.a =a
		self.b =b


	
	
if (len(sys.argv)!=3):
        print ("Incorrect arguments:  <prot.data file> <file.seq>")
        sys.exit(1)

file_prot=open(sys.argv[1],'r')
file_seq=sys.argv[2]
with open(file_seq) as f:
    seq = "".join(line.strip() for line in f)

count = 0
for line in file_prot:
        count+=1
        if not line.strip():
                continue
        else:
                line = line.strip().split()
                if len(line) > 1:
                        if line[1]=='atoms':
                                Natoms1=int(line[0])
                        elif line[1]=='bonds':
                                Nbonds1=int(line[0])
                        elif line[1]=='atom':
                                Natype1=int(line[0])
                        elif line[1]=='bond':
                                Nbtype1=int(line[0])
			elif line[1]=='angles':
                                angle1=int(line[0])
                        elif line[1]=='dihedrals':
                                dihedral1=int(line[0])
                        elif line[1]=='angle':
                                angletype1=int(line[0])
                        elif line[1]=='dihedral':
                                dihedraltype1=int(line[0])
			elif line[0] == 'Bond':
                        	Lcoeff1 = count+1
		elif line[0]=='Atoms':
                        Latom1 = count+1
                elif line[0]=='Bonds':
                        Lbond1 = count+1
		elif line[0]=='Masses':
                        Lmasses1 = count+1
	
################PROTEIN########################
with open(sys.argv[1]) as f:
    content1 = f.readlines()

bead1 = [Bead(i) for i in range(Natoms1)]
masses1 = [Masses(i) for i in range (Natype1) ]
coeff1 = [Coeff(i) for i in range (Nbtype1) ]
bond1 = [Bond(i) for i in range (Nbonds1) ]
for i in range(int(Natype1)):
        line=content1[Lmasses1+i].strip().split()
        masses1[i].read_coords(line[0],line[1])

for i in range(int(Nbtype1)):
        line=content1[Lcoeff1+i].strip().split()
        coeff1[i].read_coords(line[0],line[1],line[2])

for i in range(int(Nbonds1)):
        line=content1[Lbond1+i].strip().split()
        bond1[i].create_bond(line[0],line[1],line[2],line[3])

ia = -1
q=0
res=-1

for i in range(int(Natoms1)):
        line=content1[Latom1+i].strip().split()
	ia += 1
	if (ia+1)%3==0:
        	res+=1
        	if seq[res] in ('K','R'):
                	q=1.0
        	elif seq[res] in ('D','E'):
                	q=-1.0
        	else:
                	q=0.0
        bead1[i].read_coords(line[0],line[1],line[2],line[3],str(q),line[5],line[6],line[7])
	q=0.0

### Writing up to a file

output = open('data.charged', 'w')

output.write("LAMMPS protein data file with charges\n\n")
output.write('	%d atoms\n 	%d bonds\n 	%d angles\n 	%d dihedrals\n\t0 impropers\n\n 	%d atom types\n 	%d bond types\n 	%d angle types\n 	%d dihedral types\n\t0 improper types\n' %(Natoms1,Nbonds1,angle1,dihedral1,Natype1,Nbtype1,angletype1,dihedraltype1))
output.write('\n	0.0 500.0 xlo xhi \n	0.0 500.0 ylo yhi \n 	0.0 500.0 zlo zhi\n\n')

output.write('Masses\n\n')
for i in range(Natype1):
	output.write('	%d	%s\n' %(int(masses1[i].id),masses1[i].mass))

output.write('\nAtoms\n\n')
for i in range(Natoms1):
	output.write('	%s	%s	%s	%d	%s	%s	%s	%s   \n' % (bead1[i].id,bead1[i].mol,bead1[i].res,int(bead1[i].atype),bead1[i].q,bead1[i].x,bead1[i].y,bead1[i].z))
 	num_res = int(bead1[i].res)
        num_mol = int(bead1[i].mol)

output.write('\nBond Coeffs\n\n')
for i in range(Nbtype1):
	output.write('	%d	%s	%s\n ' %(int(coeff1[i].id),coeff1[i].a,coeff1[i].b))

output.write('\nBonds\n\n')
for i in range(Nbonds1):
        output.write('	%d	%d      %d   %d   \n' % (int(bond1[i].id),int(bond1[i].type),int(bond1[i].a1),int(bond1[i].a2)))

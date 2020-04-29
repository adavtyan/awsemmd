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


	
	
if (len(sys.argv)!=4):
        print "Incorrect arguments:  <prot.data file> <dna.data file> <file.seq>"
        sys.exit(1)

file_prot=open(sys.argv[1],'r')
file_DNA=open(sys.argv[2],'r')
flie_seq=open(sys.argv[3],'r')
#seq = 'SVSERPPYSYMAMIQFAINSTERKRMTLKDIYTWIEDHFPYFKHIAKPGWKNSIRHNLSLHDMFVRETSANGKVSFWTIHPSANRYLTLDSERPPYSYMAMIQFAINSTERKRMTLKDIYTWIEDHFPYFKHIAKPGWKNSIRHNLSLHDMFVRETSANGKVSFWTIHPSANRYLTLDQVFKPLD'
with open(sys.argv[3]) as f:
    seq = " ".join(line.strip() for line in f)

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
count=0
for line in file_DNA:
        count+=1
        if not line.strip():
                continue
        else:
                line = line.strip().split()
                if len(line) > 1:
                        if line[1]=='atoms':
                                Natoms2=int(line[0])
                        elif line[1]=='bonds':
                                Nbonds2=int(line[0])
                        elif line[1]=='atom':
                                Natype2=int(line[0])
                        elif line[1]=='bond':
                                Nbtype2=int(line[0])
			elif line[1]=='angles':
                                Nangles2=int(line[0])
			elif line[1]=='dihedrals':
                                Ndihedrals2=int(line[0])
			elif line[1]=='angle':
                                angletype2=int(line[0])
			elif line[1]=='dihedral':
                                dihedraltype2=int(line[0])
			elif line[0]=='Bond':
                        	Lbondco2 = count+1
			elif line[0]=='Angle':
                                Langleco2 = count+1
			elif line[0]=='Dihedral':
                                Ldihedralco2 = count+1
                elif line[0]=='Atoms':
                        Latom2 = count+1
                elif line[0]=='Bonds':
                        Lbond2 = count+1
		elif line[0]=='Masses':
                        Lmasses2 = count+1
		elif line[0]=='Angles':
			Langle2 = count+1
		elif line[0]=='Dihedrals':
                        Ldihedral2 = count+1

	
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

#####################DNA##########################
with open(sys.argv[2]) as f:
    content2 = f.readlines()

bead2 = [Bead(i) for i in range(Natoms2)]
masses2 = [Masses(i) for i in range (Natype2) ]
bond2 = [Bond(i) for i in range (Nbonds2) ]
angle2 = [Angle(i) for i in range (Nangles2) ]
dihedral2 = [Dihedral(i) for i in range (Ndihedrals2) ]
for i in range(int(Natype2)):
        line=content2[Lmasses2+i].strip().split()
        masses2[i].read_coords(line[0],line[1])

for i in range(int(Nbonds2)):
        line=content2[Lbond2+i].strip().split()
        bond2[i].create_bond(line[0],line[1],line[2],line[3])

for i in range(int(Nangles2)):
        line=content2[Langle2+i].strip().split()
        angle2[i].create_angle(line[0],line[1],line[2],line[3],line[4])

for i in range(int(Ndihedrals2)):
        line=content2[Ldihedral2+i].strip().split()
        dihedral2[i].create_dihedral(line[0],line[1],line[2],line[3],line[4],line[5])

for i in range(int(Natoms2)):
        line=content2[Latom2+i].strip().split()
        bead2[i].read_coords(line[0],line[1],line[0],line[2],line[3],line[4],line[5],line[6])


### Writing up to a file

output = open('data.merge', 'w')

output.write("LAMMPS data file for AWSEM & 3SPN CG model\n\n")
output.write('\t%d atoms\n\t%d bonds\n\t%d angles\n\t%d dihedrals\n\n\t%d atom types\n\t%d bond types\n\t%d angle types\n\t%d dihedral types\n' %(Natoms1+Natoms2,Nbonds1+Nbonds2,angle1+Nangles2,dihedral1+Ndihedrals2,Natype1+Natype2,Nbtype1+Nbtype2,angletype1+angletype2,dihedraltype1+dihedraltype2))
output.write('\n\t-500.0 500.0 xlo xhi\n\t-500.0 500.0 ylo yhi\n\t-500.0 500.0 zlo zhi\n\n')

output.write('Masses\n\n')
for i in range(Natype2):
	output.write('\t%s\t%s\n' %(masses2[i].id,masses2[i].mass))
for i in range(Natype1):
	output.write('\t%d\t%s\n' %(int(masses1[i].id)+Natype2,masses1[i].mass))

output.write('\nBond Coeffs\n\n %s' %(content2[Lbondco2]))
for i in range(Nbtype1):
	output.write('\t%d\tharmonic\t%s\t%s\n' %(int(coeff1[i].id)+1,coeff1[i].a,coeff1[i].b))

output.write('\nAngle Coeffs\n\n')
for i in range(angletype2):
	output.write('\t%s'%(content2[Langleco2+i]))

output.write('\nDihedral Coeffs\n\n')
for i in range(dihedraltype2):
        output.write('\t%s'%(content2[Ldihedralco2+i]))

output.write('\nAtoms\n\n')
for i in range(Natoms1):
	output.write('\t%s\t%s\t%s\t%d\t%s\t%s\t%s\t%s\n' % (bead1[i].id,bead1[i].mol,bead1[i].res,int(bead1[i].atype)+Natype2,bead1[i].q,bead1[i].x,bead1[i].y,bead1[i].z))
 	num_res = int(bead1[i].res)
        num_mol = int(bead1[i].mol)
k = 0
chain = 0
for i in range(Natoms2):
	if int(bead2[i].mol)!=chain:
		output.write('\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\n' % (int(bead2[i].id)+Natoms1,int(bead2[i].mol)+num_mol,roundup(int(bead2[i].res)-2*k)+num_res,int(bead2[i].atype),bead2[i].q,bead2[i].x,bead2[i].y,bead2[i].z))
#		output.write('\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\n' % (int(bead2[i+1].id)+Natoms1,int(bead2[i+1].mol)+num_mol,1+num_res,int(bead2[i+1].atype),bead2[i+1].q,bead2[i+1].x,bead2[i+1].y,bead2[i+1].z))
		chain = int(bead2[i].mol)
		num_res = num_res+1
		i = i+1
		k = k+1
	else:
		chain = int(bead2[i].mol)	
        	output.write('\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\n' % (int(bead2[i].id)+Natoms1,int(bead2[i].mol)+num_mol,roundup(int(bead2[i].res)-2*k)+num_res,int(bead2[i].atype),bead2[i].q,bead2[i].x,bead2[i].y,bead2[i].z))
#		num_res=roundup(int(bead2[i].res)-2*k)+num_res

output.write('\nBonds\n\n')
for i in range(Nbonds1):
        output.write('\t%d\t%d\t%d\t%d\n' % (int(bond1[i].id),int(bond1[i].type)+1,int(bond1[i].a1),int(bond1[i].a2)))
for i in range(Nbonds2):
        output.write('\t%d\t%s\t%d\t%d\n' % (int(bond2[i].id)+Nbonds1,bond2[i].type,int(bond2[i].a1)+Natoms1,int(bond2[i].a2)+Natoms1))
#        output.write('    %d      1	%d   %d   \n' % (int(bond1[i].id)+Nbonds2,int(bond1[i].a1),int(bond1[i].a2)))

output.write('\nAngles\n\n')
for i in range(Nangles2):
        output.write('\t%s\t%s\t%d\t%d\t%d\n' % (angle2[i].id,angle2[i].type,int(angle2[i].a1)+Natoms1,int(angle2[i].a2)+Natoms1,int(angle2[i].a3)+Natoms1))

output.write('\nDihedrals\n\n')
for i in range(Ndihedrals2):
        output.write('\t%s\t%s\t%d\t%d\t%d\t%d\n' % (dihedral2[i].id,dihedral2[i].type,int(dihedral2[i].a1)+Natoms1,int(dihedral2[i].a2)+Natoms1,int(dihedral2[i].a3)+Natoms1,int(dihedral2[i].a4)+Natoms1))


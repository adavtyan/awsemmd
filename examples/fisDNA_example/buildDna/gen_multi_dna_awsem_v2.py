#!/usr/bin/python
##########################################################################
# Author: Khoa Pham 
# Date: Jan 26 2017
# Description: Generate LAMMPS data file for multiple chains of DNA based
#              on the existing single DNA data file. In the new data file,
#              all DNAs are properly aligned according to the displacement
#              values specified in x, y, and z.
#              Generate LAMMPS list files for multiple chains of DNA
#              based on the existing single DNA list files. 
# Input: Data file of single DNA; list files for bond, angle, dihedral;
#	 number of DNAs; values of displacement in x, y, z
# Output: New data file and list files (with multiple DNAs) compatible 
#         with AWSEM;
# Usage: python gen_dna_awsem.py <dna.data file><bond><angle><dihedral>
#        <index shift><number_of_dna><x><y><z><output> 
# Example: python gen_dna_awsem.py bdna_curv_conf.in in00_bond.list \
#          in00_angl.list in00_dihe.list 2 50 0 0 2bdna_curv_conf.in
# Revision: 
#          v0: Tue 21 Feb 2017 11:52 AM by Min-Yeh (Victor), Tsai
#          v1: bugfix for atom indexing for angle and dihed list files
#              Mar 2, 2017 by Victor
#          v2: add index shift variable in order to correct indexing issue
#              in DNA list files (in the presence of protein)
#              May 18, 2017 by Victor
#########################################################################

import sys
from pylab import *
from itertools import islice

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

class Bond_list:
        def __init__(self,id):
                self.id=id
                self.id1 =0
                self.id2 =0
                self.r = 0
                self.k2 = 0
                self.k3 = 0
                self.k4 = 0

        def create_bondlist(self,id,id1,id2,r,k2,k3,k4):
                self.id = id
                self.id1 =id1
                self.id2 =id2
                self.r = r
                self.k2 = k2
                self.k3 = k3
                self.k4 = k4


class Angle_list:
        def __init__(self,id):
                self.id=id
                self.id1 =0
                self.id2 =0
                self.id3 =0
                self.k1 = 0
                self.k2 = 0

        def create_anglelist(self,id,id1,id2,id3,k1,k2):
                self.id = id
                self.id1 =id1
                self.id2 =id2
                self.id3 =id3
                self.k1 = k1
                self.k2 = k2

class Dihedral_list:
        def __init__(self,id):
                self.id=id
                self.id1 =0
                self.id2 =0
                self.id3 =0
                self.id4 =0
                self.k1 = 0
                self.k2 = 0
                self.k3 = 0
                self.k4 = 0
                self.k5 = 0

        def create_dihedrallist(self,id,id1,id2,id3,id4,k1,k2,k3,k4,k5):
                self.id = id
                self.id1 =id1
                self.id2 =id2
                self.id3 =id3
                self.id4 =id4
                self.k1 = k1
                self.k2 = k2
                self.k3 = k3
                self.k4 = k4
                self.k5 = k5

	
	
if (len(sys.argv)!=11):
        print "Incorrect arguments:  <dna.data file> <bond> <angle> <dihedral> <index shift> <number_of_dna> <x> <y> <z> <output>"
        sys.exit(1)
file_DNA=open(sys.argv[1],'r')
file_bond=open(sys.argv[2],'r')
file_angle=open(sys.argv[3],'r')
file_dihedral=open(sys.argv[4],'r')
num_atoms_prot = int(sys.argv[5])
num_dna = int(sys.argv[6])
x_coor = int(sys.argv[7])
y_coor = int(sys.argv[8])
z_coor = int(sys.argv[9])

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

	
#####################DNA##########################
with open(sys.argv[1]) as f:
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
output = open(sys.argv[10],'w')

output.write("LAMMPS data file for 3SPN CG model\n\n")
output.write('\t%d atoms\n\t%d bonds \n\t%d angles\n\t%d dihedrals\n\n\t%d atom types\n\t%d bond types\n\t%d angle types\n\t%d dihedral types\n' %(Natoms2*num_dna,Nbonds2*num_dna,Nangles2*num_dna,Ndihedrals2*num_dna,14,Nbtype2,angletype2,dihedraltype2))
output.write('\n\t-500.0 500.0 xlo xhi\n\t-500.0 500.0 ylo yhi\n \t-500.0 500.0 zlo zhi\n\n')

output.write('Masses\n\n')
for i in range(14):
	output.write('\t%s\t%s\n' %(masses2[i].id,masses2[i].mass))

output.write('\nBond Coeffs\n\n %s' %(content2[Lbondco2]))

output.write('\nAngle Coeffs\n\n')
for i in range(angletype2):
	output.write('\t%s'%(content2[Langleco2+i]))

output.write('\nDihedral Coeffs\n\n')
for i in range(dihedraltype2):
        output.write('\t%s'%(content2[Ldihedralco2+i]))
num_mol = 0
temp = 0
output.write('\nAtoms\n\n')
for j in range(num_dna):
	num_mol = temp + num_mol 
	for i in range(Natoms2):
        	output.write('\t%d\t%d\t%s\t%s\t%f\t%f\t%f\n' % (int(bead2[i].id)+j*Natoms2,int(bead2[i].mol)+num_mol,bead2[i].atype,bead2[i].q,float(bead2[i].x)+j*x_coor,float(bead2[i].y)+j*y_coor,float(bead2[i].z)+j*z_coor))
		temp = int(bead2[i].mol) 	

output.write('\nBonds\n\n')
for j in range(num_dna):
	for i in range(Nbonds2):
        	output.write('\t%d\t%s\t%d\t%d\n' % (int(bond2[i].id)+j*Nbonds2,bond2[i].type,int(bond2[i].a1)+j*Natoms2,int(bond2[i].a2)+j*Natoms2))


output.write('\nAngles\n\n')
for j in range(num_dna):
	for i in range(Nangles2):
        	output.write('\t%d\t%s\t%d\t%d\t%d\n' % (int(angle2[i].id)+j*Nangles2,angle2[i].type,int(angle2[i].a1)+j*Natoms2,int(angle2[i].a2)+j*Natoms2,int(angle2[i].a3)+j*Natoms2))

output.write('\nDihedrals\n\n')
for j in range(num_dna):
	for i in range(Ndihedrals2):
        	output.write('\t%d\t%s\t%d\t%d\t%d\t%d\n' % (int(dihedral2[i].id)+j*Ndihedrals2,dihedral2[i].type,int(dihedral2[i].a1)+j*Natoms2,int(dihedral2[i].a2)+j*Natoms2,int(dihedral2[i].a3)+j*Natoms2,int(dihedral2[i].a4)+j*Natoms2))

#####################   Bond    ##########################
countl=0
for line in file_bond:
        if not line.strip():
                continue
        else:
                countl = countl+1
with open(sys.argv[2]) as f:
    content2 = f.readlines()
bond_list2 = [Bond_list(i) for i in range (countl)]
for i in range(int(countl)):
        line=content2[i].strip().split()
        bond_list2[i].create_bondlist(i,line[0],line[1],line[2],line[3],line[4],line[5])

### Writing up bond_list 

output = open('new_bond.list', 'w')
for j in range(num_dna):
	for i in range(countl):
        	output.write('%d\t%d\t%s\t%s\t%s\t%s\n' % (int(bond_list2[i].id1)+j*Natoms2+num_atoms_prot,int(bond_list2[i].id2)+j*Natoms2+num_atoms_prot,bond_list2[i].r,bond_list2[i].k2,bond_list2[i].k3,bond_list2[i].k4))

###################     Angle   ############################
countl=0
for line in file_angle:
        if not line.strip():
                continue
        else:
                countl = countl+1
with open(sys.argv[3]) as f:
    content3 = f.readlines()
angle_list2 = [Angle_list(i) for i in range (countl)]
for i in range(int(countl)):
        line=content3[i].strip().split()
        angle_list2[i].create_anglelist(i,line[0],line[1],line[2],line[3],line[4])

### Writing up angle_list 

output = open('new_angle.list', 'w')
for j in range(num_dna):
	for i in range(countl):
        	output.write('%d\t%d\t%d\t%s\t%s\n' % (int(angle_list2[i].id1)+j*Natoms2+num_atoms_prot,int(angle_list2[i].id2)+j*Natoms2+num_atoms_prot,int(angle_list2[i].id3)+j*Natoms2+num_atoms_prot,angle_list2[i].k1,angle_list2[i].k2))

###################     Dihedral   ############################
countl=0
for line in file_dihedral:
        if not line.strip():
                continue
        else:
                countl = countl+1
with open(sys.argv[4]) as f:
    content4 = f.readlines()
dihedral_list2 = [Dihedral_list(i) for i in range (countl)]
for i in range(int(countl)):
        line=content4[i].strip().split()
        dihedral_list2[i].create_dihedrallist(i,line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7],line[8])

### Writing up dihedral_list 

output = open('new_dihedral.list', 'w')
for j in range(num_dna):
	for i in range(countl):
        	output.write('%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%s\n' % (int(dihedral_list2[i].id1)+j*Natoms2+num_atoms_prot,int(dihedral_list2[i].id2)+j*Natoms2+num_atoms_prot,int(dihedral_list2[i].id3)+j*Natoms2+num_atoms_prot,int(dihedral_list2[i].id4)+j*Natoms2+num_atoms_prot,dihedral_list2[i].k1,dihedral_list2[i].k2,dihedral_list2[i].k3,dihedral_list2[i].k4,dihedral_list2[i].k5))


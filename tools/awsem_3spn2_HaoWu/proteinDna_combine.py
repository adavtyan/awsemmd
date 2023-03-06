#!/usr/local/bin/python

#######################################################################

# Copyright (2017) Hao Wu
# Papoian Research Group
# University of Maryland, College Park
# http://papoian.chem.umd.edu/
# Modified from Davit Potoyan's original code
# Last update: Nov 20 2017

#----------------------------------------------------------------------

# Write data.merge for AWSEM and 3SPN2 simulation

#######################################################################

import sys
#from pylab import *
#from itertools import islice

import math
def roundup(x):
    return int(math.ceil(x / 3.0))

#######################################################################
# CLASSES
#######################################################################

# Bead
class Bead:

    def __init__(self, id):
        self.id = id
        self.mol = 0
        self.res = 0
        self.atype = 0
        self.q = 0
        self.x = 0
        self.y = 0
        self.z = 0

    def read_coords(self, id, mol, res, atype, q, x, y, z):
        self.id = id
        self.mol = mol
        self.res = res
        self.atype = atype
        self.q = q
        self.x = x
        self.y = y
        self.z = z

# Bond
class Bond:

    def __init__(self, id):
        self.id = id
        self.type = 0
        self.a1 = 0
        self.a2 = 0

    def create_bond(self, id, type, a1, a2):
        self.id = id
        self.type = type
        self.a1 = a1
        self.a2 = a2

# Angle
class Angle:

    def __init__(self, id):
        self.id = id
        self.type = 0
        self.a1 = 0
        self.a2 = 0
        self.a3 = 0

    def create_angle(self, id, type, a1, a2, a3):
        self.id = id
        self.type = type
        self.a1 = a1
        self.a2 = a2
        self.a3 = a3

# Dihedral
class Dihedral:

    def __init__(self, id):
        self.id = id
        self.type = 0
        self.a1 = 0
        self.a2 = 0
        self.a3 = 0
        self.a4 = 0

    def create_dihedral(self, id, type, a1, a2, a3, a4):
        self.id = id
        self.type = type
        self.a1 = a1
        self.a2 = a2
        self.a3 = a3
        self.a4 = a4

# Masses
class Masses:

    def __init__(self, id):
        self.id = id
        self.mass = 0

    def read_coords(self, id, mass):
        self.id = id
        self.mass = mass

# Coeff
class Coeff:

    def __init__(self, id):
        self.id = id
        self.a = 0
        self.b = 0

    def read_coords(self, id, a, b):
        self.id = id
        self.a = a
        self.b = b

#######################################################################
# READING INPUT FILES
#######################################################################

#----------------------------------------------------------------------
# read input data
#----------------------------------------------------------------------

if (len(sys.argv) != 4):
    print ("Incorrect arguments: <prot.data file> <dna.data file> <file.seq>")
    sys.exit(1)

file_prot = open(sys.argv[1], 'r')
file_DNA = open(sys.argv[2], 'r')
flie_seq = open(sys.argv[3], 'r')
#seq = 'SVSERPPYSYMAMIQFAINSTERKRMTLKDIYTWIEDHFPYFKHIAKPGWKNSIRHNLSLHDMFVRETSANGKVSFWTIHPSANRYLTLDSERPPYSYMAMIQFAINSTERKRMTLKDIYTWIEDHFPYFKHIAKPGWKNSIRHNLSLHDMFVRETSANGKVSFWTIHPSANRYLTLDQVFKPLD'
# this definition of sequence will assign wrong charge!
# with open(sys.argv[3]) as f:
#     seq = " ".join(line.strip() for line in f)
with open(sys.argv[3]) as f:
    seq = "".join(line.strip() for line in f)

#----------------------------------------------------------------------
# read protein data
#----------------------------------------------------------------------

count = 0
# read protein data file
for line in file_prot:
    count += 1
    if not line.strip():
        continue
    else:
        line = line.strip().split()
        # overall information & coefficients
        if len(line) > 1:
            if line[1] == 'atoms':
                Natoms1 = int(line[0])
            elif line[1] == 'bonds':
                Nbonds1 = int(line[0])
            elif line[1] == 'atom':
                Natype1 = int(line[0])
            elif line[1] == 'bond':
                Nbtype1 = int(line[0])
            elif line[1] == 'angles':
                Nangles1 = int(line[0])
            elif line[1] == 'dihedrals':
                Ndihedrals1 = int(line[0])
            elif line[1] == 'angle':
                angletype1 = int(line[0])
            elif line[1] == 'dihedral':
                dihedraltype1 = int(line[0])
            elif line[0] == 'Bond':
                Lcoeff1 = count + 1
        # atoms, bonds, masses, angles and dihedrals header
        elif line[0] == 'Atoms':
            Latom1 = count + 1
        elif line[0]=='Bonds':
            Lbond1 = count + 1
        elif line[0]=='Masses':
            Lmasses1 = count + 1

#----------------------------------------------------------------------
# read dna data
#----------------------------------------------------------------------
count = 0
for line in file_DNA:
    count += 1
    if not line.strip():
        continue
    else:
        line = line.strip().split()
        # overall information & coefficients
        if len(line) > 1:
            if line[1] == 'atoms':
                Natoms2 = int(line[0])
            elif line[1] == 'bonds':
                Nbonds2 = int(line[0])
            elif line[1] == 'atom':
                Natype2 = int(line[0])
            elif line[1] == 'bond':
                Nbtype2 = int(line[0])
            elif line[1] == 'angles':
                Nangles2 = int(line[0])
            elif line[1] == 'dihedrals':
                Ndihedrals2 = int(line[0])
            elif line[1] == 'angle':
                angletype2 = int(line[0])
            elif line[1] == 'dihedral':
                dihedraltype2 = int(line[0])
            elif line[0]=='Bond':
                Lbondco2 = count + 1
            elif line[0]=='Angle':
                Langleco2 = count + 1
            elif line[0]=='Dihedral':
                Ldihedralco2 = count + 1
        # atoms, bonds, masses, angles and dihedrals header
        elif line[0] == 'Atoms':
            Latom2 = count + 1
        elif line[0] == 'Bonds':
            Lbond2 = count + 1
        elif line[0] == 'Masses':
            Lmasses2 = count + 1
        elif line[0] == 'Angles':
            Langle2 = count + 1
        elif line[0] == 'Dihedrals':
            Ldihedral2 = count + 1

#----------------------------------------------------------------------
# assign protein values
#----------------------------------------------------------------------

with open(sys.argv[1]) as f:
    content1 = f.readlines()

bead1 = [Bead(i) for i in range(Natoms1)]
masses1 = [Masses(i) for i in range (Natype1)]
coeff1 = [Coeff(i) for i in range (Nbtype1)]
bond1 = [Bond(i) for i in range (Nbonds1)]

# assign masses
for i in range(int(Natype1)):
    line = content1[Lmasses1 + i].strip().split()
    masses1[i].read_coords(line[0], line[1])

# assign bond coeff
for i in range(int(Nbtype1)):
    line = content1[Lcoeff1 + i].strip().split()
    coeff1[i].read_coords(line[0], line[1], line[2])

# assign bonds
for i in range(int(Nbonds1)):
    line = content1[Lbond1 + i].strip().split()
    bond1[i].create_bond(line[0], line[1], line[2], line[3])

# assign beads
ia = -1
q = 0.0
res = -1
for i in range(int(Natoms1)):
    line = content1[Latom1 + i].strip().split()
    ia += 1
    if (ia + 1) % 3 == 0:
        res += 1
        if seq[res] in ('K','R'):
            q = 1.0
        elif seq[res] in ('D','E'):
            q = -1.0
        else:
            q = 0.0
    bead1[i].read_coords(line[0], line[1], line[2], line[3], str(q), line[5], line[6], line[7])
    q = 0.0

#----------------------------------------------------------------------
# assign dna values
#----------------------------------------------------------------------

with open(sys.argv[2]) as f:
    content2 = f.readlines()

bead2 = [Bead(i) for i in range(Natoms2)]
masses2 = [Masses(i) for i in range (Natype2)]
bond2 = [Bond(i) for i in range (Nbonds2)]
angle2 = [Angle(i) for i in range (Nangles2)]
dihedral2 = [Dihedral(i) for i in range (Ndihedrals2)]

# assign masses
for i in range(int(Natype2)):
    line = content2[Lmasses2 + i].strip().split()
    masses2[i].read_coords(line[0], line[1])

# assign bonds
for i in range(int(Nbonds2)):
    line = content2[Lbond2 + i].strip().split()
    bond2[i].create_bond(line[0], line[1], line[2], line[3])

# assign angles
for i in range(int(Nangles2)):
    line = content2[Langle2 + i].strip().split()
    angle2[i].create_angle(line[0], line[1], line[2], line[3], line[4])

# assign dihedrals
for i in range(int(Ndihedrals2)):
    line = content2[Ldihedral2 + i].strip().split()
    dihedral2[i].create_dihedral(line[0], line[1], line[2], line[3], line[4], line[5])

# read correct residue id from crd file
with open('in00_conf.crd') as f:
    content_crd = f.readlines()

# assign beads
L_crd = 2 # read crd from the third line
for i in range(int(Natoms2)):
    line = content2[Latom2 + i].strip().split()
    line_crd = content_crd[L_crd + i].strip().split()
    bead2[i].read_coords(line[0], line[1], line_crd[1], line[2], line[3], line[4], line[5], line[6])

#######################################################################
# WRITE MERGED DATA
#######################################################################

output = open('data.merge', 'w')

# write general info
output.write("LAMMPS data file for AWSEM & 3SPN CG model\n\n")
output.write('	%d atoms \n 	%d bonds \n 	%d angles \n 	%d dihedrals \n\n 	%d atom types \n 	%d bond types \n 	%d angle types \n 	%d dihedral types \n'
             %(Natoms1 + Natoms2, Nbonds1 + Nbonds2, Nangles1 + Nangles2, Ndihedrals1 + Ndihedrals2,
               Natype1 + Natype2, Nbtype1 + Nbtype2, angletype1 + angletype2, dihedraltype1 + dihedraltype2))
output.write('\n	-200.0 200.0 xlo xhi \n	-200.0 200.0 ylo yhi \n 	-200.0 200.0 zlo zhi\n\n')

# write masses
output.write('Masses\n\n')
for i in range(Natype2):
    output.write('	%s	%s\n' %(masses2[i].id, masses2[i].mass))
for i in range(Natype1):
    output.write('	%d	%s\n' %(int(masses1[i].id) + Natype2, masses1[i].mass))

# write bond coeffs
output.write('\nBond Coeffs\n\n %s' %(content2[Lbondco2])) # 3SPN2 only has one bond coeff: list
for i in range(Nbtype1):
    output.write('	%d	harmonic	%s	%s\n ' %(int(coeff1[i].id) + 1, coeff1[i].a, coeff1[i].b))

# write angle coeffs
output.write('\nAngle Coeffs\n\n')
for i in range(angletype2):
    output.write('%s'%(content2[Langleco2 + i]))

# write dihedral coeffs
output.write('\nDihedral Coeffs\n\n')
for i in range(dihedraltype2):
    output.write('%s'%(content2[Ldihedralco2 + i]))

# write Atoms
output.write('\nAtoms\n\n')
for i in range(Natoms1):
    output.write('	%s      %s      %s      %d      %s      %s      %s      %s   \n'
                 % (bead1[i].id, bead1[i].mol, bead1[i].res, int(bead1[i].atype) + Natype2,
                    bead1[i].q, bead1[i].x, bead1[i].y, bead1[i].z))
    num_res = int(bead1[i].res)
    num_mol = int(bead1[i].mol)
for i in range(Natoms2):
    output.write('	%d      %d      %d      %d      %s      %s      %s      %s   \n'
                 % (int(bead2[i].id) + Natoms1, int(bead2[i].mol) + num_mol, int(bead2[i].res) + num_res,
                    int(bead2[i].atype), bead2[i].q, bead2[i].x, bead2[i].y, bead2[i].z))

# write Bonds
output.write('\nBonds\n\n')
for i in range(Nbonds2):
    output.write('	%s      %s      %d      %d      \n'
                 % (bond2[i].id, bond2[i].type, int(bond2[i].a1) + Natoms1, int(bond2[i].a2) + Natoms1))
for i in range(Nbonds1):
    output.write('	%d      %d      %d      %d      \n'
                 % (int(bond1[i].id) + Nbonds2, int(bond1[i].type) + 1, int(bond1[i].a1), int(bond1[i].a2)))

# write Angles
output.write('\nAngles\n\n')
for i in range(Nangles2):
    output.write('  %s      %s      %d      %d      %d      \n'
                 % (angle2[i].id, angle2[i].type, int(angle2[i].a1) + Natoms1,
                    int(angle2[i].a2) + Natoms1, int(angle2[i].a3) + Natoms1))

# write Angles
output.write('\nDihedrals\n\n')
for i in range(Ndihedrals2):
    output.write('  %s      %s      %d      %d      %d      %d      \n'
                 % (dihedral2[i].id, dihedral2[i].type, int(dihedral2[i].a1) + Natoms1,
                    int(dihedral2[i].a2) + Natoms1, int(dihedral2[i].a3) + Natoms1,
                    int(dihedral2[i].a4) + Natoms1))

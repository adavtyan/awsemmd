#!/usr/bin/python

# ----------------------------------------------------------------------
# Copyright (2010) Aram Davtyan and Garegin Papoian

# Papoian's Group, University of Maryland at Collage Park
# http://papoian.chem.umd.edu/

# Last Update: 03/04/2011
# ----------------------------------------------------------------------

import sys
from math import *
from VectorAlgebra import *

class Atom:
    ty = ''
    x = 0.0
    y = 0.0
    z = 0.0
    
    def __init__(self, No, ch, ty, x, y, z, desc=''):
        self.No = No
        self.ch = ch
        self.ty = ty
        self.x = x
        self.y = y
        self.z = z
        self.desc = desc

    def print_(self):
        print self.No, self.ch, self.ty , self.x, ',', self.y, ',', self.z, self.desc

    def write_(self, f):
        f.write(str(self.No))
        f.write('\t')
        f.write(str(self.ch))
        f.write('\t')
        f.write(self.ty)
        f.write('  ')
        f.write( ("               "+str(round(self.x,8)))[-15:] )
        f.write('\t')
        f.write( ("               "+str(round(self.y,8)))[-15:] )
        f.write('\t')
        f.write( ("               "+str(round(self.z,8)))[-15:] )
        if self.desc!='':
            f.write('\t\t')
            f.write( self.desc )
        f.write('\n')

def print_array(a):
    for ia in a:
        print ia

filename = ""
if len(sys.argv)>1: filename = sys.argv[1]

if filename=="":
    print "\nZ-MatrixToCoordinates.py input_file [output_file]\n"
    exit()

output_filename = ""
if len(sys.argv)>2: output_filename = sys.argv[2]

inp = open(filename)
atoms = []
no = 0
for l in inp:
    no += 1
    l = l.strip().split()
    if len(l)==0: continue
    if len(l)%2==1:
        atom = Atom(no, l[1], l[2], 0.0, 0.0, 0.0)
        di = 0
    else:
        atom = Atom(no, l[1], l[2], 0.0, 0.0, 0.0, l[3])
        di = 1
    if len(l)==5 or len(l)==6:
        i1 = int(l[4+di])-1
        atom.x = round(atoms[i1].x + float(l[3+di]), 12)
    elif len(l)==7 or len(l)==8:
        i1 = int(l[4+di])-1
        i2 = int(l[6+di])-1
        d = float(l[3+di])
        a = pi*float(l[5+di])/180
        dx = atoms[i1].x - atoms[i2].x
        dy = atoms[i1].y - atoms[i2].y
        d0 = sqrt(pow(dx,2) + pow(dy,2))
        atom.x = round(atoms[i1].x + d*(dy*sin(a) - dx*cos(a))/d0, 12)
        atom.y = round(atoms[i1].y - d*(dx*sin(a) + dy*cos(a))/d0, 12)
        atom.z = atoms[i1].z
    elif len(l)==9 or len(l)==10:
        i1 = int(l[4+di])-1
        i2 = int(l[6+di])-1
        i3 = int(l[8+di])-1
        d = float(l[3+di])
        a = pi*float(l[5+di])/180
        b = -pi*float(l[7+di])/180
        v1 = [atoms[i2].x-atoms[i3].x, atoms[i2].y-atoms[i3].y, atoms[i2].z-atoms[i3].z]
        v2 = [atoms[i2].x-atoms[i1].x, atoms[i2].y-atoms[i1].y, atoms[i2].z-atoms[i1].z]
        # must check if v1 & v2 are coliniar or not
        nz = v2
        ny = vcross_product(v1, nz)
        nx = vcross_product(ny, nz)

        nx = vproduct(d*sin(a)*cos(b)/vabs(nx), nx)
        ny = vproduct(d*sin(a)*sin(b)/vabs(ny), ny)
        nz = vproduct(d*cos(a)/vabs(nz), nz)

        atom.x = round(atoms[i1].x + nx[0] + ny[0] + nz[0], 12)
        atom.y = round(atoms[i1].y + nx[1] + ny[1] + nz[1], 12)
        atom.z = round(atoms[i1].z + nx[2] + ny[2] + nz[2], 12)
    
    atoms.append(atom)
inp.close()

if output_filename=="":
    for iAtm in atoms:
        iAtm.print_()
else:
    out = open(output_filename,'w')
    i = 0
    for iAtm in atoms:
        iAtm.write_(out)
        i+=1
    out.close()

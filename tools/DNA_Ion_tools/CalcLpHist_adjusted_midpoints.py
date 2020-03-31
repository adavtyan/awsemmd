#!/usr/bin/python

# ----------------------------------------------------------------------
# Copyright (2010) Aram Davtyan and Garegin Papoian

# Papoian's Group, University of Maryland at Collage Park
# http://papoian.chem.umd.edu/

# Last Update: 03/04/2011
# ----------------------------------------------------------------------

import sys
import numpy as np
from VectorAlgebra import *
from scipy.optimize import curve_fit


def func(x, a):
        return np.exp(-a * x)

atom_type = {'1' : 'C', '2' : 'N', '3' : 'O', '4' : 'C', '5' : 'H', '6' : 'C'}
#atom_desc = {'1' : 'C-Alpha', '2' : 'N', '3' : 'O', '4' : 'C-Beta', '5' : 'H-Beta', '6' : 'C-Prime'}
atom_desc = {'1' : 'DNA', '2' : "Na", '3' : 'Cl'}
PDB_type = {'1' : 'CA', '2' : 'N', '3' : 'O', '4' : 'CB', '5' : 'HB', '6' : 'C' }

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

if len(sys.argv)!=3:
	print "\nCalcQValue.py Input_file Output_file\n"
	exit()

input_file = sys.argv[1]

output_file = ""
if len(sys.argv)>2: output_file = sys.argv[2]


n_atoms = 0
i_atom = 0
item = ''
step = 0
ca_atoms = []
box = []
A = []

offset = 15
sep = 10
uij_sum = []
uij_cout = []

out = open(output_file, 'w')
hist = open("hist_cl.dat", 'w')

def calcMidPoint(p1, p2):
	half_period_x = (A[0][1] - A[0][0])/2
	half_period_y = (A[1][1] - A[1][0])/2
	half_period_z = (A[2][1] - A[2][0])/2
	mp = [0.0,0.0,0.0]

	mp[0] = (p1[0]+p2[0])/2
	mp[1] = (p1[1]+p2[1])/2
	mp[2] = (p1[2]+p2[2])/2
	if p1[0]-p2[0]>half_period_x:
		mp[0] += half_period_x
	elif p1[0]-p2[0]<-half_period_x:
		mp[0] -= half_period_x
	if p1[1]-p2[1]>half_period_y:
		mp[1] += half_period_y
	elif p1[1]-p2[1]<-half_period_y:
		mp[1] -= half_period_y
	if p1[2]-p2[2]>half_period_z:
		mp[2] += half_period_z
	elif p1[2]-p2[2]<-half_period_z:
		mp[2] -= half_period_z

	return mp

def calcDistVec(p1, p2):
	period_x = (A[0][1] - A[0][0])
	period_y = (A[1][1] - A[1][0])
	period_z = (A[2][1] - A[2][0])
	half_period_x = period_x/2
	half_period_y = period_y/2
	half_period_z = period_z/2

	v = [0.0, 0.0, 0.0]
	v[0] = p2[0]-p1[0]
	v[1] = p2[1]-p1[1]
	v[2] = p2[2]-p1[2]

	if v[0]<=-half_period_x:
		v[0] += period_x
	elif v[0]>half_period_x:
		v[0] -= period_x
	if v[1]<=-half_period_y:
		v[1] += period_y
	elif v[1]>half_period_y:
		v[1] -= period_y
	if v[2]<=-half_period_z:
		v[2] += period_z
	elif v[2]>half_period_z:
		v[2] -= period_z

	vm = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])

	return v

def calcUVec(p1, p2):
	period_x = (A[0][1] - A[0][0])
	period_y = (A[1][1] - A[1][0])
	period_z = (A[2][1] - A[2][0])
	half_period_x = period_x/2
	half_period_y = period_y/2
	half_period_z = period_z/2

	v = [0.0, 0.0, 0.0]
	v[0] = p2[0]-p1[0]
	v[1] = p2[1]-p1[1]
	v[2] = p2[2]-p1[2]

	if v[0]<=-half_period_x:
		v[0] += period_x
	elif v[0]>half_period_x:
		v[0] -= period_x
	if v[1]<=-half_period_y:
		v[1] += period_y
	elif v[1]>half_period_y:
		v[1] -= period_y
	if v[2]<=-half_period_z:
		v[2] += period_z
	elif v[2]>half_period_z:
		v[2] -= period_z

	vm = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])

	v[0] = v[0]/vm
	v[1] = v[1]/vm
	v[2] = v[2]/vm
	
	return v

def computeLp():
	global uij_sum
	global uij_count

	if len(ca_atoms)==0:
                print "Error. Empty snapshot"
                sys.exit()
	if len(ca_atoms)%2!=0:
		print "Odd number of atoms"
                sys.exit()

	N = len(ca_atoms)
	mid_points = []
	for i in range(N):
#		mid = calcMidPoint(ca_atoms[i],ca_atoms[i+N])
		mid_points.append(ca_atoms[i])

#	for i in range(N-1):
#		v12=vector(mid_points[i], mid_points[i+1])
#		d=sqrt(v12[0]*v12[0]+v12[1]*v12[1]+v12[2]*v12[2])
#		if d<100.0: hist.write(str(round(d,3))+"\n")

	for i in range(offset, N-offset):
		v12=vector(mid_points[i], mid_points[i+1])
		d=sqrt(v12[0]*v12[0]+v12[1]*v12[1]+v12[2]*v12[2])
		if d<100.0: hist.write(str(round(d,3))+"\n")
		
	U = []
	for i in range(offset, N-offset-sep+1):
		U.append(calcUVec(mid_points[i], mid_points[i+sep]))

	if len(uij_sum)==0:
		uij_sum = [0.0]*(len(U)-1)
		uij_count = [0]*(len(U)-1)

	for i in range(len(U)-1):
		ui = U[i]
		for j in range(i+1, len(U)):
			uj = U[j]
			uij = ui[0]*uj[0]+ui[1]*uj[1]+ui[2]*uj[2]
			uij_sum[j-i-1] += uij
			uij_count[j-i-1] += 1

step = 0	
lfile = open(input_file)
for l in lfile:
	l = l.strip()
	if l[:5]=="ITEM:":
		item = l[6:]
	else:
		if item == "TIMESTEP":
			if len(ca_atoms)>0:
				computeLp()
				n_atoms = len(ca_atoms)
			step = int(l)
			if step%100000==0 and step!=0 and len(uij_sum)>0:
				print step, uij_sum[0], uij_count[0]
			ca_atoms = []
			box = []
			A = []
		elif item == "NUMBER OF ATOMS":
			n_atoms = int(l)
		elif item[:10] == "BOX BOUNDS":
			box.append(l)
			l = l.split()
			A.append([float(l[0]), float(l[1])])
		elif item[:5] == "ATOMS":
			if step<1000000: continue
#			if step>100000000: break
			l = l.split()
			i_atom = l[0]
			x = float(l[2])
			y = float(l[3])
			z = float(l[4])
			x = (A[0][1] - A[0][0])*x + A[0][0]
			y = (A[1][1] - A[1][0])*y + A[1][0]
			z = (A[2][1] - A[2][0])*z + A[2][0]
#			desc = atom_desc[l[1]]
#			if desc=='DNA':
			if int(l[1])==2:
#				atom = Atom(i_atom, atom_type[l[1]], l[1], x, y, z, desc)
				atom = [x,y,z]
				ca_atoms.append(atom)
lfile.close()

if len(ca_atoms)>0:
	computeLp()
	n_atoms = len(ca_atoms)

x = []
y = []
for i in range(len(uij_sum)):
        avg = uij_sum[i]/uij_count[i]
        x.append(i+1)
        y.append(avg)
#       out.write(str(round(avg,3)))
#       out.write(' ')

x = np.array(x)
y = np.array(y)

np.savetxt('Uij_cl.dat', np.column_stack((x,y)), fmt='%g')

print x
print
print y
print

# Fit
popt, pcov = curve_fit(func, x, y)

std = 0.0
for i in range(len(x)):
        std += np.power(y[i] - func(x[i], *popt), 2)
std /= float(len(x)-1)
std = np.sqrt(std)

Lp = 1/popt[0]
print "Lp:", Lp
print "STD:", std

out.write(str(round(Lp,3)))

out.close()
hist.close()


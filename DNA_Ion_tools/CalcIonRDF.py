#!/usr/bin/python

# ----------------------------------------------------------------------
# Copyright (2012) Aram Davtyan and Garegin Papoian

# Papoian's Group, University of Maryland at Collage Park
# http://papoian.chem.umd.edu/

# Last Update: 08/01/2012
# ----------------------------------------------------------------------

import sys
from math import sqrt
from math import pi

class Atom:
	no = 0
	type = 0
	x = 0.0
	y = 0.0
	z = 0.0
	desc = ''

	def __init__(self, No, ty, x, y, z, desc):
		self.no = No
		self.type = ty
		self.x = x
		self.y = y
		self.z = z
		self.desc = desc

if len(sys.argv)!=3:
        print "\n"+sys.argv[0]+" Input_file Output_file\n"
        sys.exit()

infile = sys.argv[1]
outfile = sys.argv[2]

n_atoms = 0
n_frames = 0
i_atom = 0
item = ''
step = 0
atoms = []
box = []
A = []
V_Avg = 0.0

dr = 0.5
rmax = 20.0
nbins = int(round(rmax/dr))
hist_sums = []
hist_npairs = []
type_count = {}
hist_map = {}
hist2types_map = []
nhist = 0
offset = 10000000

def distance(i1, i2):
	p1 = [atoms[i1].x, atoms[i1].y, atoms[i1].z]
	p2 = [atoms[i2].x, atoms[i2].y, atoms[i2].z]

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

	return vm

def calc_rdfs():
	global n_frames
	global V_Avg
	global nhist

	n_frames += 1
	V_Avg += (A[0][1] - A[0][0])*(A[1][1] - A[1][0])*(A[2][1] - A[2][0])

	N = len(atoms)
	for i in range(N):
		ity = atoms[i].type
		if type_count.has_key(ity):
			type_count[ity] += 1
		else:
			type_count[ity] = 1
		for j in range(i+1,N):
			if i==j: continue
			jty = atoms[j].type
			d = distance(i,j)
			if d<=rmax:
				bin = int(d/dr)
				if not hist_map.has_key(ity):
					hist_map[ity] = {}
				if not hist_map.has_key(jty):
					hist_map[jty] = {}
				if not hist_map[ity].has_key(jty):
					hist_map[ity][jty] = nhist
					hist_map[jty][ity] = nhist
					hist_sums.append([0.0]*nbins)
					hist_npairs.append(0)
					hist2types_map.append([ity, jty])
					nhist += 1
				
				ih = hist_map[ity][jty]
				
				hist_sums[ih][bin] += 1
				hist_npairs[ih] += 1
				

atom_desc = {'1' : "Na+", '2' : "Cl-"}

lfile = open(infile)
for l in lfile:
        l = l.strip()
        if l[:5]=="ITEM:":
                item = l[6:]
        else:
                if item == "TIMESTEP":
                        if len(atoms)>0:
				if step%100000==0: print "Processing Frame # ", step
                                calc_rdfs()
                        step = int(l)
                        atoms = []
                        box = []
                        A = []
#			if step>=11000001: break
                elif item == "NUMBER OF ATOMS":
                        n_atoms = int(l)
                elif item[:10] == "BOX BOUNDS":
                        box.append(l)
                        l = l.split()
                        A.append([float(l[0]), float(l[1])])
                elif item[:5] == "ATOMS":
			if step<offset: continue 
                        l = l.split()
                        i_atom = l[0]
                        x = float(l[2])
                        y = float(l[3])
                        z = float(l[4])
                        x = (A[0][1] - A[0][0])*x + A[0][0]
                        y = (A[1][1] - A[1][0])*y + A[1][0]
                        z = (A[2][1] - A[2][0])*z + A[2][0]
                        desc = atom_desc[l[1]]
                        atom = Atom(i_atom, int(l[1]), x, y, z, desc)
                        atoms.append(atom)
lfile.close()

if len(atoms)>0:
	if step%100000==0: print "Processing Frame # ", step
	calc_rdfs()

if n_frames!=0: V_Avg /= n_frames

print
print "Number of Frames:", n_frames
print "Avarage Volume:", V_Avg
print "Number of Histogrames:", nhist
print "Atom Descriptions: ", atom_desc
print "Histogram Map:", hist_map
print "Number of pairs:", hist_npairs
print "Counted types:", type_count
print "Histogram:", hist_sums

for i in range(len(hist_sums)):
	itys = hist2types_map[i]
	for j in range(len(hist_sums[i])):
		dV = 4*pi*dr*dr*dr*(3*j*j+3*j+1)/3 # 4pi/3 * [r(i+1)^3 - r(i)^3] = 4pi/3 * dr^3 * [(i+1)^3 - i^3]
		hist_sums[i][j] *= n_frames*V_Avg/(type_count[itys[0]]*type_count[itys[1]]*dV)
		if itys[0]==itys[1]: hist_sums[i][j] *=2

print
print "Histogram:", hist_sums

out = open(outfile, 'w')
out.write('r\t')
for ih in range(len(hist_sums)):
	itys = hist2types_map[ih]
	atom_ty1 = atom_desc[str(itys[0])]
	atom_ty2 = atom_desc[str(itys[1])]
	out.write(atom_ty1)
	out.write('-')
	out.write(atom_ty2)
	if ih<len(hist_sums)-1: out.write('\t')
out.write('\n')
for i in range(nbins):
	r = i*dr+dr/2
	out.write(str(round(r,4)))
	out.write('\t')
	for ih in range(len(hist_sums)):
		out.write(str(round(hist_sums[ih][i],4)))
		if ih<len(hist_sums)-1: out.write('\t')
	if i<nbins-1: out.write('\n')
out.close()

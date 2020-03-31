#--------------------------------------------------------------------------------
# This routine calculates residue specific contact q values for a
# given dump.lammpstrj file.
#
# Notes: This routine assumes that the PDB file of interest contains a single chain.
#
# Created by Bobby Kim, adapted from code by Aram Davtyan.
#-------------------------------------------------------------------------------- 

import sys
import numpy
from VectorAlgebra import *
from Bio.PDB.PDBParser import PDBParser

if len(sys.argv) != 4 and len(sys.argv) != 5:
	print "\n" + str(sys.argv[0]) + " PDB_ID dump_file output_file [snapshot]\n"
	sys.exit()

struct_id = sys.argv[1]
filename = struct_id + ".pdb"
lammps_file = sys.argv[2] #'dump.lammpstrj'
output = sys.argv[3]
snapshot = -1
if len(sys.argv) > 4:
	snapshot = int(sys.argv[4])

p = PDBParser(PERMISSIVE=1)
s = p.get_structure(struct_id, filename)
chains = s[0].get_list()
atom_desc = {'1' : 'C-Alpha', '2' : 'N', '3' : 'O', '4' : 'C-Beta', '5' : 'H-Beta', '6' : 'C-Prime'}
cb_atoms = []
native_coords = []
native_distances = []
cutoff = 9.5
qis_total = []

#--------------------------------------------------------------------------------
def calc_dis(p1, p2):
	v = vector(p1, p2)
	return vabs(v)  

def compute_qis():
	if len(cb_atoms)!=len(native_distances):
		print "Error: length mismatch!"
		print "Pdb: ", len(native_distances), "trj: ", len(cb_atoms)
		exit()	
   	N = len(cb_atoms)
	qis = numpy.zeros([len(cb_atoms)])

	for i in range(0, N):
		for j in range(i+3, N):
			if native_contacts[i][j] == 1:
				r = vabs(vector(cb_atoms[i], cb_atoms[j]))
				if (r < 9.5):
					qis[i] += 1
					qis[j] += 1
	for i in range(0, N):
		qis[i] = qis[i]/norm[i]
	return qis
#--------------------------------------------------------------------------------

#import pdb file
for chain in chains:
	dis = []
	all_res = []
   	for res in chain:
		is_regular_res = res.has_id('CA') and res.has_id('O')
		res_id = res.get_id()[0]
		if (res.get_resname()=='GLY'):
			native_coords.append(res['CA'].get_coord())
		elif (res_id==' ' or res_id=='H_MSE' or res_id=='H_M3L' or res_id=='H_CAS') and is_regular_res:
			native_coords.append(res['CB'].get_coord())
		else:
			print 'ERROR: irregular residue at %s!' % res
			exit()

#calculate native distances from pdb file			
for i in range(0,len(native_coords)):
	native_distances.append([])
	for j in range(0, len(native_coords)):
		rij = calc_dis(native_coords[i],native_coords[j])
		native_distances[i].append(rij)

#create contact map
native_contacts = numpy.zeros([len(native_distances),len(native_distances)],int)
norm = numpy.zeros([len(native_distances)],int)
for i in range(0,len(native_coords)):
	for j in range(i+3,len(native_coords)):
		if (native_distances[i][j] < cutoff):
			native_contacts[i][j] = 1
			norm[i]+=1
			norm[j]+=1
		else:
			native_contacts[i][j] = 0
	if norm[i]==0:
		norm[i]=1
		print 'Warning: Residue %s has no native contacts. norm set to 1. The residue will always appear "unfolded".' % i

#read in dump.lammpstrj and calculate qi's
nFrame = 0
found = False
Qis_array = []
box = []
A = []
qis_array = []
lfile = open(lammps_file)
if snapshot < 0:
	for l in lfile:
		l = l.strip()
		if l[:5]=="ITEM:":
			item = l[6:]
		else:
			if item == "TIMESTEP":
				if len(cb_atoms)>0:
					qis_snapshot = compute_qis()
					qis_array.append(qis_snapshot)
					step = int(l)
				cb_atoms = []
				box = []
				A = []
			elif item[:10] == "BOX BOUNDS":
				box.append(l)
				l = l.split()
				A.append([float(l[0]), float(l[1])])
			elif item[:5] == "ATOMS":
				l = l.split()
				i_atom = l[0]
				x = float(l[2])
				y = float(l[3])
				z = float(l[4])
				x = (A[0][1] - A[0][0])*x + A[0][0]
				y = (A[1][1] - A[1][0])*y + A[1][0]
				z = (A[2][1] - A[2][0])*z + A[2][0]
				desc = atom_desc[l[1]]
				if desc=='C-Alpha':
					ca_atom = [x,y,z]
				elif desc=='C-Beta':
					cb_atom = [x,y,z]
					cb_atoms.append(cb_atom)
				elif desc=='H-Beta':
					cb_atoms.append(ca_atom)

        # calculate qis for last snapshot
	if len(cb_atoms)>0:
		qis_snapshot = compute_qis()
		qis_array.append(qis_snapshot)
		qis_total.append(qis_array)

else:
	for l in lfile:
		l = l.strip()
		if l[:5]=="ITEM:":
			item = l[6:]
			if item == "TIMESTEP":
				if found: break
				elif nFrame==snapshot: found = True
				nFrame = nFrame + 1
		elif found:			
			if item == "TIMESTEP":
				step = int(l)
			elif item[:10] == "BOX BOUNDS":
				box.append(l)
				l = l.split()
				A.append([float(l[0]), float(l[1])])
			elif item[:5] == "ATOMS":
				l = l.split()
				i_atom = l[0]
				x = float(l[2])
				y = float(l[3])
				z = float(l[4])
				x = (A[0][1] - A[0][0])*x + A[0][0]
				y = (A[1][1] - A[1][0])*y + A[1][0]
				z = (A[2][1] - A[2][0])*z + A[2][0]
				desc = atom_desc[l[1]]
				if desc=='C-Alpha':
					ca_atom = [x,y,z]
				elif desc=='C-Beta':
					cb_atom = [x,y,z]
					cb_atoms.append(cb_atom)
				elif desc=='H-Beta':
					cb_atoms.append(ca_atom)
        # calculate qis for last snapshot
	if len(cb_atoms)>0:
		qis_snapshot = compute_qis()
		qis_array.append(qis_snapshot)
		qis_total.append(qis_array)
		
lfile.close()


snapshotfile = open(output, 'w')

for i in range(0,len(qis_array)):
	for j in range(0,len(qis_array[i])):
		snapshotfile.write('%.2f' % qis_array[i][j])
		snapshotfile.write(' ')
	snapshotfile.write('\n')

snapshotfile.close()

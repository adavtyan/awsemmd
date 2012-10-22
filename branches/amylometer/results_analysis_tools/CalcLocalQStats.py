#--------------------------------------------------------------------------------
# This routine calculates residue q contacts given a list of dump.lammpstrj files.  
# The residue q contacts for a particular snapshot are then binned into 1 of 50 bins
# in the range [0,1] based on the global Qw.
# Within each bin, the average and variance is then calculated for each residue q contact.
#
# A file, 'directory_list', containing the full paths of the directories where the dump.lammpstrj and qw.dat
# files are located must be included.  qw.dat is a file containing a column of Qw's for each dump.lammpstrj
# snapshot excluding the first dump.lammpstrj snapshot (the zeroth snapshot).
#
# Notes: This routine assumes that the PDB file of interest contains a single chain.
#
#-------------------------------------------------------------------------------- 

import sys
import numpy
from VectorAlgebra import *
from pylab import *
from Bio.PDB.PDBParser import PDBParser

if len(sys.argv)!=4:
	print "\nCalcQValue.py PDB_ID directory_list output_file \n"
	exit()

struct_id = sys.argv[1]
filename = struct_id + ".pdb"
directory_list = sys.argv[2]
output = sys.argv[3]

p = PDBParser(PERMISSIVE=1)
s = p.get_structure(struct_id, filename)
chains = s[0].get_list()
atom_desc = {'1' : 'C-Alpha', '2' : 'N', '3' : 'O', '4' : 'C-Beta', '5' : 'H-Beta', '6' : 'C-Prime'}
lammps_file = 'dump.lammpstrj'
qw_file = 'qw.dat'
directories = []
cb_atoms = []
native_coords = []
native_distances = []
sigma = []
sigma_sq = []
sigma_exp = 0.15
cutoff = 9.5
Qw_data = []
qis_total = []
num_bin = 50
qis_new_bin = []
bin_spacing = 1.0/num_bin
average_array = []

#initialize qis_new_bin list
for i in range(0,num_bin):
	qis_new_bin.append([])

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
	qis = zeros([len(cb_atoms)])

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
native_contacts = zeros([len(native_distances),len(native_distances)],int)
norm = zeros([len(native_distances)],int)
for i in range(0,len(native_coords)):
	for j in range(i+3,len(native_coords)):
		if (native_distances[i][j] < cutoff):
			native_contacts[i][j] = 1
			norm[i]+=1
			norm[j]+=1
		else:
			native_contacts[i][j] = 0
	if norm[i]==0:
		print 'res %s has no native contacts' % i
		exit()

#calculate sigma table
for i in range(0,len(native_distances)+1):
	sigma.append((1+i)**sigma_exp)
	sigma_sq.append(sigma[-1]*sigma[-1])

#read in 'directories_list'
directories_file = open(directory_list,'r')
for line in directories_file:
	line=line.split()
	directories.append(line)

#read in qw.dat files
counter = 0
for path in directories:
	path = path[0] + '/' + qw_file
	qfile = open(path)
	Qw_data.append([])
	for l in qfile:
		l = l.strip()
		Qw_data[counter].append(l)
   	counter += 1
qfile.close()

#open mean and variance output files
out_mean = output + 'mean'
out_var = output + 'var'
mean_file = open(out_mean, 'w')
var_file = open(out_var, 'w')

#read in dump.lammpstrj and calculate qi's
for path in directories:
	qis_array = []
	path = path[0] + '/' + lammps_file
	print path
	lfile = open(path)
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
	lfile.close()

#calculate qis for last snapshot
	if len(cb_atoms)>0:
		qis_snapshot = compute_qis()
		qis_array.append(qis_snapshot)

#reset cb_atoms for next dump.lammpstrj file 
	cb_atoms = []	

	qis_total.append(qis_array)

#append snapshot qi data to appropriate bin
for i in range (0,len(Qw_data)):
	for j in range(0,len(Qw_data[i])):
		index = int(floor(float(Qw_data[i][j])/bin_spacing))
		qis_new_bin[index].append(qis_total[i][j+1])

#calculate average
for i in range(0,num_bin):
	average_array.append([])
	for j in range(0,len(native_distances)):
		sum = 0
		for k in range(0,len(qis_new_bin[i])):
			sum += qis_new_bin[i][k][j]
		if (len(qis_new_bin[i])!=0):
			average = sum/len(qis_new_bin[i])
		else:
			average = 0
		average_array[i].append(average)
		mean_file.write(str(round(average,4)))
		mean_file.write(' ')
	mean_file.write('\n')	
							
#calculate variance
for i in range(0,num_bin):
	for j in range(0,len(native_distances)):
		sum = 0
		for k in range(0,len(qis_new_bin[i])):
			sum += pow(qis_new_bin[i][k][j]-average_array[i][j],2)
		if (len(qis_new_bin[i])!=0):
			var = sum/len(qis_new_bin[i])
	   	else:
			var = 0
		var_file.write(str(round(var,7)))
	   	var_file.write(' ')
	var_file.write('\n')

mean_file.close()
var_file.close()

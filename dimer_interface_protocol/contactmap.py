import matplotlib.pyplot as plt
import sys
import numpy as np

from Bio.PDB import Vector

def add_frame():
	global infile
	global chain_boundary
	global end
	global average_img
	chain1 = []
	chain2 = []
	line = []
	while line[:11] !='ITEM: ATOMS':
		line = next(infile)
	line = next(infile)
	while line[:5] not in ['ITEM:', '\n']:
		line_split = line.split()
		if line_split[1] != '1':
			pass
		else:
			if int(line_split[0]) <=chain_boundary:
				atom = Vector(line_split[-3:])
				atom = Vector((atom._ar * np.array(400.0))+np.array(-200.0))
				chain1.append(atom)
			else:
				atom = Vector(line_split[-3:])
				atom = Vector((atom._ar * np.array(400.0))+np.array(-200.0))
				chain2.append(atom)
				
		try: line = next(infile)
		except StopIteration:
			end = True
			break
	#if line == '\n':
	#	end = True
	for i, atom1 in enumerate(chain1):
		row = []
		for j, atom2 in enumerate(chain2):
			d = (atom1-atom2).norm()
			average_img[i][j]+=d

def skip_frame():
	global infile
	global chain_boundary
	global end
	global average_img
	chain1 = []
	chain2 = []
	line = []
	while line[:11] !='ITEM: ATOMS':
		line = next(infile)
	line = next(infile)
	while line[:5] not in ['ITEM:', '\n']:
		try: line = next(infile)
		except StopIteration:
			end = True
			break
	#if line == '\n':
	#	end = True



def contactmap_draw(parametersobject):
	pd = parametersobject.parameterdic
	dd = parametersobject.deriveddic
	skip = pd['Skip_initial_frames']
	residues1 = dd['first_chain_length']
	residues2 = dd['second_chain_length']
	chain_boundary = dd['first_chain_max_id']
	f_data = open("contactmap.dat", "r")
	contact_min = pd['Contact_map_min_distance']
	contact_max = pd['Contact_map_max_distance']
	
	#average_img = [x[:] for x in [[0.0] * residues2] * residues1]
	average_img = []
	
	for i in range(residues1):
		line = next(f_data).strip()
		line_split = line.split()
		average_img.append([])
		for j in range(residues2):
			average_img[i].append(float(line_split[j]))
	f_data.close()
	plt.imshow(average_img, origin = 'lower', vmin = contact_min, vmax = contact_max, interpolation ='nearest', cmap = 'rainbow')
	plt.colorbar()
	plt.savefig('contactmap.png', bbox_inches = 'tight')
	plt.savefig('contactmap.svg', bbox_inches = 'tight')
	plt.close()
	plt.imshow(average_img, origin = 'lower', interpolation ='nearest', cmap = 'rainbow')
	plt.colorbar()
	plt.savefig('contactmap_free.png', bbox_inches = 'tight')
	plt.savefig('contactmap_free.svg', bbox_inches = 'tight')
	plt.close()
	
	
def contactmap_getdata(parametersobject):
	global infile
	global chain_boundary
	global end
	global average_img
	
	pd = parametersobject.parameterdic
	dd = parametersobject.deriveddic
	skip = pd['Skip_initial_frames']
	residues1 = dd['first_chain_length']
	residues2 = dd['second_chain_length']
	chain_boundary = dd['first_chain_max_id']
	f_data = open("contactmap.dat", "w+")
	orientations = pd['Number_of_orientations']
	
	reading_atoms = False
	end = False
	average_img = [x[:] for x in [[0.0] * residues2] * residues1]
	
	count = 0
	for i in range(1,orientations+1):
		end = False
		infile = open("r_"+str(i).zfill(3)+".lammpstrj", "r")
		frames = 0
		while not end:
			frames+=1
			if frames>skip:
				add_frame()
				count+=1
			else:
				skip_frame()
			if frames%10 == 0:
				sys.stdout.write('%d frames processed in file %d out of %d.\r' % (frames, i, orientations) )
				sys.stdout.flush()
		print('%d frames processed in file %d out of %d.' % (frames, i, orientations))
		infile.close()
	f_data = open("contactmap.dat", "w+")
	for i in range(residues1):
		for j in range(residues2):
			average_img[i][j] /= count
			f_data.write("%.4f\t" %(average_img[i][j]))
		f_data.write("\n")
	f_data.close()
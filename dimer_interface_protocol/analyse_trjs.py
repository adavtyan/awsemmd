import sys, os
import Bio.PDB
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import Vector
from Bio.PDB import PDBIO
import numpy as np
from Bio.PDB.Vector import rotmat
from Bio.SVDSuperimposer import SVDSuperimposer

def analyse(input_file_name, refer_file_name, moved_chain_id, fixed_chain_id, r_moved_chain_id, r_fixed_chain_id, output_file, r_model_number = 0):

	structure = PDBParser(PERMISSIVE=1).get_structure('to_analyse', input_file_name)
	reference = PDBParser(PERMISSIVE=1).get_structure('reference', refer_file_name)

	r_chain_moved = reference[r_model_number][r_moved_chain_id]
	r_chain_fixed = reference[r_model_number][r_fixed_chain_id]

	theta = []
	phi = []
	theta_x = []
	theta_y = []
	theta_z = []
	d=[]

	for model_number, model in enumerate(structure):
		chain_moved = structure[model_number][moved_chain_id]
		chain_fixed = structure[model_number][fixed_chain_id]
		com_denominator = 0.0
		com_numerator = Vector(0,0,0)
		for atom in chain_moved.get_atoms():
			position = atom.get_vector()
			com_numerator += Vector(position._ar*np.array(atom.mass))
			com_denominator +=atom.mass
		
		moved_centre = com_numerator.__div__(com_denominator)
		com_denominator=0.0
		com_numerator = Vector(0,0,0)
		for atom in chain_fixed.get_atoms():
			position = atom.get_vector()
			com_numerator += Vector(position._ar*np.array(atom.mass))
			com_denominator +=atom.mass
		
		fixed_centre = com_numerator.__div__(com_denominator)
		com_denominator=0.0
		com_numerator = Vector(0,0,0)
		
		
		
		
		reference_set = np.asarray([coord for coord in [atom.get_coord() for atom in r_chain_fixed.get_atoms()]])
		coordinate_set = np.asarray([coord for coord in [atom.get_coord() for atom in chain_fixed.get_atoms()]])
		sup = SVDSuperimposer()
		sup.set(reference_set, coordinate_set)
		sup.run()
		R, V = sup.get_rotran()
		for atom in model.get_atoms():
			atom.transform(R, V)
		for atom in chain_moved.get_atoms():
			com_numerator += Vector((atom.get_vector())._ar*np.array(atom.mass))
			com_denominator +=atom.mass
		moved_centre = com_numerator.__div__(com_denominator)
		com_denominator=0.0
		com_numerator = Vector(0,0,0)
		
		for atom in chain_fixed.get_atoms():
			com_numerator += Vector((atom.get_vector())._ar*np.array(atom.mass))
			com_denominator +=atom.mass
		fixed_centre = com_numerator.__div__(com_denominator)
		if fixed_centre.norm() > 0.5:
			print("Fixed chain norm is "+str(fixed_centre.norm())+" in model "+str(model_number)+". Should have been at the origin. Check code...")
		com_denominator=0.0
		com_numerator = Vector(0,0,0)
		
		d.append((moved_centre - fixed_centre).norm())
		if moved_centre.norm() > 1e-6:
			theta.append(moved_centre.angle(Vector(0,0,1)))
			
			x = moved_centre._ar[0]
			y = moved_centre._ar[1]
			norm = np.sqrt(x*x + y*y)
			if norm > 1e-6:
				'''
				if y>0:
					phi.append(np.arccos(x/norm))
				else:
					phi.append(np.arccos(x/norm) + np.pi)
				'''
				phi.append(np.arctan2(y,x))
		else:
			theta.append(0.0)

		
		
		reference_set = np.asarray([coord for coord in [atom.get_coord() for atom in r_chain_moved.get_atoms()]])
		coordinate_set = np.asarray([coord for coord in [atom.get_coord() for atom in chain_moved.get_atoms()]])
		sup = SVDSuperimposer()
		sup.set(reference_set, coordinate_set)
		sup.run()
		R, V = sup.get_rotran()
		theta_x.append(np.arctan2(R[2][1], R[2][2]))
		theta_y.append(np.arctan2(-R[2][0], np.sqrt(R[2][1]*R[2][1]+R[2][2]*R[2][2])))
		theta_z.append(np.arctan2(R[1][0], R[0][0]))

	f_results = open(output_file, "w+")
	for frame in range(0,len(structure)):
		f_results.write(str(frame)+'\t'+str(d[frame])+'\t'+str(theta[frame])+'\t'+str(phi[frame])+'\t'+str(theta_x[frame])+'\t'+str(theta_y[frame])+'\t'+str(theta_z[frame])+'\n')
	f_results.close()
	
def analyse_trjs(parametersobject):
	pd = parametersobject.parameterdic
	dd = parametersobject.deriveddic
	Path_to_awsem = pd['Path_to_awsem']
	Number_of_orientations = pd['Number_of_orientations']
	name = pd['Initial_dimer_pdb'][:-4]
	for i in range(1,1+Number_of_orientations):
		print("Analysing trajectory number\t"+str(i))
		location = os.path.normpath(Path_to_awsem+'/dimer_interface_protocol/Build_pdb_with_models.py')
		os.system('python2 '+location+" r_"+str(i).zfill(3)+".lammpstrj t_"+str(i).zfill(3)+" "+name+"_recentred"+".seq")
		input_file_name = 't_'+str(i).zfill(3)+'.pdb'
		output_file = 'angles_'+str(i).zfill(3)+".txt"
		if dd['first_chain_is_bigger']:
			fixed_chain_id = 'A'
			r_fixed_chain_id = 'A'
			moved_chain_id = 'B'
			r_moved_chain_id = 'B'
		else:
			fixed_chain_id = 'B'
			r_fixed_chain_id = 'B'
			moved_chain_id = 'A'
			r_moved_chain_id = 'A'
			
		analyse(
			input_file_name = input_file_name,
			refer_file_name = 'refpdb.pdb',
			fixed_chain_id = fixed_chain_id,
			moved_chain_id = moved_chain_id,
			r_fixed_chain_id = r_fixed_chain_id,
			r_moved_chain_id = r_moved_chain_id,
			r_model_number = 0,
			output_file = output_file
		)

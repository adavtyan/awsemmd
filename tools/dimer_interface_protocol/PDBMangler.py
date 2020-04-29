import sys
import Bio.PDB
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import Vector
from Bio.PDB import PDBIO
import random
import numpy as np
from Bio.PDB.Vector import rotmat


def COM(object):
	com_n = Vector(0,0,0)
	com_d = 0.0
	for atom in object.get_atoms():
		position = atom.get_vector()
		com_n += Vector(position._ar*np.array(atom.mass))
		com_d += atom.mass
	com = com_n.__div__(com_d)
	return com
	


def generate_3d():
	#generate uniform random SO(3) matrix
	#from randrot package sourcecode
	#website https://github.com/qobilidop/randrot/blob/master/randrot/__init__.py

	x1, x2, x3 = np.random.rand(3)
	R = np.matrix([[np.cos(2 * np.pi * x1), np.sin(2 * np.pi * x1), 0],
                   [-np.sin(2 * np.pi * x1), np.cos(2 * np.pi * x1), 0],
                   [0, 0, 1]])
	v = np.matrix([[np.cos(2 * np.pi * x2) * np.sqrt(x3)],
                   [np.sin(2 * np.pi * x2) * np.sqrt(x3)],
                   [np.sqrt(1 - x3)]])
	H = np.eye(3) - 2 * v * v.T
	M = -H * R
	return np.array(M)


def create_random_pdb(separation_distance, move_chain_id, fix_chain_id, input_file_name, output_pdb_name, model_number = 0):
	results = {}
	structure = PDBParser(PERMISSIVE=1).get_structure('whatever', input_file_name)
	chain_moved = structure[model_number][move_chain_id]
	chain_fixed = structure[model_number][fix_chain_id]

	old_fixed_centre = COM(chain_fixed)
	old_moved_centre = COM(chain_moved)
	com_denominator=0.0
	com_numerator = Vector(0,0,0)
		
	for atom in chain_moved.get_atoms():
		position = atom.get_vector()
		atom.set_coord(position - old_fixed_centre)

	#first step is to move origin to the com of fixed_chain.
	#So far the atoms in the moved_chain have been relocated.

	for atom in chain_fixed.get_atoms():
		position = atom.get_vector()
		atom.set_coord(position - old_fixed_centre)
	#now fixed_chain has been relocated. All coordinates are now wrt com of fixed_chain
	
	moved_centre = old_moved_centre - old_fixed_centre
	fixed_centre = Vector(0,0,0)


	d = (old_fixed_centre - old_moved_centre).norm()
	results["1_Input_Separation"] = d
	results["1_Old_fixed_chain_com"]=old_fixed_centre
	results["1_Old_moved_chain_com"]= old_moved_centre
	results["0_Intended_Output_Separation"] = separation_distance

	R1 = generate_3d()
	R2 = generate_3d()

	max_distance = 0.0
	com_numerator = Vector(0,0,0)
	com_denominator=0.0
	
	#Now we scale the separation distance and also rotate the chain_moved
	for atom in chain_moved.get_atoms():
		position = atom.get_vector()
		a = moved_centre.normalized()._ar * np.array(separation_distance)
		atom.set_coord((position - moved_centre).left_multiply(R2) + Vector(a))
		max_distance = max(max_distance, (atom.get_vector().norm()))
		position = atom.get_vector()
		com_numerator += Vector(position._ar*np.array(atom.mass))
		com_denominator +=atom.mass

	final_moved_centre = com_numerator.__div__(com_denominator)

	com_denominator=0.0
	com_numerator = Vector(0,0,0)	

	#Now we rotate the chain_fixed
	for atom in chain_fixed.get_atoms():
		position = atom.get_vector()
		atom.set_coord(position.left_multiply(R1))
		max_distance = max(max_distance, (atom.get_vector().norm()))
		position = atom.get_vector()
		com_numerator += Vector(position._ar*np.array(atom.mass))
		com_denominator +=atom.mass

	final_fixed_centre = com_numerator.__div__(com_denominator)
	d = (final_fixed_centre - final_moved_centre).norm()

	w = PDBIO()
	w.set_structure(structure)
	w.save(output_pdb_name)
	results["2_Output_Separation"]=d
	results["2_fixed_chain_com"]=final_fixed_centre
	results["2_moved_chain_com"]=final_moved_centre
	results["Max_distance"]=max_distance
	return results

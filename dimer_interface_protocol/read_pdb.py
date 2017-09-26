from Bio.PDB.PDBParser import PDBParser
from PDBMangler import COM
from Bio.PDB import Vector
from Bio.PDB import PDBIO
import numpy as np
import os, sys

def read_pdb(parametersobject):
	actualstdout = sys.stdout
	sys.stdout = open(os.devnull, 'w')
	pdbname = parametersobject.parameterdic['Initial_dimer_pdb']
	Path_to_awsem = parametersobject.parameterdic['Path_to_awsem']
	Python2_command = parametersobject.parameterdic['Python2_command']
	name = pdbname[:-4]
	structure = PDBParser(PERMISSIVE=1).get_structure('init', pdbname)
	if len(structure) > 1:
		print("More than one model found in PDB. Using model 0 only.")
	if len(structure[0]) >2:
		print("More than two chains found in PDB. Exiting.")
		sys.exit(1)
	elif len(structure[0])<2:
		print("Less than two chains found in PDB. Exiting.")
		sys.exit(1)

	chainnames = []
	chain = structure[0].get_list()
	for c in chain:
		chainnames.append(c.id)
	if len(chain[0])<len(chain[1]):
		bigger = 1
		bigid = chain[1].id
		smaller = 0
		smallid = chain[0].id
		first_chain_is_bigger = False
	else:
		bigger = 0
		bigid = chain[0].id
		smaller = 1
		smallid = chain[1].id
		first_chain_is_bigger = True
	#now chain[0] is the bigger chain
	#next steps are
	#recentre
	#get .data file.
	#convert to lammpstrj
	#convert back to pdb
	#yeah, it's a bit ridiculous but I don't know how the weirdly mangled pdb is created by awsem and I need it to be exactly the same as the simulation.
	#remove useless files
	#write information e
	centre = COM(chain[bigger])
	for atom in structure[0].get_atoms():
		atom.set_coord(atom.coord - centre._ar)
	w = PDBIO()
	w.set_structure(structure)
	w.save(name+'_recentred.pdb')
	
	
	
	
	
	cwd = os.getcwd()
	
	directorynames = ['md_input', 'md_output', 'analysis', 'results_main', 'results_individual', 'pdb_trajectories']
	
	for d in directorynames:
		directory = os.path.normpath(cwd+'/'+d)
		try:
			os.makedirs(directory)
		except OSError as e:
			pass
	
	
	
	
	os.system(Python2_command+" "+Path_to_awsem+"/create_project_tools/PDBToCoordinates.py "+name+"_recentred "+name+"_recentred"+".coord")
	os.system(Python2_command+" "+Path_to_awsem+"/create_project_tools/CoordinatesToWorkLammpsDataFile.py "+name+"_recentred"+".coord "+name+"_recentred"+".data -b")
	os.system(Python2_command+" "+Path_to_awsem+"/frag_mem_tools/Pdb2Gro.py "+name+"_recentred "+" md_input/chain1.gro "+chain[0].id)
	os.system(Python2_command+" "+Path_to_awsem+"/frag_mem_tools/Pdb2Gro.py "+name+"_recentred "+" md_input/chain2.gro "+chain[1].id)
	f_data = open(name+"_recentred"+".data", "r")
	f_lammps = open(name+"_recentred"+".lammpstrj", "w+")
	f_lammps.write("ITEM: TIMESTEP\n0\nITEM: BOX BOUNDS ff ff ff\n")
	for _ in range(3):
		f_lammps.write("-2.0000000000000000e+02 2.0000000000000000e+02\n")
	f_lammps.write("ITEM: ATOMS id type xs ys zs\n")
	firstchain = True
	for linecount, line in enumerate(f_data):
		if linecount<28:
			continue
		linesplit = line.strip().split()
		if len(linesplit)<1:
			break
		x = (float(linesplit[5])+200)/400
		y = (float(linesplit[6])+200)/400
		z = (float(linesplit[7])+200)/400
		f_lammps.write(linesplit[0]+' '+linesplit[3])
		f_lammps.write(' %.9f %.9f %.9f\n' % (x, y, z))
		if firstchain:
			if linesplit[1]=='2':
				firstchain=False
				first_chain_max_id = int(linesplit[0])-1
				
	f_lammps.close()
	f_data.close()
	
	
	location = os.path.normpath(Path_to_awsem+"results_analysis_tools/BuildAllAtomsFromLammps_seq_multichain.py "+name+"_recentred"+".lammpstrj")
	os.system(Python2_command+" "+location+" refpdb "+name+"_recentred"+".seq")
	sys.stdout = actualstdout
	os.remove("refpdb.psf")
	os.remove(name+"_recentred"+".lammpstrj")
	os.remove(name+"_recentred"+".data")
	os.remove(name+"_recentred"+".coord")
	d = parametersobject.deriveddic
	d['first_chain'] = chainnames[0]
	d['second_chain'] = chainnames[1]
	d['first_chain_length'] = len(chain[0])
	d['second_chain_length'] = len(chain[1])
	d['bigger_chain'] = bigid
	d['smaller_chain'] = smallid
	d['first_chain_max_id'] = first_chain_max_id
	d['first_chain_is_bigger'] = first_chain_is_bigger
	parametersobject.save_derived()
	
	
	
	
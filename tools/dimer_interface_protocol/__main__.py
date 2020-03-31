import sys
import os
from shutil import copy

from Parameterclass import Parameterclass
from read_pdb import read_pdb
from create_project import create_project
from analyse_trjs import analyse_trjs
from get_energies import get_energies
from graph_angles import graph_angles
from contactmap import contactmap_draw
from contactmap import contactmap_getdata
from delete_extras import delete_extras
from cluster import *


def write_sample_parameters():
	cwd = os.getcwd()
	dest = os.path.normpath(cwd+'/sample_parameters.dat')
	dir_path = os.path.dirname(os.path.realpath(__file__))
	src = os.path.normpath(dir_path+'/sample_parameters.dat')
	copy(src, dest)
	
def execute1(parametersobject):
	read_pdb(parametersobject)
	create_project(parametersobject)
def execute2(parametersobject):
	parametersobject.read_derived()
	if not parametersobject.parameterdic['Replot_only']:
		analyse_trjs(parametersobject)
		if parametersobject.parameterdic['Plot_energy']:
			get_energies(parametersobject)
		contactmap_getdata(parametersobject)
		cluster(parametersobject)
	graph_angles(parametersobject)
	contactmap_draw(parametersobject)
	cluster_dotplot(parametersobject)


if len(sys.argv) < 2:
	print("Please specify the parameter file. eg: python dimer_interface_protocol.py parameters.dat\nTo write a sample parameters file in the current directory, type 1. This will create the file 'sample_parameters.dat', overwriting any existing file of the same name.\n")
	line = sys.stdin.readline().strip()
	if line == '1':
		write_sample_parameters()
		print("File written")
	else:
		print("File not written")
	sys.exit(0)
else:
	parameterfilename = sys.argv[1]
	parametersobject = Parameterclass()
	parametersobject.read_parameters(parameterfilename)
	'''except InputError:
		print("Parameters file is not in expected format")
		sys.exit(1)
	'''
if len(sys.argv) == 3:
	instruction = sys.argv[2]
else:
	print("What would you like to do? Enter the corresponding number.")
	print("1: Create project from the initial PDB")
	print("\t1.1 Read the initial PDB only")
	print("\t1.2 Create the project files only")
	print("2: Analyse trajectory files and create graphs")
	print("\t2.1 Analyse trajectory orientations only")
	print("\t2.2 Analyse trajectory energies only")
	print("\t2.3 Analyse trajectory contacts only")
	print("\t2.4 Analyse clusters only")
	print("\t2.5 Draw orientation (and energy) plots only")
	print("\t2.6 Draw contact map only")
	print("\t2.7 Draw cluster coded dot plot only")
	print("3: Remove PDB trajectories, keeping only lammps trajectories.")
	instruction = sys.stdin.readline().strip()
	
instructionlist = {'1':execute1, '2':execute2, '1.1':read_pdb, '1.2':create_project, '2.1':analyse_trjs, '2.2':get_energies, '2.3':contactmap_getdata, '2.4':cluster, '2.5':graph_angles, '2.6':contactmap_draw, '2.7':cluster_dotplot, '3':delete_extras}
if instruction not in instructionlist:
	print("Invalid instructions entered.")
	sys.exit(0)

elif instruction in ['1.1', '1']:
	instructionlist[instruction](parametersobject)
else:
	parametersobject.read_derived()
	instructionlist[instruction](parametersobject)


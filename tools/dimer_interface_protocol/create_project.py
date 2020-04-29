import numpy as np
import os
from shutil import copy
from PDBMangler import create_random_pdb
from cd import cd

def create_project(parametersobject):
	dd = parametersobject.deriveddic
	pd = parametersobject.parameterdic
	name = pd['Initial_dimer_pdb']
	steps = pd['Timesteps']
	hours = pd['Simulation_hours']
	minutes = pd['Simulation_minutes']
	number_of_orientations = pd['Number_of_orientations']
	separation_distance = pd['COM_separation']
	Boundary_margin = pd['Boundary_margin']
	Initial_dimer_pdb = pd['Initial_dimer_pdb']
	move_chain_id = dd['smaller_chain']
	fix_chain_id = dd['bigger_chain']
	Jobname = pd['Jobname']
	Path_to_awsem = pd['Path_to_awsem']
	Path_to_lmp_serial = pd['Path_to_lmp_serial']
	Python2_command = pd['Python2_command']
	model_number = 0
	
	f_orientation_details = open("create_project_data.txt", "w+")
	max_radius = 0.0
	
	for orientation in range(1,number_of_orientations+1):
		output_pdb_name = 'r_'+str(orientation).zfill(3)+'.pdb'
		input_file_name = os.path.normpath(os.getcwd()+'/'+Initial_dimer_pdb)
		with cd("md_input"):
			results = create_random_pdb(separation_distance = separation_distance, move_chain_id = move_chain_id, fix_chain_id = fix_chain_id, input_file_name = input_file_name, output_pdb_name = output_pdb_name, model_number = model_number)
		max_radius = max(max_radius, results["Max_distance"])
		f_orientation_details.write('orientation number\t'+str(orientation)+'\n\n')
		for key in results:
			f_orientation_details.write(key+'\t\t'+str(results[key])+'\n')
		f_orientation_details.write('-------------------------\n\n')

	max_radius+=Boundary_margin
	f_orientation_details.close()
	
	
	group_names = []
	name = pd['Initial_dimer_pdb'][:-4]
	f_in = open(name+"_recentred"+".in", "r")
	for line in f_in:
		if line.strip().split()[:1] == ['group']:
			group_names.append(line)
	f_in.close()
	os.remove(name+"_recentred"+".in")
	
	copy(name+"_recentred"+".seq", 'md_input')
	
	
	
	stride = os.path.normpath(Path_to_awsem+'dimer_interface_protocol/stride/stride')
	os.system(stride + ' '+pd['Initial_dimer_pdb']+' > ssweight.stride')
	location = os.path.normpath(Path_to_awsem+"create_project_tools/stride2ssweight.py")
	os.system(Python2_command+' '+location+' > md_input/ssweight')
	
	
	
	with cd("md_input"):	
		for orientation in range(1,number_of_orientations+1):
			file_name_start = 'r_'+str(orientation).zfill(3)
			location = os.path.normpath(Path_to_awsem+"create_project_tools/PDBToCoordinates.py")
			os.system(Python2_command+" "+location+' '+file_name_start+" "+file_name_start+".coord")
			location = os.path.normpath(Path_to_awsem+"create_project_tools/CoordinatesToWorkLammpsDataFile.py")
			os.system(Python2_command+" "+location+" "+file_name_start+".coord "+file_name_start+".data -b")
			os.remove(file_name_start+".in")
			os.remove(file_name_start+".seq")
			f = open(file_name_start+".pbs", "w+")
			f.write("#!/bin/bash\n")
			f.write("#PBS -S /bin/bash\n")
			f.write("#PBS -l pmem=512mb\n")
			f.write("#PBS -l nodes=1:ppn=1\n")
			f.write("#PBS -l walltime="+str(hours).zfill(2)+':'+str(minutes).zfill(2)+':00\n')
			f.write("#PBS -N "+Jobname+str(orientation).zfill(3)+'\n')
			f.write("cd $PBS_O_WORKDIR\n")
			f.write(Path_to_lmp_serial+" < r_"+str(orientation).zfill(3)+".in\n")
			f.close()


		f_submit_all_pbs = open("submitall.sh", "w+")
		for orientation in range(1,number_of_orientations+1):
			f_submit_all_pbs.write('qsub '+'r_'+str(orientation).zfill(3)+'.pbs >> submited.txt\n')
		f_submit_all_pbs.close()
		
		src = os.path.normpath(Path_to_awsem + '/dimer_interface_protocol/files')
		src_files = os.listdir(src)
		dest = os.getcwd()
		for file_name in src_files:
			full_file_name = os.path.join(src, file_name)
			if (os.path.isfile(full_file_name)):
				copy(full_file_name, dest)
		
		
		
		
		
		first_chain_length = dd['first_chain_length']
		second_chain_length = dd['second_chain_length']
		f_fragmem = open("fragsLAMW.mem", "w+")
		f_gro = open("chain1.gro", "r")
		line = next(f_gro)
		line = next(f_gro)
		line = next(f_gro)
		g1 = int(line.strip().split()[0])
		f_gro.close()
		f_gro = open("chain2.gro", "r")
		line = next(f_gro)
		line = next(f_gro)
		line = next(f_gro)
		g2 = int(line.strip().split()[0])
		f_gro.close()
		f_fragmem.write("[Target]\nquery\n\n[Memories]\n")
		f_fragmem.write("chain1.gro %d  %d %d 1\n" %(1, g1, first_chain_length))
		f_fragmem.write("chain2.gro %d %d %d 1" %(1+first_chain_length, g2, second_chain_length))
		f_fragmem.close()
		
		
		
		
		
		
		
		
		
		
		for orientation in range(1,number_of_orientations+1):
			random_integer = np.random.randint(low = 1000, high = 9999999)
			file_name_start = 'r_'+str(orientation).zfill(3)
			
			
			f = open(file_name_start+".in", "w+")

			f.write('# 3d protein simulation\n')
			f.write('\n')
			f.write('units real\n')
			f.write('\n')
			f.write('timestep 5\n')
			f.write('\n')
			f.write('dimension\t3\n')
			f.write('\n')
			f.write('boundary f f f\n')
			f.write('\n')
			f.write('log\t'+file_name_start+'.log\t  \n')
			f.write('neighbor\t10 bin\n')
			f.write('neigh_modify\tdelay 5\n')
			f.write('\n')
			f.write('atom_modify sort 0 0.0\n')
			f.write('\n')
			f.write('special_bonds fene\n')
			f.write('\n')
			f.write('region\tr1 sphere 0.0 0.0 0.0 {0:.2f} side in \n'.format(max_radius))
			f.write('\n')
			f.write('atom_style\tawsemmd\n')
			f.write('\n')
			f.write('\n')
			f.write('bond_style harmonic\n')
			f.write('\n')
			f.write('pair_style vexcluded 2 3.5 3.5\n')
			f.write('\n')
			f.write('read_data '+file_name_start+'.data\n')
			f.write('\n')
			f.write('pair_coeff * * 0.0\n')
			f.write('pair_coeff 1 1 20.0 3.5 4.5\n')
			f.write('pair_coeff 1 4 20.0 3.5 4.5\n')
			f.write('pair_coeff 4 4 20.0 3.5 4.5\n')
			f.write('pair_coeff 3 3 20.0 3.5 3.5\n')
			f.write('\n')
			f.write('\n')
			f.write('velocity\tall create 300.0 '+str(random_integer)+'\n')
			f.write('\n')
			for line in group_names:
				f.write(line)
				f.write('\n')
			n = dd['first_chain_max_id']
			Dump_time = pd['Dump_time']
			Restart_time = pd['Restart_time']
			f.write('group\t\tchain_1 id <= %d\n' % (n))
			f.write('group\t\tchain_2 id >= %d\n' % (n+1))
			f.write('\n')
			f.write('fix\t\t  1 all nvt temp 300.0 300.0 10.0\n')
			f.write('fix\t\t  2 alpha_carbons backbone beta_atoms oxygens fix_backbone_coeff.data '+name+"_recentred"+'.seq\n')
			f.write('fix\t\t  3 all wall/region r1 harmonic 10.0 1.0 5.0\n')
			if dd['first_chain_is_bigger']:
				f.write('fix\t\t  4 chain_1 recenter 0.0 0.0 0.0 \n')
			else:
				f.write('fix\t\t  4 chain_2 recenter 0.0 0.0 0.0 \n')
			f.write('\n')
			f.write('\n')
			f.write('\n')
			f.write('\n')
			f.write('thermo_style\tcustom step etotal pe ke temp evdwl enthalpy eangle epair emol\n')
			f.write('thermo\t\t5000\n')
			f.write('dump\t\t1 all atom '+str(Dump_time)+' '+file_name_start+'.lammpstrj\n')
			f.write('\n')
			f.write('dump_modify\t1 sort id\n')
			f.write('\n')
			f.write('restart\t\t%d '% (5000)+file_name_start+'.restarttemp1 '+file_name_start+'.restarttemp2\n' )
			f.write('restart\t\t%d '% (Restart_time)+file_name_start+'.restart\n' )
			f.write('\n')
			f.write('variable E_bond  equal emol\n')
			f.write('variable E_chain equal f_2[1]\n')
			f.write('variable E_excl  equal epair\n')
			f.write('variable E_chi   equal f_2[3]\n')
			f.write('variable E_rama  equal f_2[4]\n')
			f.write('variable E_dssp  equal f_2[6]\n')
			f.write('variable E_pap   equal f_2[7]\n')
			f.write('variable E_water equal f_2[8]\n')
			f.write('variable E_helix equal f_2[10]\n')
			f.write('variable E_fmem  equal f_2[12]\n')
			f.write('variable E_P     equal v_E_chain+v_E_chi+v_E_rama+v_E_water+v_E_helix+v_E_fmem+v_E_excl+v_E_bond+v_E_dssp+v_E_pap\n')
			f.write('variable E_K     equal ke\n')
			f.write('variable E_total equal v_E_P+v_E_K\n')
			f.write('variable e_total equal etotal\n')
			f.write('variable Step equal step\n')
			f.write('variable p_e equal pe\n')
			f.write('fix energy all print 5000 "${Step} ${e_total} ${p_e} ${E_K} ${E_chain} ${E_bond} ${E_chi} ${E_rama} ${E_excl} ${E_dssp} ${E_pap} ${E_water} ${E_helix} ${E_fmem} ${E_P} ${E_total} " file '+file_name_start+'_energy.log screen no\n')
			f.write('\n')
			f.write('\n')
			f.write('\n')
			f.write('\n')
			f.write('reset_timestep\t0\n')
			f.write('run\t\t'+str(steps)+'\n')
			f.close()
				
		
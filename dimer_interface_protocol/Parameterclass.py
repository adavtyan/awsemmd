import sys

class Parameterclass(object):
	def __init__(self):
		self.parameterdic = {}
		self._values = []
		self.deriveddic = {}
	
	def read_parameters(self, filename):
		with open(filename) as f:
			content = f.readlines()
			self._values = [x.strip() for x in content[1::2]]
		if len(self._values)!=22:
			#raise InputError(expression = 'filename', message = 'Invalid parameters file')
			print("Invalid parameters file")
			sys.exit(0)
		else:
			self.set()
	
	def set(self):
		parameternames = [
			"Jobname",
			"Initial_dimer_pdb",
			"Timesteps",
			"Simulation_hours",
			"Simulation_minutes",
			"COM_separation",
			"Boundary_margin",
			"Number_of_orientations",
			"Plot_x_range",
			"Plot_energy",
			"Replot_only",
			"Delete_PDB_trajs",
			"Contact_map_max_distance",
			"Contact_map_min_distance",
			"Energy_graph_max",
			"Energy_graph_min",
			"Skip_initial_frames",
			"Path_to_awsem",
			"Path_to_lmp_serial",
			"Python2_command",
			"Dump_time",
			"Restart_time"
		]	
		self.parameterdic = dict(zip(parameternames, self._values))
		
		d = self.parameterdic
		d['Timesteps'] = int(d['Timesteps'])
		d['Simulation_hours'] = int(d['Simulation_hours'])
		d['Simulation_minutes'] = int(d['Simulation_minutes'])
		d['COM_separation'] = float(d['COM_separation'])
		d['Boundary_margin'] = float(d['Boundary_margin'])
		d['Number_of_orientations'] = int(d['Number_of_orientations'])
		d['Plot_x_range'] = int(d['Plot_x_range'])
		d['Contact_map_max_distance'] = float(d['Contact_map_max_distance'])
		d['Contact_map_min_distance'] = float(d['Contact_map_min_distance'])
		d['Energy_graph_max'] = float(d['Energy_graph_max'])
		d['Energy_graph_min'] = float(d['Energy_graph_min'])
		d['Skip_initial_frames'] = int(d['Skip_initial_frames'])
		d['Dump_time'] = int(d['Dump_time'])
		d['Restart_time'] = int(d['Restart_time'])
		if d['Plot_energy'].lower() == 'yes':
			d['Plot_energy'] = True
		elif d['Plot_energy'].lower() == 'no':
			d['Plot_energy'] = False
		else:
			print("Error in parameters file")
			sys.exit(1)
		if d['Replot_only'].lower() == 'yes':
			d['Replot_only'] = True
		elif d['Replot_only'].lower() == 'no':
			d['Replot_only'] = False
		else:
			print("Error in parameters file")
			sys.exit(1)
		if d['Delete_PDB_trajs'].lower() == 'yes':
			d['Delete_PDB_trajs'] = True
		elif d['Delete_PDB_trajs'].lower() == 'no':
			d['Delete_PDB_trajs'] = False
		else:
			print("Error in parameters file")
			sys.exit(1)
	def save_derived(self):
		f_out = open("_config.dat", "w+")
		for k, v in self.deriveddic.items():
			f_out.write(k+'\n'+str(v)+'\n')
		f_out.close()
	def read_derived(self):
		with open('_config.dat', 'r') as f:
			content = f.readlines()
			v = [x.strip() for x in content[1::2]]
			k = [x.strip() for x in content[0::2]]
			self.deriveddic = dict(zip(k, v))
			d = self.deriveddic
		
			d['first_chain_max_id'] = int(d['first_chain_max_id'])
			d['second_chain_length'] = int(d['second_chain_length'])
			d['first_chain_length'] = int(d['first_chain_length'])
			d['first_chain_is_bigger'] = (d['first_chain_is_bigger']=='True')
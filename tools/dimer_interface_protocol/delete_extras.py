import os

def delete_extras(parametersobject):
	bool = parametersobject.parameterdic['Delete_PDB_trajs']
	number_of_orientations = parametersobject.parameterdic['Number_of_orientations']
	
	if bool:
		for i in range(1, 1+number_of_orientations):
			filenamestart = 'pdb_trajectories/t_'+str(i).zfill(3)
			os.remove(filenamestart+'.pdb')
			os.remove(filenamestart+'.psf')
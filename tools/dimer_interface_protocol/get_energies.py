def get_energies(parametersobject):
	pd = parametersobject.parameterdic
	number_of_orientations = pd['Number_of_orientations']
	for i in range(1, 1+number_of_orientations):
		f_in = open('md_output/r_'+str(i).zfill(3)+'_energy.log', "r")
		f_out = open('analysis/e_'+str(i).zfill(3)+'.txt', "w+")
		next(f_in)
		for j, line in enumerate(f_in):
			data = line.strip().split()
			if data[0].isdigit():
				f_out.write(data[0]+'\t'+data[1]+'\n')
				
			else:
				break
		f_in.close()
		f_out.close()
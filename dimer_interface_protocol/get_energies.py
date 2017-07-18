def get_energies(parametersobject):
	pd = parametersobject.parameterdic
	number_of_orientations = pd['Number_of_orientations']
	for i in range(1, 1+number_of_orientations):
		f_in = open('r_'+str(i).zfill(3)+'_energy.log', "r")
		f_out = open('e_'+str(i).zfill(3)+'.txt', "w+")
		next(f_in)
		for j, line in enumerate(f_in):
			data = line.strip().split()
			if data[0].isdigit():
				f_out.write(line[0]+'\t'+line[1]+'\n')
				
			else:
				break
		f_in.close()
		f_out.close()
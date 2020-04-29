import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import matplotlib.gridspec as gridspec
import os, sys




def graph_angles(parametersobject):
	pd = parametersobject.parameterdic
	dd = parametersobject.deriveddic
	x_range = pd['Plot_x_range']
	include_energy = pd['Plot_energy']
	number_of_orientations = pd['Number_of_orientations']
	Dump_time = pd['Dump_time']
	e_max = pd['Energy_graph_max']
	e_min = pd['Energy_graph_min']
	alpha = 0.01
	data = []
	energy=[]
	
	
	if include_energy:
		gs = gridspec.GridSpec(3, 3)
	else:
		gs = gridspec.GridSpec(2, 3)
		
	gs.update(hspace=0.3)
	y_tick = np.arange(-2*np.pi, 2*np.pi+0.1, 0.5*np.pi)
	y_label = [r"$" + format(r/np.pi, ".2g")+ r"\pi$" for r in y_tick]
	y_tick2 = np.arange(0, np.pi+0.1, 0.25*np.pi)
	y_label2 = [r"$" + format(r/np.pi, ".2g")+ r"\pi$" for r in y_tick2]
	y_tick3 = np.arange(-0.5*np.pi, 0.5*np.pi+0.1, 0.25*np.pi)
	y_label3 = [r"$" + format(r/np.pi, ".2g")+ r"\pi$" for r in y_tick3]

	for i in range(1,1+number_of_orientations):
		sys.stdout.write('Part 1 of 2. Progress: %d out of %d\r' % (i, number_of_orientations))
		sys.stdout.flush()
		data.append(np.loadtxt('analysis/angles_'+str(i).zfill(3)+'.txt', unpack=True))
		if include_energy:
			energy.append(np.loadtxt('analysis/e_'+str(i).zfill(3)+'.txt', unpack=True))
		fig=plt.figure(i, figsize=(15,8))

		subplot = plt.subplot(gs[0,0])
		plt.plot(data[i-1][0], data[i-1][1], '.')
		plt.axis([0,x_range,24,46])#distance
		subplot.set_title("Distance")
		subplot = plt.subplot(gs[0,1])
		plt.plot(data[i-1][0], (data[i-1][2]), '.')
		plt.axis([0,x_range,0,3.2])#theta
		subplot.set_yticks(y_tick2)
		subplot.set_yticklabels(y_label2)
		subplot.set_title(r"$\theta$")
		subplot = plt.subplot(gs[0,2])
		plt.plot(data[i-1][0], (data[i-1][3]), '.')
		plt.plot(data[i-1][0], (data[i-1][3]+2*np.pi), '.')
		plt.plot(data[i-1][0], (data[i-1][3]-2*np.pi), '.')
		plt.axis([0,x_range,-6,6])#phi
		subplot.set_yticks(y_tick)
		subplot.set_yticklabels(y_label)
		subplot.set_title(r"$\phi$")
		subplot = plt.subplot(gs[1,0])
		plt.plot(data[i-1][0], (data[i-1][4]), '.')
		plt.plot(data[i-1][0], (data[i-1][4]+2*np.pi), '.')
		plt.plot(data[i-1][0], (data[i-1][4]-2*np.pi), '.')
		plt.axis([0,x_range,-6,6])
		subplot.set_yticks(y_tick)
		subplot.set_yticklabels(y_label)
		subplot.set_title(r"$\theta_x$")
		subplot = plt.subplot(gs[1,1])
		plt.plot(data[i-1][0], (data[i-1][5]), '.')
		plt.axis([0,x_range,-1.6,1.6])
		subplot.set_yticks(y_tick3)
		subplot.set_yticklabels(y_label3)
		subplot.set_title(r"$\theta_y$")
		subplot = plt.subplot(gs[1,2])
		plt.plot(data[i-1][0], (data[i-1][6]), '.')
		plt.plot(data[i-1][0], (data[i-1][6]+2*np.pi), '.')
		plt.plot(data[i-1][0], (data[i-1][6]-2*np.pi), '.')
		plt.axis([0,x_range,-6,6])
		subplot.set_yticks(y_tick)
		subplot.set_yticklabels(y_label)
		subplot.set_title(r"$\theta_z$")
		if include_energy:
			subplot = plt.subplot(gs[2,:])
			plt.plot(energy[i-1][0]/Dump_time, energy[i-1][1], '.')
			plt.axis([0,x_range,e_min,e_max])
		plt.savefig('results_individual/dot_'+str(i).zfill(3)+'.png', bbox_inches = 'tight')
		plt.close(fig)
		#plt.show()
		

	for i in range(1,1+number_of_orientations):
		sys.stdout.write('Part 2 of 2: Consolidating...\r')
		sys.stdout.flush()
		#data.append(np.loadtxt('analysis/angles_'+str(i).zfill(3)+'.txt', unpack=True))
		#energy.append(np.loadtxt('analysis/e_'+str(i).zfill(3)+'.txt', unpack=True))
		fig=plt.figure(100, figsize=(15,8))
		subplot = plt.subplot(gs[0,0])
		plt.plot(data[i-1][0], data[i-1][1], 'b.', lw=1, alpha = alpha)
		plt.axis([0,x_range,24,46])
		subplot.set_title("Distance")
		subplot = plt.subplot(gs[0,1])
		plt.plot(data[i-1][0], (data[i-1][2]), 'b.', lw=1, alpha = alpha)
		plt.axis([0,x_range,0,3.2])
		subplot.set_yticks(y_tick2)
		subplot.set_yticklabels(y_label2)
		subplot.set_title(r"$\theta$")
		subplot = plt.subplot(gs[0,2])
		plt.plot(data[i-1][0], (data[i-1][3]), 'b.', lw=1, alpha = alpha)
		plt.plot(data[i-1][0], (data[i-1][3]+2*np.pi), 'b.', lw=1, alpha = alpha)
		plt.plot(data[i-1][0], (data[i-1][3]-2*np.pi), 'b.', lw=1, alpha = alpha)
		plt.axis([0,x_range,-6,6])
		subplot.set_yticks(y_tick)
		subplot.set_yticklabels(y_label)
		subplot.set_title(r"$\phi$")
		subplot = plt.subplot(gs[1,0])
		plt.plot(data[i-1][0], (data[i-1][4]), 'b.', lw=1, alpha = alpha)
		plt.plot(data[i-1][0], (data[i-1][4]+2*np.pi), 'b.', lw=1, alpha = alpha)
		plt.plot(data[i-1][0], (data[i-1][4]-2*np.pi), 'b.', lw=1, alpha = alpha)
		plt.axis([0,x_range,-6,6])
		subplot.set_yticks(y_tick)
		subplot.set_yticklabels(y_label)
		subplot.set_title(r"$\theta_x$")
		subplot = plt.subplot(gs[1,1])
		plt.plot(data[i-1][0], (data[i-1][5]), 'b.', lw=1, alpha = alpha)
		plt.axis([0,x_range,-1.6,1.6])
		subplot.set_yticks(y_tick3)
		subplot.set_yticklabels(y_label3)
		subplot.set_title(r"$\theta_y$")
		subplot = plt.subplot(gs[1,2])
		plt.plot(data[i-1][0], (data[i-1][6]), 'b.', lw=1, alpha = alpha)
		plt.plot(data[i-1][0], (data[i-1][6]+2*np.pi), 'b.', lw=1, alpha = alpha)
		plt.plot(data[i-1][0], (data[i-1][6]-2*np.pi), 'b.', lw=1, alpha = alpha)
		plt.axis([0,x_range,-6,6])
		subplot.set_yticks(y_tick)
		subplot.set_yticklabels(y_label)
		subplot.set_title(r"$\theta_z$")
		if include_energy:
			subplot = plt.subplot(gs[2,:])
			plt.plot(energy[i-1][0]/Dump_time, energy[i-1][1], 'b.', lw=1, alpha=alpha)
			plt.axis([0,x_range,e_min,e_max])
	plt.savefig('results_main/dot_all.png', bbox_inches = 'tight')
	plt.close(fig)
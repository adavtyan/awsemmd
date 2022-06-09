import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from sklearn.cluster import MeanShift
import matplotlib.gridspec as gridspec
from Bio.PDB import Vector
import os, sys

def cluster(parametersobject):
	number_of_orientations = parametersobject.parameterdic['Number_of_orientations']
	skip = parametersobject.parameterdic['Skip_initial_frames']
	data = []
	
	f_framelist = open('analysis/frames_read.txt', 'w+')
	framelist = []
	
	for i in range(1, 1+number_of_orientations):
		with open('analysis/coord_matrix_'+str(i).zfill(3)+'.txt', 'r') as f:
			l = f.readlines()
			data.extend([[float(p) for p in line.strip().split()[1:10]] for line in l[skip:]])
			f_framelist.write(str(len(l))+'\n')
			framelist.append(len(l)-skip)
	f_framelist.close()
	
	
	ms = MeanShift(n_jobs = -2, cluster_all = True)
	ms.fit(data)
	labels = ms.labels_
	names, numbers = np.unique(labels, return_counts = True)
	cluster_centers = ms.cluster_centers_
	fout = open('analysis/clusters.txt', 'w+')
	fout.write('label\tcount\ttheta\tphi\ttheta_x\ttheta_y\ttheta_z\tx\ty\tz\tR[0][0]\tR[0][1]\tR[0][2]\tR[1][0]\tR[1][1]\tR[1][2]\tR[2][0]\tR[2][1]\tR[2][2]\n')
	phi = 0
	for i, line in enumerate(cluster_centers):
		R = [_[:] for _ in [[]]*9]
		x, y, z = line[0], line[1], line[2]
		R[0] = line[3:6]
		R[1] = line[6:9]
		R[2] = np.cross(R[0], R[1])
		V = Vector(x, y, z)
		if V.norm() > 1e-6:
			theta = V.angle(Vector(0,0,1))
			norm = np.sqrt(x*x + y*y)
			if norm > 1e-6:
				phi = np.arctan2(y,x)
				#otherwise phi isn't updated and the previous value is copied. Keeps it from jumping near the poles.
		else:
			theta = 0.0
		
		theta_x = np.arctan2(R[2][1], R[2][2])
		theta_y = np.arctan2(-R[2][0], np.sqrt(R[2][1]*R[2][1]+R[2][2]*R[2][2]))
		theta_z = np.arctan2(R[1][0], R[0][0])
		
		fout.write(str(names[i])+'\t'+str(numbers[i])+'\t')
		fout.write(str(theta)+'\t'+str(phi)+'\t')
		fout.write(str(theta_x)+'\t'+str(theta_y)+'\t'+str(theta_z)+'\t')
		
		
		for value in line:
			fout.write(str(value)+'\t')
		fout.write(str(R[2][0])+'\t'+str(R[2][1])+'\t'+str(R[2][2])+'\t')
		fout.write('\n')
	fout.close()
	n_clusters_ = len(np.unique(labels))
	print("Number of estimated clusters:", n_clusters_)
	
	
	classification = open('analysis/frame_cluster_types.txt', 'w+')
	frame_counter = 0
	framelist_counter = 0
	for label in labels:
		if frame_counter<framelist[framelist_counter]:
			classification.write(str(label)+'\t')
			frame_counter+=1
		else:
			framelist_counter+=1
			classification.write('\n'+str(label)+'\t')
			frame_counter=1
	
	'''
	e = enumerate(labels)
	print(framelist)
	for p in framelist:
		count = -1
		while count<p:
			count, v = next(e)
			print(count)
			print(v)
			classification.write(str(v)+'\t')
		classification.write('\n')
	'''
		
def cluster_dotplot(parametersobject):
	pd = parametersobject.parameterdic
	dd = parametersobject.deriveddic
	x_range = pd['Plot_x_range']
	include_energy = pd['Plot_energy']
	#number_of_orientations = pd['Number_of_orientations']
	number_of_orientations = 5
	Dump_time = pd['Dump_time']
	e_max = pd['Energy_graph_max']
	e_min = pd['Energy_graph_min']

	skip = pd['Skip_initial_frames']
	alpha = 0.01
	data = []
	energy=[]
	
	labels = []
	with open('analysis/frame_cluster_types.txt', 'r') as f_guide:
		l = f_guide.readlines()
		for line in l:
			labels.append([int(p) for p in line.strip().split()])
		
	colors = 10*['r.','g.','k.','m.','c.','y.','b.']
	
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
		sys.stdout.write('Progress: %d out of %d\r' % (i, number_of_orientations))
		sys.stdout.flush()
		data=(np.loadtxt('analysis/angles_'+str(i).zfill(3)+'.txt', unpack=True))
		if include_energy:
			energy=(np.loadtxt('analysis/e_'+str(i).zfill(3)+'.txt', unpack=True))
		fig=plt.figure(100, figsize=(15,8))
		
		subplot = plt.subplot(gs[0,0])
		plt.plot(data[0][:skip], data[1][:skip], 'b.', lw=1, alpha = alpha)
		for j,p in enumerate(labels[i-1]):
			plt.plot(data[0][skip+j], data[1][skip+j], colors[labels[i-1][j]], lw=1, alpha = alpha)
		plt.axis([0,x_range,24,46])
		subplot.set_title("Distance")
		
		subplot = plt.subplot(gs[0,1])
		plt.plot(data[0][:skip], (data[2][:skip]), 'b.', lw=1, alpha = alpha)
		for j,p in enumerate(labels[i-1]):
			plt.plot(data[0][skip+j], (data[2][skip+j]), colors[labels[i-1][j]], lw=1, alpha = alpha)
		plt.axis([0,x_range,0,3.2])
		subplot.set_yticks(y_tick2)
		subplot.set_yticklabels(y_label2)
		subplot.set_title(r"$\theta$")
		
		subplot = plt.subplot(gs[0,2])
		plt.plot(data[0][:skip], (data[3][:skip]), 'b.', lw=1, alpha = alpha)
		for j,p in enumerate(labels[i-1]):
			plt.plot(data[0][skip+j], (data[3][skip+j]), colors[labels[i-1][j]], lw=1, alpha = alpha)
		plt.plot(data[0][:skip], (data[3][:skip]+2*np.pi), 'b.', lw=1, alpha = alpha)
		for j,p in enumerate(labels[i-1]):
			plt.plot(data[0][skip+j], (data[3][skip+j]+2*np.pi), colors[labels[i-1][j]], lw=1, alpha = alpha)
		plt.plot(data[0][:skip], (data[3][:skip]-2*np.pi), 'b.', lw=1, alpha = alpha)
		for j,p in enumerate(labels[i-1]):
			plt.plot(data[0][skip+j], (data[3][skip+j]-2*np.pi), colors[labels[i-1][j]], lw=1, alpha = alpha)
		plt.axis([0,x_range,-6,6])
		subplot.set_yticks(y_tick)
		subplot.set_yticklabels(y_label)
		subplot.set_title(r"$\phi$")
		
		subplot = plt.subplot(gs[1,0])
		plt.plot(data[0][:skip], (data[4][:skip]), 'b.', lw=1, alpha = alpha)
		for j,p in enumerate(labels[i-1]):
			plt.plot(data[0][skip+j], (data[4][skip+j]), colors[labels[i-1][j]], lw=1, alpha = alpha)
		plt.plot(data[0][:skip], (data[4][:skip]+2*np.pi), 'b.', lw=1, alpha = alpha)
		for j,p in enumerate(labels[i-1]):
			plt.plot(data[0][skip+j], (data[4][skip+j]+2*np.pi), colors[labels[i-1][j]], lw=1, alpha = alpha)
		plt.plot(data[0][:skip], (data[4][:skip]-2*np.pi), 'b.', lw=1, alpha = alpha)
		for j,p in enumerate(labels[i-1]):
			plt.plot(data[0][skip+j], (data[4][skip+j]-2*np.pi), colors[labels[i-1][j]], lw=1, alpha = alpha)
		plt.axis([0,x_range,-6,6])
		subplot.set_yticks(y_tick)
		subplot.set_yticklabels(y_label)
		subplot.set_title(r"$\theta_x$")
		
		subplot = plt.subplot(gs[1,1])
		plt.plot(data[0][:skip], (data[5][:skip]), 'b.', lw=1, alpha = alpha)
		for j,p in enumerate(labels[i-1]):
			plt.plot(data[0][skip+j], (data[5][skip+j]), colors[labels[i-1][j]], lw=1, alpha = alpha)
		plt.axis([0,x_range,-1.6,1.6])
		subplot.set_yticks(y_tick3)
		subplot.set_yticklabels(y_label3)
		subplot.set_title(r"$\theta_y$")
		
		subplot = plt.subplot(gs[1,2])
		plt.plot(data[0][:skip], (data[6][:skip]), 'b.', lw=1, alpha = alpha)
		for j,p in enumerate(labels[i-1]):
			plt.plot(data[0][skip+j], (data[6][skip+j]), colors[labels[i-1][j]], lw=1, alpha = alpha)
		plt.plot(data[0][:skip], (data[6][:skip]+2*np.pi), 'b.', lw=1, alpha = alpha)
		for j,p in enumerate(labels[i-1]):
			plt.plot(data[0][skip+j], (data[6][skip+j]+2*np.pi), colors[labels[i-1][j]], lw=1, alpha = alpha)
		plt.plot(data[0][:skip], (data[6][:skip]-2*np.pi), 'b.', lw=1, alpha = alpha)
		for j,p in enumerate(labels[i-1]):
			plt.plot(data[0][skip+j], (data[6][skip+j]-2*np.pi), colors[labels[i-1][j]], lw=1, alpha = alpha)
		plt.axis([0,x_range,-6,6])
		subplot.set_yticks(y_tick)
		subplot.set_yticklabels(y_label)
		subplot.set_title(r"$\theta_z$")
		
		if include_energy:
			subplot = plt.subplot(gs[2,:])
			plt.plot(energy[0][:skip-1]/Dump_time, energy[1][:skip-1], 'b.', lw=1, alpha=alpha)
			for j,p in enumerate(labels[i-1]):
				plt.plot(energy[0][skip-1+j]/Dump_time, energy[1][skip-1+j], colors[labels[i-1][j]], lw=1, alpha=alpha)
			plt.axis([0,x_range,e_min,e_max])
	plt.savefig('results_main/dot_all_clustered.png', bbox_inches = 'tight')
	plt.close(fig)

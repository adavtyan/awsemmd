# input is the TMHMM result from http://www.cbs.dtu.dk/services/TMHMM/
# output is a zim file where each line is either a 1,2,3 where 1 and 3 represent
# outermembrane/hydrophilic regions and 2 represents the transmembrane region
#  example:
#  Sequence TMHMM2.0 outside     1     9
#  Sequence TMHMM2.0 TMhelix    10    32
#  Sequence TMHMM2.0 inside    33    52
#  Sequence TMHMM2.0 TMhelix    53    75
#  Sequence TMHMM2.0 outside    76    79


import sys

if len(sys.argv)!=2:
	print "\nPrepzimfile.py zimfile_input\n"
	exit()

zimfile_input = sys.argv[1]

def zim_writer(var):
	for i in range(0,num_res):
		out.write(var)
		out.write('\n')


file = open(zimfile_input)
out = open('zim','w')
zim_array = []

for line in file:
	line = line.strip()
	line = line.split()
	zim_array.append(line)

for line in zim_array:
	num_res = int(line[4])-int(line[3])+1
	if line[2]=='outside':
		zim_writer('1')
	elif line[2]=='TMhelix':
		zim_writer('2')
	elif line[2]=='inside':
		zim_writer('3')
	else:
		print "\nINVALID REGION ASSIGNMENT TO RESIDUES\n"
		exit()

file.close()
out.close()

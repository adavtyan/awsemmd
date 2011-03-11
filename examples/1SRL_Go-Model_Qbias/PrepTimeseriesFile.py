import sys

if len(sys.argv)!=3 and len(sys.argv)!=4:
	print "PrepTimeseriesFile.py input_file output_file [start_from]"
	exit()

input_filename = sys.argv[1]
output_filename = sys.argv[2]

istart = 0
if len(sys.argv)==4:
	istart = int(sys.argv[3])

file = open(input_filename, 'r')
q = file.read()
file.close()

q = q.split()

out = open(output_filename, 'w')
for i in range(istart, len(q)):
	out.write(str(i*100))
	out.write('\t')
	out.write(q[i])
	out.write('\n')
out.close()

import sys
from pylab import *

cols=[]
files=[]
legends=[]
desc=""
fpng = ""

# Read command parameters

if len(sys.argv)<4:
	print ("\n>", sys.argv[0], "file1 [file2 file3 ...] column1 [column2 column3 ...] legend1 [legend2 legend3 ...] [desc_text] [image_file_name]\n")
	sys.exit()

# Read files
nplots=0
for iarg in sys.argv[1:]:
	if iarg.isdigit(): 
		break
	files.append(iarg)
	nplots += 1

if len(sys.argv)<3*nplots+1:
	print ("\n>", sys.argv[0], "file1 [file2 file3 ...] column1 [column2 column3 ...] legend1 [legend2 legend3 ...] [desc_text] [image_file_name]\n")
	sys.exit()
	

# Read coulumns to be ploted
for i in range(nplots):
	iarg = sys.argv[nplots + i + 1]
	if not iarg.isdigit():
		print ("Column number was specified in a wrong format!\n")
		sys.exit()
	cols.append(int(iarg))

# Read legends
for i in range(nplots):
	iarg = sys.argv[2*nplots + i + 1]
	legends.append(iarg)

if len(sys.argv)>3*nplots+1:
	desc = sys.argv[3*nplots + 1]

if len(sys.argv)>3*nplots+2:
	fpng = sys.argv[3*nplots + 2]


figure()
ymax = 0.0
for i in range(nplots):
	ifile=files[i]

	y=[]
	x=[]

	file=open(ifile, 'r')

	for l in file:
		if l[0]=='#': continue
		l = l.strip().split()
	
		x.append(float(l[0]))
		icol=cols[i]
		y.append(float(l[icol]))
		if y[-1]>ymax: ymax=y[-1]

	file.close()

	plot(x, y, label=legends[i], linewidth=1.0)

xlabel('r')
ylabel('RDF')
title('Plot of Radial Distribution Function')
if desc!="": text((x[-1]-x[0])/2.0, ymax*0.98, desc, va='top', ha='center',  fontsize=18)
grid(True)
legend()

if fpng!="":
	print ("Saving plot to", fpng, "...")
	savefig(fpng)

show()


#!/usr/bin/python

# ----------------------------------------------------------------------
# Copyright (2010) Aram Davtyan and Garegin Papoian

# Papoian's Group, University of Maryland at Collage Park
# http://papoian.chem.umd.edu/

# Last Update: 03/04/2011
# ----------------------------------------------------------------------

import sys

if len(sys.argv)!=3:
	print "\n" + sys.argv[0] + " inpute_file output_file\n"
	exit()

input_file = sys.argv[1]
output_file = sys.argv[2]

inp = open(input_file, 'r')
st = inp.read().strip()
inp.close()

out =  open(output_file, 'w')
for s in st:
	a = 0.0
	b = 0.0
	if s=='H':
		a = 1.0
	elif s=='E':
		b = 1.0
	out.write(str(round(a,1)))
	out.write(' ')
	out.write(str(round(b,1)))
	out.write('\n')
out.close()

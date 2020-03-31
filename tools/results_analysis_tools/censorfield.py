#!python

import sys
import math

file=open(sys.argv[1])

for line in file:
    splitline = line.split()

    vectornorm = math.sqrt(float(splitline[3])**2+float(splitline[4])**2+float(splitline[5])**2)
    if(vectornorm>float(sys.argv[2])):
        print line
    


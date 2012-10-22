#!python

import sys
import math

file=open(sys.argv[1])
distancerange=int(sys.argv[2])

sfn1 = [0.0]*distancerange
sfn2 = [0.0]*distancerange
sfn3 = [0.0]*distancerange
skewness = [0.0]*distancerange

numpts = [0]*distancerange

# read in file, add structure function data to appropriate distance range
for line in file:
    line = line.strip()
    if line == "": continue
    line = line.split()
    # bin the data by rounding the distance to the nearest number of contacts
    i = int(round(float(line[0])))
    numpts[i] = numpts[i] + 1
    sfn1[i] = sfn1[i] + float(line[1])
    sfn2[i] = sfn2[i] + float(line[2])
    sfn3[i] = sfn3[i] + float(line[3])

# normalize by the number of points and calculate skewness
for i in range(len(numpts)):
    n = float(numpts[i])
    if(n == 0): continue
    sfn1[i]=sfn1[i]/n
    sfn2[i]=sfn2[i]/n
    sfn3[i]=sfn3[i]/n
    skewness[i]=sfn3[i]/math.pow(sfn2[i],1.5)

# print out results
for i in range(len(numpts)):
    # include the number of points in each distance range to indicate
    # how well sampled each region is
    print i,sfn1[i],sfn2[i],sfn3[i],skewness[i],numpts[i]


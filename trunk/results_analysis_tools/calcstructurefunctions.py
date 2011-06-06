#!python

import sys
import math

file=open(sys.argv[1])

vecpts = []

# read in file, store vector field as list of 6D points
for line in file:
    line = line.strip()
    if line == "": continue
    splitline = line.split()
    vecpts.append([splitline[0],splitline[1],splitline[2],splitline[3],splitline[4],splitline[5]])
    
# loop over each pair of points
for i in range(len(vecpts)):
    for j in range(i+1,len(vecpts)):
        # calculate displacement vector
        lvec = [(float(vecpts[i][0])-float(vecpts[j][0])),(float(vecpts[i][1])-float(vecpts[j][1])),(float(vecpts[i][2])-float(vecpts[j][2]))]        
        # calculate distance difference
        distance = math.sqrt(lvec[0]**2+lvec[1]**2+lvec[2]**2)
        print distance,
        # calculate vector difference
        vecdiff = [(float(vecpts[i][3])-float(vecpts[j][3])),(float(vecpts[i][4])-float(vecpts[j][4])),(float(vecpts[i][5])-float(vecpts[j][5]))]
        # print vecdiff,
        # calculate each structure function and print line with distance and each structure function
        sfn1 = abs((vecdiff[0]*lvec[0]+vecdiff[1]*lvec[1]+vecdiff[2]*lvec[2])/distance)
        sfn2 = sfn1**2
        sfn3 = sfn1*sfn2
        print sfn1,sfn2,sfn3    


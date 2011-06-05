#!python

from math import atan2
import sys

def sign(x):
    if(x > 0):
        return 1
    elif(x < 0):
        return -1
    else:
        return 0
        


filename=sys.argv[1]
dimension=int(sys.argv[2])
x=[int(sys.argv[3]),int(sys.argv[4]),int(sys.argv[5])]

# Initialize vector field
flowfield=[]
for i in xrange(dimension):
    flowfield.append([])
    for j in xrange(dimension):
        flowfield[i].append([])
        for k in xrange(dimension):
            flowfield[i][j].append([0,0,0])

# Go and get the first point
firstone = 1

# For all of the rest of the points
for line in file(filename):
    # If this is the first line, get it and go on to the next line
    if(firstone):
        point1 = line.split()
        firstone=0
        continue

    # Get the next line
    point2 = line.split()
    # Loop over all possible points where you might want to add to the vector field
    for xt in range(min(int(point1[x[0]]),int(point2[x[0]])),max(int(point1[x[0]]),int(point2[x[0]]))+1):
        for yt in range(min(int(point1[x[1]]),int(point2[x[1]])),max(int(point1[x[1]]),int(point2[x[1]]))+1):
            for zt in range(min(int(point1[x[2]]),int(point2[x[2]])),max(int(point1[x[2]]),int(point2[x[2]]))+1):
                # For each test point
                tp=[xt,yt,zt]
                # Loop over all components of the vector x, y and z
                for coord in [0,1,2]:
                    # Make sure you are not going to divide by zero
                    if(int(point1[x[coord]]) != int(point2[x[coord]])):
                        # Find where the transition line crosses the plane defined by the test point and the component
                        test1 = float(((float(point2[x[(coord+1)%3]])-float(point1[x[(coord+1)%3]]))/(float(point2[x[coord]])-float(point1[x[coord]])))*(float(tp[coord])-float(point1[x[coord]]))+float(point1[x[(coord+1)%3]]))
                        if((test1 > float(tp[(coord+1)%3])-0.5) and (test1 < float(tp[(coord+1)%3])+0.5)):
                            test2 = float(((float(point2[x[(coord+2)%3]])-float(point1[x[(coord+2)%3]]))/(float(point2[x[coord]])-float(point1[x[coord]])))*(float(tp[coord])-float(point1[x[coord]]))+float(point1[x[(coord+1)%3]]))
                            if((test2 > float(tp[(coord+2)%3])-0.5) and (test2 < float(tp[(coord+2)%3])+0.5)):
                                # If it passes nearby the test point, add the appropriate component to the vector field
                                flowfield[xt][yt][zt][coord]=flowfield[xt][yt][zt][coord]+sign(int(point2[x[coord]])-int(point1[x[coord]]))
                                

        
    point1 = point2

# print the results in a Mathematica friendly way
print "{",
for i in xrange(dimension):
    for j in xrange(dimension):
        for k in xrange(dimension):
            print "{{",i+1,",",j+1,",",k+1,"},{",flowfield[i][j][k][0],",",flowfield[i][j][k][1],",",flowfield[i][j][k][2],"}}",
            if(k != dimension-1): print ",",

print "}"

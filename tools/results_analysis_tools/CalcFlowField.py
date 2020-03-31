#!python

from math import atan2
import sys

# define a function which returns the sign of an argument
# this is used when adding vector components to the field
def sign(x):
    if(x > 0):
        return 1
    elif(x < 0):
        return -1
    else:
        return 0
        
# this is the file that contains all the transition information
# it can, and should, contain data from multiple simulations
# each simulation is separated by a blank line so that nothing
# is added to the vector field for going from the end of one simulation
# to the start of another
filename=sys.argv[1]
# the dimension should be greater than the maximum number of contacts
# formed in any given foldon, it is used to size arrays
dimension=int(sys.argv[2])
# x is a vector that stores the column numbers for the intrafoldon
# or interfoldon number of contacts that you are interested in
x=[int(sys.argv[3]),int(sys.argv[4]),int(sys.argv[5])]
# at the end of calculating the field from all the transitions, each
# vector component is multiplied by a scale factor which should be chosen
# such that the maximum component is of order one, for plotting purposes
scalefactor=sys.argv[6]

# Initialize vector field, it is a 3D array of 3D vectors
flowfield=[]
for i in xrange(dimension):
    flowfield.append([])
    for j in xrange(dimension):
        flowfield[i].append([])
        for k in xrange(dimension):
            flowfield[i][j].append([0,0,0]) # all components are zero initially

# The loop below needs to know if the state it reads in is the first one for
# a given simulation; if it is, it will break out of the loop after reading it in
firstone = 1

# For all of the rest of the points
for line in file(filename):
    # If this is the first line, get it and go on to the next line
    line = line.strip()
    if line == "": # blank lines are used to indicate the start of a new simulation
        firstone = 1 # so the next line it reads in should be a "first one"
        continue
    if(firstone): # as mentioned above, if this is the first state in a simulation,
        point1 = line.split() # read it in and then go immediately to the next state
        firstone=0            # without adding anything to the vector field
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
                    # Make sure you are not going to divide by zero, if so skip because those transitions won't contribute to that component anyway
                    if(int(point1[x[coord]]) != int(point2[x[coord]])):
                        # Find where the transition line crosses the plane defined by the test point and the component
                        # in the first case
                        i = x[coord]       # i.e. x then y then z
                        j = x[(coord+1)%3] #      y then z then x
                        k = x[(coord+2)%3] #      z then x then y
                        p1i = float(point1[i]) # starting point i
                        p1j = float(point1[j]) #                j
                        p1k = float(point1[k]) #                k
                        p2i = float(point2[i]) # ending point   i
                        p2j = float(point2[j]) #                j
                        p2k = float(point2[k]) #                k
                        tpi = float(tp[coord]) # test point     i
                        tpj = float(tp[(coord+1)%3]) #          j
                        tpk = float(tp[(coord+2)%3]) #          k
                        # find the point where the transition line crosses the i=tpi (const) plane
                        test1 = ((p2j-p1j)/(p2i-p1i))*(tpi-p1i)+p1j 
                        if(test1 > tpj-0.5 and test1 < tpj+0.5): # if it passes close by the test point, go to the next test
                           test2 = ((p2k-p1k)/(p2i-p1i))*(tpi-p1i)+p1k # find the other coordinate for the same point
                           if(test2 > tpk-0.5 and test2 < tpk+0.5):
                                # If it passes nearby the test point, add the appropriate component to the vector field
                                flowfield[xt][yt][zt][coord]=flowfield[xt][yt][zt][coord]+sign(int(p2i-p1i))
                                        
    point1 = point2 # make the ending point the new starting point then loop back to get a new ending point

# print the results in a Mathematica friendly way
# print "{",
# for i in xrange(dimension):
#     for j in xrange(dimension):
#         for k in xrange(dimension):
#             print "{{",i,",",j,",",k,"},{",flowfield[i][j][k][0],",",flowfield[i][j][k][1],",",flowfield[i][j][k][2],"}}",
#             if(k != dimension-1): print ",",

# print "}"

# print the results in a gnuplot friendly way
for i in xrange(dimension):
    for j in xrange(dimension):
        for k in xrange(dimension):
            # scale the field
            flowfield[i][j][k][0]=flowfield[i][j][k][0]*float(scalefactor)
            flowfield[i][j][k][1]=flowfield[i][j][k][1]*float(scalefactor)
            flowfield[i][j][k][2]=flowfield[i][j][k][2]*float(scalefactor)
            # print out a line with coordinates and field
            print i,j,k,flowfield[i][j][k][0],flowfield[i][j][k][1],flowfield[i][j][k][2]


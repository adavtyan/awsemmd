# This script was written by Nick Schafer on 4/7/12
# It can be used to determine the fragment coverage for a given
# list of fragments in the format of a ".mem" file, e.g., fragsLAMW.mem.
# The output is a single row of integers that indicates how many
# other residues (summed over all memories) that a particular residue
# interacts with.

import sys

minSep=3
maxSep=12

# Finds all the lines after the [Memories] header and returns them
def getMemoryLines():
    memoryfile = open(sys.argv[1])

    lines = []

    foundMemories = False    

    for line in memoryfile:
        splitline = line.split()
        if len(splitline) == 0:
            continue
        
        if line[0] == '#':
            continue

        if foundMemories == False:
            if splitline[0] == '[Memories]':
                foundMemories = True
                
            continue
            
        lines.append(splitline)
    
    return lines

# Finds the size of the protein by assuming that it is equal to
# the last residue where an associative memory interaction is present
def getNumberOfResidues():
    memoryLines = getMemoryLines()
    
    numberOfResidues = 0
    maxResidueOfLine = 0
    for line in memoryLines:
        maxResidueOfLine = int(line[1]) + int(line[3]) - 1
        if maxResidueOfLine > numberOfResidues:
            numberOfResidues = maxResidueOfLine

    return numberOfResidues

# Main program
if len(sys.argv) != 3:
    print "Syntax:"
    print "\ncalculateFragmentCoverage.py input output\n"
    exit()

numberOfResidues = getNumberOfResidues()

# The array which will contain the number of interacting residues per residue
fragmentCoverage = []

# Zero the array
for i in range(numberOfResidues):
    fragmentCoverage.append(0)

# Loop over all memories
for line in getMemoryLines():
    # Find the starting residue (subtract 1 b/c python arrays start at 0)
    startingResidue = int(line[1])-1
    # Find the fragment length
    fragmentLength = int(line[3])
    # Find the ending residue
    endingResidue = startingResidue + fragmentLength-1

    # Add the appropriate number of interactions to all positions in the fragment
    for position in range(startingResidue,endingResidue+1):
        for otherposition in range(position+minSep,min(position+maxSep,endingResidue+1)):
            fragmentCoverage[position] += 1
            fragmentCoverage[otherposition] += 1

# Write the array out to a file
fragmentCoverageFile = open(sys.argv[2],"w")
for i in range(numberOfResidues):
    fragmentCoverageFile.write(str(fragmentCoverage[i])+" ")



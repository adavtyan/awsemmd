# This script was started by Nick Schafer on 4/14/12.
# The goal is to be able to calculate formulate a coarse-grained
# kinetic model based on foldon, connectivity and rate assumptions
# using umbrella sampled data. A definition of a contact and the
# native structure of the protein is also required.

# The input is a set of snapshots from simulation. For each snapshot,
# the coordinates, energy (including the biasing energy) and global Q 
# value are required.
# {X_i}, V_i, Q_i

# Each protein is, prior to the simulation, divided into foldons,
# which are sets of residues. The only requirements are that this
# set of residues have at least one native contact, either within
# the foldon or with residues in another foldon, and that every
# residue in the protein belong to exactly one foldon.

# Foldon definition input file format:
# 1 2 3 4 5 <-- first foldon
# 6 7 8     <-- second foldon
# 9 10      <-- third foldon
# The above foldons happen to be contiguous in sequence, but this
# isn't required.

# The contact definition to be used will be as simple as possible:
# A threshold distance will be defined and two residues are in contact
# if their C-alpha atoms are within that threshold distance. The same
# threshold is used on the native structure to define the native contacts.

# The number of foldons defines the possible coarse-grained microstates.
# If there are three foldons, the possible coarse grained microstates are:
# 000
# 100
# 010
# 001
# 110
# 011
# 101
# 111
# There are 2^(# of foldons) possible coarse-grained microstates. If the
# fraction of native contacts for a foldon exceeds a predefined threshold,
# then that foldon in the microstate is given a "1"; if not, it is "0". So,
# a conformation where the threshold is exceeded for the first and third
# foldons but not the second will be represented by "101".

# The free energy of each microstate can be calculated by:
# First, summing the Boltzmann factor e^[(V_i-F(Q_i))/kT] over all the
# structures to obtain Z then the probability of any particular structure
# is given by e^[(V_i-F(Q_i))/kT]/Z. Finally, by summing these probabilities
# over all structures that correspond to a certain microstate and taking the
# negative logarithm, F=-kTlog(P), the free energy of all microstates can
# be obtained.

# The coarse-grained microstates are then connected using a Hamming distance
# of 1 condition, i.e., only two states that can be interconverted by one 1->0 
# or 0->1 change are considered directly connected.

# For those pairs of coarse-grained microstates that are found to be present
# in the simulation data and connected by the Hamming distance condition,
# a rate is assigned using the following assumption:
# The rate of going from state i to state j is 
# kij = k0 if F_i > F_j or kij = k0e^[(F_i-F_j)/kT] if F_i < F_j.
# k0 is also predefined.

# A (perhaps sparse) rate matrix is then formulated using kij as the off
# diagonal elements and the negative sum of the off diagonal elements as
# the diagonal elements (to satisfy detailed balance). The eigenvalues
# of this matrix should be real and non-positive, with one 0 eigenvalue
# whose eigenvector corresponds to the equilibrium probability distribution
# of the model.

# The global free energy F(Q) is temperature dependent, and so therefore
# are the microstate probabilities and the rate matrix. By performing this
# calculation for a range of temperatures, and assuming that the relaxation
# process is exponential (which is satisfied so long as the first non-zero
# eigenvalue is well separated from all the rest), one can obtain "thermal
# chevron plots" by plotting (the negative of) the first non-zero eigenvalue
# as a function of temperature.

# Files
# Native structure file (in PDB format)
nativeStructureFile = './1nor.pdb'
# The metadata file, containing links to the dump files and Qw/Potential energy
# files. The format is: dumpfile qw-pot-file
metadataFile = './metadata'
foldonFile = './foldons'

# Parameters
# Native contact threshold, in Angstroms, to be applied to Ca-Ca distances:
nativeContactThreshold = 10 

# Variables and arrays
# Native contact list: {[residue1 residue2],[...],}
nativeContactList = []

class foldon:
    numRes = 0
    residues = []

    def __init__(self,residues):
        self.residues = residues
        self.numRes = len(self.residues)

    def display(self):
        print self.residues

class ustate:
    code = ""
    freeEnergy = 0.0
    
    def __init__(self,code):
        self.code = code

    def display(self):
        print self.code
    
    def isFolded(self,i):
        if(self.code[i-1] == '1'):
            return True
        else:
            return False
    
    def distanceTo(self,otherustate):
        distance = 0
        for i in range(len(self.code)):
            distance += abs(int(self.code[i])-int(otherustate.code[i]))

        return distance

    def isConnectedTo(self,otherustate):
        if(self.distanceTo(otherustate)==1):
            return True
        else:
            return False

    def isDownhillFrom(self,otherustate):
        if(self.freeEnergy < otherustate.freeEnergy):
            return True
        else:
            return False

class Residue:
    pos = 0
    x = 0.0
    y = 0.0
    z = 0.0
    
    def __init__(self,pos,x,y,z):
        self.pos = pos
        self.x = x
        self.y = y
        self.z = z
        
    def display(self):
        print "pos: " + str(pos)
        print "coords: " + str(self.x) + " " + str(self.y) + " " + str(self.z)

class Snapshot:
    numRes = 0
    residues = []
    energy = 0.0
    Q = 0.0
    ustate = 0
    
    def __init__(self, numRes):
        self.numRes = numRes
        self.residues = []
        for i in range(self.numRes):
            residue = Residue(0.0, 0.0, 0.0)
            self.residues.append(residue)
                                           
        self.energy = 0.0
        self.Q = 0.0
        self.ustate = 0

    def display(self):
        print "Ca coordinates:"
        for i in range(self.numRes):
            print "residue " + str(i) + ": " + str(self.residues[i].x) + " " + str(self.residues[i].y) + " " + str(self.residues[i].z)
            
        print "Energy: " + str(self.energy)
        print "Qvalue: " + str(self.Q)
        print "ustate: " + str(self.ustate)

class Trajectory:
    metadataFile = ""
    snapshots = []
    numSnapshots = 0
    numRes = 0
    
    def __init__(self, numRes, numSnapshots):
        self.numSnapshots = numSnapshots
        self.numRes = numRes
        for i in range(self.numSnapshots):
            snapshot = Snapshot(self.numRes)
            self.snapshots.append(snapshot)

    def display(self):
        for i in range(self.numSnapshots):
            self.snapshots[i].display()

# Read in foldon definition file and create foldon objects
foldons = []
f = open(foldonFile,"r")
for line in f:
    residues = []
    line = line.strip()

    if line == "" or line[0] == "#":
        continue

    for i in range(len(line)):
        if line[i] == " ":
            continue
        else:
            residues.append(line[i])

    foldons.append(foldon(residues))

for i in range(len(foldons)):
    foldons[i].display()
    print foldons[i].numRes

# End reading in foldon file


# This script was started by Nick Schafer on 4/14/12.
# The goal is to be able to calculate formulate a coarse-grained
# kinetic model based on foldon, connectivity and rate assumptions
# using umbrella sampled data. A definition of a contact and the
# native structure of the protein is also required.

# The input is a set of snapshots from simulation. For each snapshot,
# the coordinates, biasing energy, and global Q value are required.
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

#########################
# Classes and Functions #
#########################
class Foldon:
    numRes = 0
    residuelist = []
    numNativeContacts = 0

    def __init__(self,residuelist):
        self.residuelist = residuelist
        self.numRes = len(self.residuelist)
        self.numNativeContacts = self.calculateNumNativeContacts(nativeSnapshot,nativeSnapshot)

    def display(self):
        print self.residuelist

    def calculateNumNativeContacts(self, snapshot, native):
        numNativeContacts = 0
        for residue1 in self.residuelist:
            for residue2 in range(len(native.residues)):
                if snapshot.residues[residue1].isInContactWith(snapshot.residues[residue2]):
                    if native.residues[residue1].isInContactWith(native.residues[residue2]):
                        numNativeContacts += 1

        return numNativeContacts

    def isFoldedIn(self, snapshot):
        currentContacts = self.calculateNumNativeContacts(snapshot,nativeSnapshot)
        if float(currentContacts)/float(self.numNativeContacts) > foldonThreshold:
            return True
        else:
            return False

class ustate:
    code = ""
    freeEnergy = 0.0
    
    def __init__(self,code):
        self.code = code

    def display(self):
        print self.code
    
    def isFolded(self,i):
        if(self.code[i] == '1'):
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
        print "pos: " + str(self.pos)
        print "coords: " + str(self.x) + " " + str(self.y) + " " + str(self.z)

    def isInContactWith(self,otherresidue):
        distance = math.sqrt(pow(self.x-otherresidue.x,2)+pow(self.y-otherresidue.y,2)+pow(self.z-otherresidue.z,2))
        if abs(self.pos - otherresidue.pos) >= minSeqSep and distance < nativeContactThreshold:
            return True
        else:
            return False

class Snapshot:
    numRes = 0
    residues = []
    energy = 0.0
    Q = 0.0
    contactmap = []
    ustate = ustate("")
    probability = 0.0
    
    def __init__(self, residues):
        self.numRes = len(residues)
        self.residues = residues                                           
        self.energy = 0.0
        self.Q = 0.0
        self.ustate = 0
        self.contactmap = []
        for i in range(self.numRes):
            self.contactmap.append([])
            for j in range(self.numRes):
                self.contactmap[i].append([])
                self.contactmap[i][j] = 0
        self.calculateContactMap()
        self.ustate = ustate("")
        self.probability = 0.0

    def display(self):
        print "Ca coordinates:"
        for i in range(self.numRes):
            print "residue " + str(i) + ": " + str(self.residues[i].x) + " " + str(self.residues[i].y) + " " + str(self.residues[i].z)
            
        print "Energy: " + str(self.energy)
        print "Qvalue: " + str(self.Q)
        print "Microstate code:"
        self.ustate.display()
        print "Probability: " + str(self.probability)

    def calculateContactMap(self):
        for i in range(self.numRes):
            for j in range(self.numRes):
                if self.residues[i].isInContactWith(self.residues[j]):
                    self.contactmap[i][j] = 1

    def assignUstate(self):
        microstatecode = ""
        for i in range(len(foldons)):
            if foldons[i].isFoldedIn(self):
                microstatecode += '1'
            else:
                microstatecode += '0'

        self.ustate = ustate(microstatecode)
        return microstatecode

    def assignProbability(self, temperature):
        self.probability = boltzmannFactor(self,temperature)/Z

class Trajectory:
    dumpFile = ""
    snapshotDataFile = ""
    snapshots = []
    numRes = 0
    
    def __init__(self, dumpFile, snapshotDataFile):
        self.dumpFile = dumpFile
        self.snapshotDataFile = snapshotDataFile
        self.snapshots = readDumpFile(dumpFile)
        ssdata = open(snapshotDataFile, 'r')
        snapshotindex = 0
        for dataline in ssdata:
            dataline = dataline.split()
            self.snapshots[snapshotindex].Q = dataline[0]
            self.snapshots[snapshotindex].energy = dataline[1]
            snapshotindex += 1

    def display(self):
        for i in range(len(self.snapshots)):
            self.snapshots[i].display()
    
    def assignAllUstates(self):
        for i in range(len(self.snapshots)):
            self.snapshots[i].assignUstate()

    def assignAllProbabilities(self,temperature):
        for i in range(len(self.snapshots)):
            self.snapshots[i].assignProbability(temperature)

def readFoldonFile(foldonFile):
    foldons = []
    f = open(foldonFile,"r")
    for line in f:
        residues = []
        line = line.split()

        if line == "" or line[0] == "#":
            continue

        for i in range(len(line)):
            if line[i] == " ":
                continue
            else:
                residues.append(int(line[i]))

        foldons.append(Foldon(residues))

    return foldons

def readDumpFile(dumpFile):
    foundAtoms = False
    residues = []
    residuePosition = 0
    snapshots = []
    firstTimestep = True
    
    foundBoxBounds = False
    boundsIndex = 0
    bounds = []

    for i in range(3):
        bounds.append([0.0, 0.0])

    f = open(dumpFile,"r")
    for line in f:
        line = line.split()

        if len(line) < 2:
            continue

        if line[1] == "ATOMS":
            foundBoxBounds = False
            boundsIndex = 0
            foundAtoms = True
            continue

        if line[1] == "BOX":
            foundBoxBounds = True
            continue

        if(foundBoxBounds):
            bounds[boundsIndex] = [float(line[0]), float(line[1])]
            boundsIndex += 1

        if line[1] == "TIMESTEP":
            if(not firstTimestep):
                snapshot = Snapshot(residues)
                snapshots.append(snapshot)
                residues = []

            foundAtoms = False
            residuePosition = 0
            continue 

        if(foundAtoms):
            firstTimestep = False
            if(int(line[1]) == 1):
                x = (bounds[0][1]-bounds[0][0])*float(line[2])+bounds[0][0]
                y = (bounds[1][1]-bounds[1][0])*float(line[3])+bounds[1][0]
                z = (bounds[2][1]-bounds[2][0])*float(line[4])+bounds[2][0]
                residues.append(Residue(residuePosition, x, y, z))
                residuePosition += 1

    snapshot = Snapshot(residues)
    snapshots.append(snapshot)

    return snapshots

def readAllTrajectories(metadataFile):
    f = open(metadataFile, 'r')
    for line in f:
        line=line.split()
        trajectory = Trajectory(line[0],line[1])
        trajectories.append(trajectory)

def readAllFreeEnergies(freeEnergyDirectory):
    print "Reading free energy files in: " + str(freeEnergyDirectory)
    path = freeEnergyDirectory
    listing = os.listdir(path)
    mintemp = 99999
    for infile in listing:
        temp = int(infile)/10
        temperaturearray.append(int(infile)/10)
        if temp < mintemp:
            mintemp = temp

    temperaturearray.sort()

    for i in range(len(temperaturearray)):
        fofqandt.append([])

    for infile in listing:
        f = open(freeEnergyDirectory+'/'+infile, 'r')
        for line in f:
            line=line.split()
            q = line[0]
            fofq = line[1]
            index = int(infile)/10-mintemp
            fofqandt[index].append([q, fofq])

    # fofqandt[temperature-mintemp][q/fofqindex][0=qvalue,1=fofqvalue]
    return mintemp

def FofQandT(Q,T):
    fofq = fofqandt[T-mintemp]
    index = 0
    for i in range(len(fofq)):
        tempQ = fofq[i][0]
        if float(tempQ) > float(Q):
            index = i
            break
    
    interpFvalue = float(fofq[index-1][1]) + ((float(fofq[index][1])-float(fofq[index-1][1]))/(float(fofq[index][0])-float(fofq[index-1][0])))*(float(Q)-float(fofq[index-1][0]))

    return interpFvalue

def boltzmannFactor(snapshot,temperature):
    exponent = (float(snapshot.energy)-float(FofQandT(snapshot.Q,temperature)))/(kb*float(temperature))
    return math.exp(exponent)

def calcZ(temperature):
    Z = 0.0
    for i in range(len(trajectories)):
        for j in range(len(trajectories[i].snapshots)):
            Z += boltzmannFactor(trajectories[i].snapshots[j],temperature)

    return Z

def findAllUstates():
    for i in range(len(trajectories)):
        for j in range(len(trajectories[i].snapshots)):
            tempustatecode = trajectories[i].snapshots[j].ustate.code
            if len(microstatecodes) == 0:
                microstatecodes.append(tempustatecode)
            for k in range(len(microstatecodes)):
                if tempustatecode == microstatecodes[k]:
                    continue
                if tempustatecode != microstatecodes[k] and k == len(microstatecodes)-1:
                    microstatecodes.append(tempustatecode)

    return microstatecodes

def sortAllSnapshots():
    sortedsnapshots = []
    for i in range(len(microstatecodes)):
        sortedsnapshots.append([])
    for i in range(len(trajectories)):
        for j in range(len(trajectories[i].snapshots)):
            for k in range(len(microstatecodes)):
                if trajectories[i].snapshots[j].ustate.code == microstatecodes[k]:
                    sortedsnapshots[k].append(trajectories[i].snapshots[j])
                    break

    return sortedsnapshots

def calculateUstateProbabilities():
    for i in range(len(sortedsnapshots)):
        ustateprob = 0.0
        for j in range(len(sortedsnapshots[i])):
            ustateprob += sortedsnapshots[i][j].probability

        microstateprobabilities.append(ustateprob)

    return microstateprobabilities

def calculateUstateFreeEnergies(temperature):
    microstatefreeenergies = []
    for i in range(len(microstatecodes)):
        microstatefreeenergies.append(-kb*temperature*math.log(microstateprobabilities[i]))

    return microstatefreeenergies

#############
# Libraries #
#############
import math
import os

#########
# Files #
#########
# The metadata file, containing links to the dump files and Qw/Potential energy
# files. The format is: dumpfile qw-pot-file
metadataFile = './metadata'
# Foldon file: each line contains the residues in a foldon
# each residue in the protein should be included once and only once
foldonFile = './foldons'
# The dump file (LAMMPS format) of the native structure coordinates
nativeDumpFile = './dump.native'
# The directory containing free energy files in the output format of UltimateWHAM
freeEnergyFileDirectory = './freeenergyfiles/'

#############
# Constants #
#############
# Boltzmann constant
kb = 0.001987 # kcal/mol/K

##############
# Parameters #
##############
# Native contact threshold, in Angstroms, to be applied to Ca-Ca distances:
nativeContactThreshold = 10
# Minimum sequence separation for two residues in contact
minSeqSep = 3
# Foldon foldedness threshold
foldonThreshold = 0.8

########################
# Variables and arrays #
########################
trajectories = []     # a list of Trajectory objects
temperaturearray = [] # a list of temperature values for which free energies are available
fofqandt = []         # a list containing all values of the free energy at all temperatures
microstatecodes = []  # a list containing all (sampled) microstate codes
sortedsnapshots = []  # a 2D list containing all snapshots corresponding to a given microstate code
microstateprobabilities = [] # a list containing all (sampled) microstate probabilities
microstatefreeenergies = []  # a list containing all (sampled) microstate free energies
Z = 0.0 # the partition sum

################
# Main program #
################
# read in native coordinates for the purposes of computing contacts
nativeSnapshot = readDumpFile(nativeDumpFile)[0]
# read in the foldon definitions
foldons = readFoldonFile(foldonFile)
# read all the trajectory information from the metadata file
readAllTrajectories(metadataFile)
# assign microstate to all snapshots
for i in range(len(trajectories)):
    trajectories[i].assignAllUstates()
# read in all the free energy information, return minimum temperature (for indexing purposes)
mintemp = readAllFreeEnergies(freeEnergyFileDirectory)
# calculate partition function for a certain temperature
Z = calcZ(500)
# assign probabilities to all snapshots in all trajectories for a certain temperature
for i in range(len(trajectories)):
    trajectories[i].assignAllProbabilities(500)
# find all the microstates present
microstatecodes = findAllUstates()
# sort all snapshots
sortedsnapshots = sortAllSnapshots()
# calculate microstate probabilities (the temperature is implied by the snapshot probabilities)
microstateprobabilities = calculateUstateProbabilities()
# calculate microstate free energies for a certain temperature
microstatefreeenergies = calculateUstateFreeEnergies(500)





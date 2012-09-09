#########################
# Classes and Functions #
#########################
class Foldon:
    numRes = 0
    residuelist = []
    numNativeContacts = 0
    nativefoldoncontactlist = []

    def __init__(self,residuelist):
        self.residuelist = residuelist
        self.numRes = len(self.residuelist)
        self.numNativeContacts = self.calculateNumNativeContacts(nativeSnapshot)
        self.nativefoldoncontactlist = []
        self.samplingtemperature = 0.0

    def display(self):
        print self.residuelist

    def calculateNumNativeContacts(self, snapshot):
        numNativeContacts = 0.0
        for nativepairindex in range(len(self.nativefoldoncontactlist)):
            residueindex1 = self.nativefoldoncontactlist[nativepairindex][0]
            residueindex2 = self.nativefoldoncontactlist[nativepairindex][1]
            residue1 = snapshot.residues[residueindex1]
            residue2 = snapshot.residues[residueindex2]
            if qType == 'QC':
                if residue1.isInContactWith(residue2):
                    numNativeContacts += 1
            elif qType == 'QW':
                qvalue = numpy.exp(-pow(residue1.distanceTo(residue2)-nativedistances[residue1.pos][residue2.pos],2)/(2*pow(abs(residue1.pos-residue2.pos),0.3)))
                numNativeContacts += qvalue

        return numNativeContacts

    def isFoldedIn(self, snapshot):
        currentContacts = self.calculateNumNativeContacts(snapshot)
        if float(currentContacts)/float(len(self.nativefoldoncontactlist)) > foldonThreshold:
            return True
        else:
            return False
    
    def containsResidue(self, residue):
        containsresidue = False
        for i in range(len(self.residuelist)):
            if residue == self.residuelist[i]:
                containsresidue = True
                break

        return containsresidue

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

    def isNativeContactWith(self,otherresidue):
        distance = math.sqrt(pow(self.x-otherresidue.x,2)+pow(self.y-otherresidue.y,2)+pow(self.z-otherresidue.z,2))
        if abs(self.pos - otherresidue.pos) >= minSeqSep and distance < nativeContactThreshold:
            return True
        else:
            return False

    def isInContactWith(self,otherresidue):
        distance = math.sqrt(pow(self.x-otherresidue.x,2)+pow(self.y-otherresidue.y,2)+pow(self.z-otherresidue.z,2))
        if abs(self.pos - otherresidue.pos) >= minSeqSep and distance < contactFactor*nativedistances[self.pos][otherresidue.pos]:
            return True
        else:
            return False

    def distanceTo(self,otherresidue):
        distance = math.sqrt(pow(self.x-otherresidue.x,2)+pow(self.y-otherresidue.y,2)+pow(self.z-otherresidue.z,2))
        return distance

class Snapshot:
    numRes = 0
    residues = []
    biasingenergy = 0.0
    internalenergy = 0.0
    reducedenergy = 0.0
    Q = 0.0
    ustate = ustate("")
    
    def __init__(self, residues):
        self.numRes = len(residues)
        self.residues = residues                                           
        self.biasingenergy = 0.0
        self.internalenergy = 0.0
        self.reducedenergy = 0.0
        self.Q = 0.0
        self.ustate = 0
        self.ustate = ustate("")

    def display(self):
        print "Ca coordinates:"
        for i in range(self.numRes):
            print "residue " + str(i) + ": " + str(self.residues[i].x) + " " + str(self.residues[i].y) + " " + str(self.residues[i].z)
            
        print "Biasing energy: " + str(self.biasingenergy)
        print "Internal energy: " + str(self.internalenergy)
        print "Reduced energy: " + str(self.reducedenergy)
        print "Qvalue: " + str(self.Q)
        print "Microstate code:"
        self.ustate.display()

    def assignUstate(self):
        microstatecode = ""
        for i in range(len(foldons)):
            if foldons[i].isFoldedIn(self):
                microstatecode += '1'
            else:
                microstatecode += '0'

        self.ustate = ustate(microstatecode)
        return microstatecode

    def calculatePairwiseDistances(self):
        pairwisedistances = []
        for i in range(self.numRes):
            pairwisedistances.append([])
            for j in range(self.numRes):
                xi = self.residues[i].x
                yi = self.residues[i].y
                zi = self.residues[i].z
                xj = self.residues[j].x
                yj = self.residues[j].y
                zj = self.residues[j].z
                dist = numpy.sqrt(pow(xi-xj,2)+pow(yi-yj,2)+pow(zi-zj,2))
                pairwisedistances[i].append(dist)

        return pairwisedistances
    
    def computeReducedEnergy(self):
        return (float(self.internalenergy)+float(self.biasingenergy))/kb*float(self.samplingtemperature)
                           
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
        lineindex = 0
        for dataline in ssdata:
            dataline = dataline.split()
            if dataline == "" or dataline[0] == "#":
                continue
            if lineindex % snapshotFreq == 0:
                self.snapshots[snapshotindex].Q = dataline[1]
                self.snapshots[snapshotindex].internalenergy = dataline[2]
                self.snapshots[snapshotindex].biasingenergy = dataline[3]
                self.snapshots[snapshotindex].samplingtemperature = dataline[4]
                self.snapshots[snapshotindex].reducedenergy = self.snapshots[snapshotindex].computeReducedEnergy()
                snapshotindex += 1
            lineindex += 1

    def display(self):
        for i in range(len(self.snapshots)):
            self.snapshots[i].display()
    
    def assignAllUstates(self):
        for i in range(len(self.snapshots)):
            self.snapshots[i].assignUstate()

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
                residues.append(int(line[i])-1)

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

    snapshotIndex = 0

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
                if snapshotIndex % snapshotFreq == 0:
                    snapshot = Snapshot(residues) # create snapshot with residue coordinates
                    snapshot.assignUstate()       # assign microstate based on coordinates
                    del snapshot.residues         # delete residue coordinates to save memory
                    snapshots.append(snapshot)    # append snapshot to snapshots list
                snapshotIndex += 1
                residues = []                 # zero out residues for next snapshot        

            foundAtoms = False
            residuePosition = 0
            continue 

        if(foundAtoms):
            firstTimestep = False
            if atomType == 'CA':
                if(int(line[1]) == 1):
                    x = (bounds[0][1]-bounds[0][0])*float(line[2])+bounds[0][0]
                    y = (bounds[1][1]-bounds[1][0])*float(line[3])+bounds[1][0]
                    z = (bounds[2][1]-bounds[2][0])*float(line[4])+bounds[2][0]
                    residues.append(Residue(residuePosition, x, y, z))
                    residuePosition += 1
            elif atomType == 'CB':
                if(int(line[1]) == 4 or int(line[1]) == 5):
                    x = (bounds[0][1]-bounds[0][0])*float(line[2])+bounds[0][0]
                    y = (bounds[1][1]-bounds[1][0])*float(line[3])+bounds[1][0]
                    z = (bounds[2][1]-bounds[2][0])*float(line[4])+bounds[2][0]
                    residues.append(Residue(residuePosition, x, y, z))
                    residuePosition += 1
            else:
                print "Wrong atom type: " + str(atomType)
                sys.exit()

    if snapshotIndex % snapshotFreq == 0:
        snapshot = Snapshot(residues) # create snapshot with residue coordinates 
        snapshot.assignUstate()       # assign microstate based on coordinates   
        del snapshot.residues         # delete residue coordinates to save memory
        snapshots.append(snapshot)    # append snapshot to snapshots list        
                                  
    return snapshots

def readNativeDumpFile(dumpFile):
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
            if atomType == 'CA':
                if(int(line[1]) == 1):
                    x = (bounds[0][1]-bounds[0][0])*float(line[2])+bounds[0][0]
                    y = (bounds[1][1]-bounds[1][0])*float(line[3])+bounds[1][0]
                    z = (bounds[2][1]-bounds[2][0])*float(line[4])+bounds[2][0]
                    residues.append(Residue(residuePosition, x, y, z))
                    residuePosition += 1
            elif atomType == 'CB':
                if(int(line[1]) == 4 or int(line[1]) == 5):
                    x = (bounds[0][1]-bounds[0][0])*float(line[2])+bounds[0][0]
                    y = (bounds[1][1]-bounds[1][0])*float(line[3])+bounds[1][0]
                    z = (bounds[2][1]-bounds[2][0])*float(line[4])+bounds[2][0]
                    residues.append(Residue(residuePosition, x, y, z))
                    residuePosition += 1
            else:
                print "Wrong atom type: " + str(atomType)
                sys.exit()

    snapshot = Snapshot(residues)
    snapshots.append(snapshot)

    return snapshots

def readAllTrajectories(metadataFile):
    f = open(metadataFile, 'r')
    for line in f:
        line=line.split()
        print "Creating trajectory for " + line[0] + " ..."
        trajectory = Trajectory(line[0],line[1])
        trajectories.append(trajectory)
        samplingtemperature.append(float(line[2]))
        Kbias.append(float(line[3]))
        biasing_value.append(float(line[4]))

def findAllUstates():
    microstatecodes = []
    binmicrostatecodes = []
    for k in range(K):
        for n in range(N_k[k]):
            tempustatecode = ustate_kn[k,n]
            microstatecodes.append(tempustatecode)

    microstatecodes = list(set(microstatecodes))
    for i in range(len(microstatecodes)):
        binmicrostatecodes.append(binaryfoldonstate(microstatecodes[i]))

    return (microstatecodes, binmicrostatecodes)

def ustateDistance(code1,code2):
    distance = 0
    for i in range(len(code1)):
        distance += abs(int(code1[i])-int(code2[i]))

    return distance

def areConnected(code1,code2):
    if(ustateDistance(code1,code2)==1):
        return True
    else:
        return False

def calculateRateMatrix():
    ratematrix = numpy.zeros((len(microstatecodes),len(microstatecodes)),dtype=numpy.float64)

    for i in range(len(microstatecodes)):
        for j in range(len(microstatecodes)):
            rate = 0.0
            if areConnected(binmicrostatecodes[i],binmicrostatecodes[j]):
                rate = k0
                if microstatefreeenergies[i] < microstatefreeenergies[j]:
                    rate = k0*numpy.exp(-(microstatefreeenergies[j]-microstatefreeenergies[i]))
                
            ratematrix[i][j] = rate

    for i in range(len(microstatecodes)):
        total = 0.0
        for j in range(len(microstatecodes)):
            total += ratematrix[i,j]

        ratematrix[i,i] = -total

    return numpy.transpose(ratematrix)

def calculateEigenvectorsEigenvalues():
    eigenvalues, eigenvectors=LA.eig(ratematrix)
    perm = numpy.argsort(-eigenvalues)  # sort in descending order
    for i in range(len(eigenvalues)):
        if abs(numpy.real(eigenvalues[i])) > 0.0:
            eigenvaluesout.write("%s\n" % str(numpy.imag(eigenvalues[i])/numpy.real(eigenvalues[i])))
        else:
            eigenvaluesout.write("0.0\n")
    return numpy.real(eigenvalues), numpy.real(eigenvectors), numpy.real(eigenvalues[perm]), numpy.real(numpy.transpose(eigenvectors[:, perm]))

def findNativeContacts(nativesnapshot):
    numNativeContacts = 0
    for residue1 in range(len(nativesnapshot.residues)):
        for residue2 in range(residue1,len(nativesnapshot.residues)):
            if nativesnapshot.residues[residue1].isNativeContactWith(nativesnapshot.residues[residue2]):
                nativecontactlist.append([residue1,residue2])

def sortNativeContactsByFoldon():
    if includeInterfaceContacts:
        for foldon in foldons:
            for nativepairindex in range(len(nativecontactlist)):
                residue1 = nativecontactlist[nativepairindex][0]
                residue2 = nativecontactlist[nativepairindex][1]
                if foldon.containsResidue(residue1) or foldon.containsResidue(residue2):
                    foldon.nativefoldoncontactlist.append([residue1, residue2])

            print "Foldon definition:"
            foldon.display()
            print "Number of native contacts:"
            print len(foldon.nativefoldoncontactlist)
            print "Native contacts being monitored:"
            print foldon.nativefoldoncontactlist

    else:
        for foldon in foldons:
            for nativepairindex in range(len(nativecontactlist)):
                residue1 = nativecontactlist[nativepairindex][0]
                residue2 = nativecontactlist[nativepairindex][1]
                if foldon.containsResidue(residue1) and foldon.containsResidue(residue2):
                    foldon.nativefoldoncontactlist.append([residue1, residue2])
                    
            print "Foldon definition:"
            foldon.display()
            print "Number of native contacts:"
            print len(foldon.nativefoldoncontactlist)
            print "Native contacts being monitored:"
            print foldon.nativefoldoncontactlist

def floatRange(a, b, inc):
    try: x = [float(a)]
    except: return False
    for i in range(1, int(math.ceil((b - a ) / inc))):
        x. append(a + i * inc)
    
    return x

def calculateRank(ustatecode):
    rank = 0;
    for i in range(len(ustatecode)):
        rank += int(ustatecode[i])

    return rank

def readNativeDistances(nativedistancefile):
    nativedistances = []
    lineIndex = 0
    f = open(nativedistancefile, 'r')
    for line in f:
        line = line.split()
        nativedistances.append([])
        for i in range(len(line)):
            nativedistances[lineIndex].append(float(line[i]))
        lineIndex += 1

    return nativedistances

def findMinRank(ranks):
    minrank = 999999
    for rank in ranks:
        if int(rank[1]) < minrank:
            minrank = int(rank[1])

    return minrank

def findMaxRank(ranks):
    maxrank = 0
    for rank in ranks:
        if int(rank[1]) > maxrank:
            maxrank = int(rank[1])

    return maxrank

def findUstateStabilityGap(minrank,maxrank):
    minrankminfreeenergy = 9999999.9
    maxrankminfreeenergy = 9999999.9
    microstateindex = 0
    for microstate in binmicrostatecodes:
        if calculateRank(microstate) == minrank:
            if microstatefreeenergies[microstateindex] < minrankminfreeenergy:
                minrankminfreeenergy = microstatefreeenergies[microstateindex]

        if calculateRank(microstate) == maxrank:
            if microstatefreeenergies[microstateindex] < maxrankminfreeenergy:
                maxrankminfreeenergy = microstatefreeenergies[microstateindex]

        microstateindex += 1

    return maxrankminfreeenergy - minrankminfreeenergy

def findUnfoldedAndFoldedStates():
    return len(foldons)*'0', len(foldons)*'1'

def subsampleTrajectories():
    for k in range(K):
        N = len(trajectories[k].snapshots)
        for n in range(N):
            snapshot = trajectories[k].snapshots[n]
            qw_kt[k,n] = snapshot.Q
            reducedU_kt[k,n] = snapshot.reducedenergy
            U_kt[k,n] = snapshot.internalenergy
            UB_kt[k,n] = snapshot.biasingenergy
            ustate_kt[k,n] = decimalfoldonstate(snapshot.ustate.code)
        # Extract timeseries.
        A_t = qw_kt[k,:]
        # Compute statistical inefficiency.
        try:
            g = timeseries.statisticalInefficiency(A_t)
        except Exception as e:
            print str(e)
            print A_t
        
        # Subsample data.
        if subsample:
            indices = timeseries.subsampleCorrelatedData(A_t, g=g)
        else:
            indices = timeseries.subsampleCorrelatedData(A_t, g=1)
        N_uncorr = len(indices) # number of uncorrlated samples
        print "k = %5d : g = %.1f, N_uncorr = %d" % (k, g, N_uncorr)
        qw_kn[k,0:N_uncorr] = qw_kt[k,indices]
        U_kn[k,0:N_uncorr] = U_kt[k,indices]
        reducedU_kn[k,0:N_uncorr] = reducedU_kt[k,indices]
        UB_kn[k,0:N_uncorr] = UB_kt[k,indices]
        ustate_kn[k,0:N_uncorr] = ustate_kt[k,indices]
        N_k[k] = N_uncorr # number of uncorrelated samples

def decimalfoldonstate(binaryfoldonstate):
    decimalfoldonstate = 0
    for i in range(len(binaryfoldonstate)):
        decimalfoldonstate += pow(2,i)*int(binaryfoldonstate[i])
    return decimalfoldonstate

def binaryfoldonstate(decimalfoldonstate):
    return bin(decimalfoldonstate)[2:].zfill(len(foldons))

def computeReducedEnergies():
    # Compute reduced potentials from all simulations in all thermodynamic states.
    print "Computing reduced potentials..."
    u_kln = numpy.zeros([K,K,N_max], numpy.float32) # u_kln[k,l,n] is reduced biased potential of uncorrelated snapshot n from simulation k in thermodynamic state l
    for k in range(K):
        N = N_k[k] # number of uncorrelated snapshots
        # Compute reduced potential in all other temperatures and biasing potentials indexed by l.
        for l in range(K):
            kT = kb * samplingtemperature[l] # thermal energy (in kcal/mol)
            beta = 1.0 / kT # inverse temperature (in 1 / (kcal/mol))
            Ubias = (Kbias[l]/2.0) * (qw_kn[k,0:N] - biasing_value[l])**2
            u_kln[k,l,0:N] = beta * (U_kn[k,0:N] + Ubias)
    
    return u_kln

#############
# Libraries #
#############
import math
import os
import numpy
numpy.set_printoptions(threshold=numpy.nan)
import sys
from numpy import linalg as LA
import cPickle
import gc
import time as timefunctions

# pymbar imports
import timeseries
import pymbar

#########
# Files #
#########
# The metadata file, containing links to the dump files and Qw/Potential energy
# files. The format is: dumpfile qw-pot-file
metadataFile = './metadatashort'
# Foldon file: each line contains the residues in a foldon
# each residue in the protein should be included once and only once
foldonFile = './foldonsfullank'
# The dump file (LAMMPS format) of the native structure coordinates
nativeDumpFile = './dump.native'
# Overall rate file
overallRateFile = './overallrates'
# Trajectories pickle file
trajectoriespicklefile = './trajectories.pkl'
# MBAR pickle file
mbarpicklefile = './mbar.pkl'
# Microstate ranks file prefix
microstateInfoFilePrefix = './microstateinfo'
# Eigenvalue debugging
eigenvaluesoutfile = './eigenvalueratio.dat'
eigenvaluesout = open(eigenvaluesoutfile,'w')

#############
# Constants #
#############
# Boltzmann constant
kb = 0.001987 # kcal/mol/K

##############
# Parameters #
##############
# Which type of Q calculation do you want to use? QC (a contact Q) or QW (a sum of gaussians)?
qType = 'QW'
# Which atom type to consider for contacts? CA or CB?
atomType = 'CA'
# Native contact threshold, in Angstroms, to be applied to Ca-Ca distances:
nativeContactThreshold = 8.0
# Contact factor: two residues are in contact if their distances is less than contactFactor*nativedistance (only used for qType = 'QC')
contactFactor = 1.2
# Include interface contacts? If False, only those native contacts for residues within the same foldon will be used
includeInterfaceContacts = False
# Minimum sequence separation for two residues in contact
minSeqSep = 3
# Foldon foldedness threshold
foldonThreshold = 0.6
# Downhill rate
k0 = 1000000
# Minimum temperature for computing overall rate
starttemp = 250
# Maximum temperature for computing overall rate
endtemp = 300
# Temperature incremement
tempinc = 1
# Temperature array
temperatures = range(starttemp,endtemp+1,tempinc)
# output microstate information for each temperature?
outputMicrostateInfo = True
# The frequency at which to accept snapshots from the dump file
snapshotFreq = 1 # WARNING: Because of recent changes, using snapshotFreq != 1 may break something (because of the subsampling)
# Always output rate matrix?
alwaysOutputRateMatrix = True
# Calculate time evolution and fluxes?
calculateTimeEvolutionAndFluxes = True
# Automatically determine folding or unfolding simulation?
autoDetermineFoldingOrUnfolding = True # folding if below the folding temperature, unfolding if above
# Folding or unfolding simulation for calulating fluxes, not used if autoDetermineFoldingOrUnfolding = True
foldingSimulation = True # if false, assume unfolding
# Time range and time step for calculating evolution and fluxes
startTime = 0 # time evolution will start at this time
endTime = 0.0001 # time evolution will end at this time and the integrated flux will be calculated
timeStep = 0.00001 # the concentration of the states will be calculated this often
# Automatically determine time evolution interval?
autoDetermineEvolutionInterval = True # if True, startTime, endTime and timeStep (above) are not used
numTimeSteps = 100 # number of points to plot on time evolution if the interval is automatically determined
numRelaxationTimes = 5 # number of relaxation times (negative inverse of the smallest nonzero eigenvalue) to integrate out to
# Show graphical evolution of concentration of states?
graphicalEvolution = False
# Calculate equilibrium flux? Otherwise, calculate net flux at end of time evolution
calculateEquilibriumFlux = True

# Time saving variables
# read trajectories from metadata? if not, load trajectories.pkl
readTrajectoriesFromMetadata = False
# Initialize MBAR? Otherwise, load from pickle file
initializeMBAR = True

# MBAR parameters
subsample = False

########################
# Variables and arrays #
########################
trajectories = []     # a list of Trajectory objects
temperaturearray = [] # a list of temperature values for which free energies are available
fofqandt = []         # a list containing all values of the free energy at all temperatures
microstatecodes = []  # a list containing all (sampled) microstate codes
binmicrostatecodes = []  # a list containing all (sampled) microstate codes in binary
microstateprobabilities = [] # a list containing all (sampled) microstate probabilities
microstatefreeenergies = []  # a list containing all (sampled) microstate free energies
microstateuncertainties = []  # a list containing all (sampled) microstate free energy uncertainties
Z = 0.0 # the partition sum
ratematrix = []   # the rate matrix
eigenvectors = [] # eigenvectors of the rate matrix
eigenvalues = []  # eigenvalues of the rate matrix 
sortedeigenvectors = [] # eigenvectors of the rate matrix sorted by eigenvalue
sortedeigenvalues = []  # eigenvalues of the rate matrix sorted in decreasing order
overallrates = [] # stores all temperature/overall rate pairs calculated
nativecontactlist = [] # stores all native contact pairs
foldons = [] # a list of all foldons
heatcapacity = [] # heat capacity array
foldingtemperature = 0.0 # folding temperature
unfoldedQ = 0.0 # Q of the unfolded basin
foldedQ = 0.0   # Q of the folded basin
ustateinfo = [] # a list of microstate information for each temperature
nativedistances = [] # a list of the CA-CA distances
minrank = 0 # minimum rank of a sampled microstate
maxrank = 0 # maximum rank of a samples microstate
ustatestabilitygap = 0.0 # free energy gap between microstate with minimum rank, minimum
                         # free energy and maximum rank, minimum free energy
fluxes = [] # fluxes between states
initialconcentrations = [] # initial concentration of states
unfoldedstate = "" # unfolded microstate code
foldedstate = ""   # folded microstate code
concentrations = [] # time dependence of the concentrations

# MBAR related variables
K = 18 # number of simulations/trajectories
N_max = 10001 # maximum number of snapshots per trajectory
N_k = numpy.zeros([K], numpy.int32) 
qw_kt = numpy.zeros([K,N_max], numpy.float32) 
U_kt = numpy.zeros([K,N_max], numpy.float32)
reducedU_kt = numpy.zeros([K,N_max], numpy.float32)
reducedU_kn = numpy.zeros([K,N_max], numpy.float32)
UB_kt = numpy.zeros([K,N_max], numpy.float32)
ustate_kt = numpy.zeros([K,N_max], numpy.int32)
qw_kn = numpy.zeros([K,N_max], numpy.float32)
U_kn = numpy.zeros([K,N_max], numpy.float32)
UB_kn = numpy.zeros([K,N_max], numpy.float32)
ustate_kn = -1 * numpy.ones([K,N_max], numpy.int32) 
Kbias = []
samplingtemperature = []
biasing_value = []

################
# Main program #
################
# find folding temperature
foldingtemperature = 267
print "Folding temperature: " + str(foldingtemperature)
# read in native coordinates for the purposes of computing contacts
print "Reading native dump file..."
nativeSnapshot = readNativeDumpFile(nativeDumpFile)[0]
# calculate native distances
print "Calculating native distances..."
nativedistances = nativeSnapshot.calculatePairwiseDistances()
# Find native contacts
print "Finding native contacts..."
findNativeContacts(nativeSnapshot)
# read in the foldon definitions
print "Reading foldon file..."
foldons = readFoldonFile(foldonFile)
sortNativeContactsByFoldon()
if readTrajectoriesFromMetadata:
    # read all the trajectory information from the metadata file and assign all microstates
    print "Reading all trajectories and assigning microstates..."
    readAllTrajectories(metadataFile)
    # Pickle trajectories list for reading later
    print "Pickling trajectories..."
    cPickle.dump(trajectories, open(trajectoriespicklefile, 'wb')) 
    cPickle.dump(Kbias, open('./kbias.pkl', 'wb'))
    cPickle.dump(samplingtemperature, open('./samplingtemperature.pkl', 'wb'))
    cPickle.dump(biasing_value, open('./biasingvalue.pkl', 'wb'))

else:
    # load trajectories from existing Pickle file
    print "Loading pickled trajectories..."
    trajectories = cPickle.load(open(trajectoriespicklefile, 'rb'))
    Kbias = cPickle.load(open('./kbias.pkl', 'rb'))
    samplingtemperature = cPickle.load(open('./samplingtemperature.pkl', 'rb'))
    biasing_value = cPickle.load(open('./biasingvalue.pkl', 'rb'))

# subsample data because of time correlations in trajectories
subsampleTrajectories()

# find all the microstates present
print "Finding all microstates..."
(microstatecodes, binmicrostatecodes) = findAllUstates()
print "Microstate codes: " + str(binmicrostatecodes)
# finding folded and unfolded states
unfoldedstate, foldedstate = findUnfoldedAndFoldedStates()
# Make sure all bins are populated.
print "Counting number of samples per bin..."
nbins = len(microstatecodes)
bin_counts = numpy.zeros([nbins], numpy.int32)
for i in range(nbins):
    code = microstatecodes[i]
    bin_counts[i] = (ustate_kn == code).sum()
print "Number of samples per microstate code: " + str(bin_counts)
print "Computing reduced energies of all samples in all states..."
u_kln = computeReducedEnergies()

# Initialize or load MBAR.
if initializeMBAR:
    print "Initializing MBAR for the calculation of free energies..."
    mbar = pymbar.MBAR(u_kln, N_k)
    cPickle.dump(mbar, open(mbarpicklefile, 'wb')) 
else:
    # load MBAR object
    print "Loading MBAR object..."
    mbar = cPickle.load(open(mbarpicklefile, 'rb'))

# Bin data
bin_kn = -1 * numpy.ones([K,N_max], numpy.int32) # bin_kn[k,n] is bin index of sample n from simulation k; otherwise -1
for k in range(K):
    N = N_k[k]
    # Compute bin assignment.
    for n in range(N):
        bin_kn[k,n] = microstatecodes.index(ustate_kn[k,n])

# loop over all desired (extrapolated) temperatures, compute overall rate
for temperature in temperatures:
    print "Microstate codes: " + str(binmicrostatecodes)
    print "Calculating rate for temperature: " + str(temperature)
    # Compute perturbed reduced potential at temperature of interest in absence of biasing potential.
    print "Computing perturbed reduced potential at temperature of interest..."
    u_kn = numpy.zeros([K,N_max], numpy.float32) # u_kn[k,n] is the unbiased reduced potential energy of snapshot n of umbrella simulation k at conditions of interest
    kT = kb * temperature
    beta = 1.0 / kT # reduced temperature
    for k in range(K):
        N = N_k[k]
        u_kn[k,0:N] = beta * U_kn[k,0:N] # unbiased reduced potential at desired temperature
    # calculate microstate free energies for a certain temperature
    print "Calculating microstate free energies..."
    (microstatefreeenergies, microstateuncertainties) = mbar.computePMF(u_kn, bin_kn, nbins)  # Compute PMF in unbiased potential (in units of kT).
    print "Microstate free energies: " + str(microstatefreeenergies)
    # calculate rate matrix, given a temperature
    print "Calculating rate matrix..."
    ratematrix = calculateRateMatrix()
    if(alwaysOutputRateMatrix):
        numpy.savetxt('ratematrix'+str(temperature),ratematrix)
    print "Rate Matrix: \n" + str(ratematrix)
    # calculate and append microstate info
    info = []
    for i in range(len(microstatecodes)):
        info.append([binmicrostatecodes[i],calculateRank(binmicrostatecodes[i]),microstatefreeenergies[i]])
    ustateinfo.append(info)
    minrank = findMinRank(info)
    maxrank = findMaxRank(info)
    ustatestabilitygap = findUstateStabilityGap(minrank,maxrank)
    # calculate eigenvalues and eigenvectors
    print "Calculating eigenvalues and eigenvectors"
    eigenvalues, eigenvectors, sortedeigenvalues, sortedeigenvectors = calculateEigenvectorsEigenvalues()
    print "Eigenvalues: \n" + str(sortedeigenvalues)
    print "Eigenvectors: \n " + str(sortedeigenvectors)
    overallrates.append([temperature,ustatestabilitygap,numpy.log(-1*sortedeigenvalues[1]),numpy.log(-1*sortedeigenvalues[2])])
    #############################################
    # Begin time evolution and flux calculation #
    #############################################
    if(calculateTimeEvolutionAndFluxes):
        fluxes = numpy.zeros((len(microstatecodes),len(microstatecodes)))
        if(autoDetermineEvolutionInterval):
            startTime = 0.0
            endTime = numRelaxationTimes*(1/-sortedeigenvalues[1])
            timeStep = (endTime-startTime)/numTimeSteps
        else:
            numTimeSteps = int((endTime-startTime)/timeStep)
        print "Starting time evolution and flux calculation..."
        print "Start time: %s\n" % str(startTime)
        print "End time: %s\n" % str(endTime)
        print "Time step: %s\n" % str(timeStep)
        concentrations = numpy.zeros((len(microstatecodes),numTimeSteps))
        if(autoDetermineFoldingOrUnfolding):
            print "Determining if this will be a folding or unfolding calculation based on folding temperature..."
            if(temperature < foldingtemperature):
                print "Starting folding calculation..."
                foldingSimulation = True
            else:
                print "Starting unfolding calculation..."
                foldingSimulation = False
        print "Setting initial concentrations..."
        initialconcentrations = numpy.zeros(len(microstatecodes))
        if(foldingSimulation):
            initialconcentrations[binmicrostatecodes.index(unfoldedstate)] = 1.0
        else:
            initialconcentrations[binmicrostatecodes.index(foldedstate)] = 1.0
        initialconcentrations = initialconcentrations[:,numpy.newaxis]
        print "Initial concentrations: \n%s\n" % str(initialconcentrations)
        coefficients = numpy.dot(LA.inv(eigenvectors),initialconcentrations)
        print "Coefficents: \n%s\n" % str(coefficients)
        print "Calculating time evolution..."
        timeIndex = 0
        for time in floatRange(startTime,endTime,timeStep):
            if(graphicalEvolution):
                print "Concentration of states:"
            for state in range(len(microstatecodes)):
                total = 0
                for eigenvalueindex in range(len(eigenvalues)):
                    total += coefficients[eigenvalueindex]*eigenvectors[state,eigenvalueindex]*numpy.exp(eigenvalues[eigenvalueindex]*float(time))
                concentrations[state,timeIndex] = total
                if(graphicalEvolution):
                    print str(binmicrostatecodes[state]) + ": " + int(concentrations[state,timeIndex]*100)*'#'
            timeIndex += 1
            if(graphicalEvolution):
                timefunctions.sleep(0.1)
            numpy.savetxt('concentrations'+str(temperature),concentrations.transpose())
            if(timeIndex >= numTimeSteps):
                break

        # calculate integrated flux at the end of the time evolution
        for state1 in range(len(microstatecodes)):
            for state2 in range(len(microstatecodes)):
                total = 0.0
                for eigenvalueindex in range(len(eigenvalues)):
                    ev = eigenvalues[eigenvalueindex]
                    # Skip the smallest eigenvalue because the equilibrium eigenvector
                    # doesn't contribute to the flux
                    if(ev == sortedeigenvalues[0]):
                        continue
                    c = coefficients[eigenvalueindex]
                    rm21 = ratematrix[state2][state1]
                    rm12 = ratematrix[state1][state2]
                    ev1 = eigenvectors[state1,eigenvalueindex]
                    ev2 = eigenvectors[state2,eigenvalueindex]
                    if(calculateEquilibriumFlux):
                        total += c/-ev*(rm21*ev1-rm12*ev2)
                    else:
                        total += c/-ev*(1-numpy.exp(ev*float(endTime)))*(rm21*ev1-rm12*ev2)

                fluxes[state1][state2] = total

        print "Fluxes: \n%s\n" % str(fluxes)
        numpy.savetxt('fluxes'+str(temperature),fluxes)
        cPickle.dump(fluxes, open('fluxes'+str(temperature)+'.pkl', 'wb')) 

f = open(overallRateFile, 'w')
f.write("# temperature stability-gap first-rate second-rate\n")
for i in range(len(overallrates)):
    f.write("%f %f %e %e \n" % (overallrates[i][0],overallrates[i][1],overallrates[i][2],overallrates[i][3]))
f.close()

if outputMicrostateInfo:
    for temperature in range(starttemp,endtemp+1,tempinc):
        f = open(microstateInfoFilePrefix + "." + str(temperature), 'w')
        f.write("# macrobasin rank free-energy-in-kT\n")
        ustateinformation = ustateinfo[int((temperature-starttemp)/tempinc)]
        for index in range(len(ustateinformation)):
            f.write("%s %d %f \n" % (ustateinformation[index][0], ustateinformation[index][1], ustateinformation[index][2]))
        f.close()

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
    energy = 0.0
    Q = 0.0
    ustate = ustate("")
    probability = 0.0
    
    def __init__(self, residues):
        self.numRes = len(residues)
        self.residues = residues                                           
        self.energy = 0.0
        self.Q = 0.0
        self.ustate = 0
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

    def assignUstate(self):
        microstatecode = ""
        for i in range(len(foldons)):
            if foldons[i].isFoldedIn(self):
                microstatecode += '1'
            else:
                microstatecode += '0'

        self.ustate = ustate(microstatecode)
        return microstatecode

    def assignProbability(self, extrapolatedtemperature):
        self.probability = boltzmannFactor(self,extrapolatedtemperature)/Z

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
                self.snapshots[snapshotindex].energy = dataline[2]
                self.snapshots[snapshotindex].samplingtemperature = dataline[3]
                snapshotindex += 1
            lineindex += 1

    def display(self):
        for i in range(len(self.snapshots)):
            self.snapshots[i].display()
    
    def assignAllUstates(self):
        for i in range(len(self.snapshots)):
            self.snapshots[i].assignUstate()

    def assignAllProbabilities(self,extrapolatedtemperature):
        for i in range(len(self.snapshots)):
            self.snapshots[i].assignProbability(extrapolatedtemperature)

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
        print "Creating trajectory for " + line[0] + "..."
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

def boltzmannFactor(snapshot,extrapolatedtemperature):
    exponent = (float(snapshot.energy)-float(FofQandT(snapshot.Q,extrapolatedtemperature)))/(kb*float(snapshot.samplingtemperature))
    return math.exp(exponent)

def calcZ(temperature):
    Z = 0.0
    for i in range(len(sortedsnapshots)):
        for j in range(len(sortedsnapshots[i])):
            Z += boltzmannFactor(sortedsnapshots[i][j],temperature)

    return Z

def findAllUstates():
    microstatecodes = []
    for i in range(len(trajectories)):
        for j in range(len(trajectories[i].snapshots)):
            tempustatecode = trajectories[i].snapshots[j].ustate.code
            microstatecodes.append(tempustatecode)

    microstatecodes = list(set(microstatecodes))

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

def calculateAverageQofMicrostates():
    averageqofmicrostates = []
    for i in range(len(microstatecodes)):
        averageqofmicrostates.append(0.0)
    for i in range(len(sortedsnapshots)):
        for j in range(len(sortedsnapshots[i])):
            for k in range(len(microstatecodes)):
                if sortedsnapshots[i][j].ustate.code == microstatecodes[k]:
                    averageqofmicrostates[k] += float(sortedsnapshots[i][j].Q)
                    break

    for i in range(len(microstatecodes)):
        averageqofmicrostates[i] /= len(sortedsnapshots[i])

    return averageqofmicrostates

def calculateUstateProbabilities():
    microstateprobabilities = []
    for i in range(len(sortedsnapshots)):
        ustateprob = 0.0
        for j in range(len(sortedsnapshots[i])):
            ustateprob += sortedsnapshots[i][j].probability

        microstateprobabilities.append(ustateprob)

    return microstateprobabilities

def calculateUstateFreeEnergies(temperature):
    microstatefreeenergies = []
    for i in range(len(microstatecodes)):
        microstatefreeenergies.append(-kb*float(temperature)*numpy.log(microstateprobabilities[i]))

    return microstatefreeenergies

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

def calculateRateMatrix(temperature):
    ratematrix = numpy.zeros((len(microstatecodes),len(microstatecodes)),dtype=numpy.float64)

    for i in range(len(microstatecodes)):
        for j in range(len(microstatecodes)):
            rate = 0.0
            if areConnected(microstatecodes[i],microstatecodes[j]):
                rate = k0
                if microstatefreeenergies[i] < microstatefreeenergies[j]:
                    rate = k0*numpy.exp(-(microstatefreeenergies[j]-microstatefreeenergies[i])/(kb*float(temperature)))
                
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
    return numpy.real(eigenvalues[perm]), numpy.real(numpy.transpose(eigenvectors[:, perm]))

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

def readHeatCapacityFile(heatcapacityfile):
    heatcapacity = []
    f = open(heatcapacityfile, 'r')
    for line in f:
        line = line.split()
        heatcapacity.append([line[0],line[1]])

    return heatcapacity

def findFoldingTemperature(heatcapacity):
    foldingtemperature = 0.0
    maxheatcapacity = 0.0
    for index in range(len(heatcapacity)):
        if float(heatcapacity[index][1]) > float(maxheatcapacity):
            maxheatcapacity = heatcapacity[index][1]
            foldingtemperature = heatcapacity[index][0]
    
    return int(float(foldingtemperature))

def findBasins(foldingtemperature):
    unfoldedQ = 0.0
    foldedQ = 0.0
    minfreeenergy = 100000
    for Q in floatRange(0.0,0.5,0.01):
        freeenergy = FofQandT(Q,foldingtemperature)
        if freeenergy < minfreeenergy:
            minfreeenergy = freeenergy
            unfoldedQ = Q

    minfreeenergy = 100000
    for Q in floatRange(0.5,1.0,0.01):
        freeenergy = FofQandT(Q,foldingtemperature)
        if freeenergy < minfreeenergy:
            minfreeenergy = freeenergy
            foldedQ = Q

    return unfoldedQ, foldedQ

def floatRange(a, b, inc):
    try: x = [float(a)]
    except: return False
    for i in range(1, int(math.ceil((b - a ) / inc))):
        x. append(a + i * inc)
    
    return x

def findStabilityGap(unfoldedQ,foldedQ,temperature):
    return FofQandT(foldedQ,temperature) - FofQandT(unfoldedQ,temperature)

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
    for microstate in microstatecodes:
        if calculateRank(microstate) == minrank:
            if microstatefreeenergies[microstateindex] < minrankminfreeenergy:
                minrankminfreeenergy = microstatefreeenergies[microstateindex]

        if calculateRank(microstate) == maxrank:
            if microstatefreeenergies[microstateindex] < maxrankminfreeenergy:
                maxrankminfreeenergy = microstatefreeenergies[microstateindex]

        microstateindex += 1

    return maxrankminfreeenergy - minrankminfreeenergy

def assignAllProbabilities(extrapolatedtemperature):
    for i in range(len(sortedsnapshots)):
        for j in range(len(sortedsnapshots[i])):
            sortedsnapshots[i][j].assignProbability(extrapolatedtemperature)

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

#########
# Files #
#########
# The metadata file, containing links to the dump files and Qw/Potential energy
# files. The format is: dumpfile qw-pot-file
metadataFile = './metadatashort'
# Foldon file: each line contains the residues in a foldon
# each residue in the protein should be included once and only once
foldonFile = './foldons'
# The dump file (LAMMPS format) of the native structure coordinates
nativeDumpFile = './dump.native'
# The directory containing free energy files in the output format of UltimateWHAM
freeEnergyFileDirectory = '/opt/home/ns24/ultimatewham/results/1n0rp1.5-250-400-qw/freeenergyfiles/'
# Overall rate file
overallRateFile = './overallrates'
# Microstate codes file
microstatecodesfile = './microstatecodes'
# Microstate probabilities file
microstateprobabilitiesfile = './microstateprobabilities'
# Microstate free energies file
microstatefreeenergiesfile = './microstatefreeenergies'
# Partition sum file
partitionsumfile = './partitionsum'
# Rate matrix file
ratematrixfile = './ratematrix'
# Eigenvalues and eigenvectors file
eigenvectorsfile = './eigenvectors'
# Trajectories pickle file
trajectoriespicklefile = './trajectories.pkl'
# Heat capapcity file
heatcapacityfile = '/opt/home/ns24/ultimatewham/results/1n0rp1.5-250-400-qw/cv'
# Microstate ranks file prefix
microstateInfoFilePrefix = './microstateinfo'
# Native distance file
nativedistancefile = './rnative.dat'
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
starttemp = 319
# Maximum temperature for computing overall rate
endtemp = 319
# read trajectories from metadata? if not, load trajectories.pkl
readTrajectoriesFromMetadata = False
# output microstate information for each temperature?
outputMicrostateInfo = True
# The frequency at which to accept snapshots from the dump file
snapshotFreq = 1
# Always output rate matrix?
alwaysOutputRateMatrix = True

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
ratematrix = []   # the rate matrix
eigenvectors = [] # eigenvectors of the rate matrix
eigenvalues = []  # eigenvalues of the rate matrix
overallrates = [] # stores all temperature/overall rate pairs calculated
nativecontactlist = [] # stores all native contact pairs
foldons = [] # a list of all foldons
heatcapacity = [] # heat capacity array
foldingTemperature = 0.0 # folding temperature
unfoldedQ = 0.0 # Q of the unfolded basin
foldedQ = 0.0   # Q of the folded basin
stabilitygap = 0.0 # stability gap between unfolded and folded basins
ustateinfo = [] # a list of microstate information for each temperature
nativedistances = [] # a list of the CA-CA distances
minrank = 0 # minimum rank of a sampled microstate
maxrank = 0 # maximum rank of a samples microstate
ustatestabilitygap = 0.0 # free energy gap between microstate with minimum rank, minimum
                         # free energy and maximum rank, minimum free energy
averageqofmicrostates = [] # the average Q of a given microstate

################
# Main program #
################
# read heat capacity file
heatcapacity = readHeatCapacityFile(heatcapacityfile)
# read in all the free energy information, return minimum temperature (for indexing purposes)
print "Reading free energy files..."
mintemp = readAllFreeEnergies(freeEnergyFileDirectory)
# find folding temperature
foldingtemperature = findFoldingTemperature(heatcapacity)
print "Folding temperature: " + str(foldingtemperature)
# find folded and unfolded basins
unfoldedQ, foldedQ = findBasins(foldingtemperature)
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
else:
    # load trajectories from existing Pickle file
    print "Loading pickled trajectories..."
    trajectories = cPickle.load(open(trajectoriespicklefile, 'rb'))

# find all the microstates present
print "Finding all microstates..."
microstatecodes = findAllUstates()
print "Microstate codes: " + str(microstatecodes)
# sort all snapshots according to their microstate
print "Sorting all snapshots by microstate..."
sortedsnapshots = sortAllSnapshots()
# Calculate average global Q for each microstate
print "Calculating average global Q of microstates..."
averageqofmicrostates = calculateAverageQofMicrostates()
print "Average global Q of microstates: " + str(averageqofmicrostates)

# loop over all desired (extrapolated) temperatures, compute overall rate
for temperature in range(starttemp,endtemp+1):
    print "Microstate codes: " + str(microstatecodes)
    print "Average global Q of microstates: " + str(averageqofmicrostates)
    print "Calculating rate for temperature: " + str(temperature)
    # find stability gap
    print "Finding stability gap..."
    stabilitygap = findStabilityGap(unfoldedQ,foldedQ,temperature)
    # calculate partition function for a certain temperature
    print "Calculating the partition function..."
    Z = calcZ(temperature)
    print "Z: " + str(Z)
    # assign probabilities to all snapshots in all trajectories for a certain temperature
    print "Assigning snapshot probabilities..."
    assignAllProbabilities(temperature)
    # calculate microstate probabilities (the temperature is implied by the snapshot probabilities)
    print "Calculating microstate probabilities..."
    microstateprobabilities = calculateUstateProbabilities()
    print "Microstate probabilities: " + str(microstateprobabilities)
    # calculate microstate free energies for a certain temperature
    print "Calculating microstate free energies..."
    microstatefreeenergies = calculateUstateFreeEnergies(temperature)
    print "Microstate free energies: " + str(microstatefreeenergies)
    # calculate rate matrix, given a temperature
    print "Calculating rate matrix..."
    ratematrix = calculateRateMatrix(temperature)
    if(alwaysOutputRateMatrix):
        numpy.savetxt('ratematrix'+str(temperature),ratematrix)
    print "Rate Matrix: \n" + str(ratematrix)
    # calculate and append microstate info
    info = []
    for i in range(len(microstatecodes)):
        info.append([microstatecodes[i],calculateRank(microstatecodes[i]),microstatefreeenergies[i],averageqofmicrostates[i],len(sortedsnapshots[i])])
    ustateinfo.append(info)
    minrank = findMinRank(info)
    maxrank = findMaxRank(info)
    ustatestabilitygap = findUstateStabilityGap(minrank,maxrank)
    # calculate eigenvalues and eigenvectors
    print "Calculating eigenvalues and eigenvectors"
    eigenvalues, eigenvectors = calculateEigenvectorsEigenvalues()
    print "Eigenvalues: \n" + str(eigenvalues)
    print "Eigenvectors: \n " + str(eigenvectors)
    overallrates.append([temperature,stabilitygap,ustatestabilitygap,-eigenvalues[1],-eigenvalues[2]])

f = open(overallRateFile, 'w')
for i in range(len(overallrates)):
    f.write("%f %f %f %e %e \n" % (overallrates[i][0],overallrates[i][1],overallrates[i][2],overallrates[i][3],overallrates[i][4]))
f.close()

if outputMicrostateInfo:
    for temperature in range(starttemp,endtemp+1):
        f = open(microstateInfoFilePrefix + "." + str(temperature), 'w')
        ustateinformation = ustateinfo[temperature-starttemp]
        for index in range(len(ustateinformation)):
            f.write("%s %d %f %f %d \n" % (ustateinformation[index][0], ustateinformation[index][1], ustateinformation[index][2], ustateinformation[index][3], ustateinformation[index][4]))
        f.close()

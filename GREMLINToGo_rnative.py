import os.path
import argparse

import pandas as pd
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("contactFile", help="direct path to GREMLIN contact prediction file")
parser.add_argument("sequenceFile", help="direct path to AWSEM sequence (.seq) file")
parser.add_argument("pThreshold", help="probability threshold (generally 0.5)")
args = parser.parse_args()

filter_threshold = args.pThreshold
GREMLIN_file = args.contactFile
sequence_file = args.sequenceFile

#read in median distances for pairwise interactions (obtained from bioinformatic analysis of the pdb)
distancesCACB=pd.read_table('CACBmediandist.dat', delim_whitespace=True, header=None)
distancesCACA=pd.read_table('CACAmediandist.dat', delim_whitespace=True, header=None)
distancesCBCB=pd.read_table('CBCBmediandist.dat', delim_whitespace=True, header=None)
distancesCACB.columns = ['i', 'j', 'dist']
distancesCACA.columns = ['i', 'j', 'dist']
distancesCBCB.columns = ['i', 'j', 'dist']

#if you want to filter the contact prediction data, adjust the parameters below
probability_column=5

#make sure that there is a sequence file for the protein and the downloaded gremlin data in the proper directories
f = open(sequence_file, "r")
n=len(f.readlines()[0].strip())
f.close()

#load downloaded gremlin file
gremlin_data=np.loadtxt(GREMLIN_file, dtype=str, skiprows=1)
rnative_matrixCACB=np.ones([n,n])*99
rnative_matrixCACA=np.ones([n,n])*99
rnative_matrixCBCB=np.ones([n,n])*99

for pair in gremlin_data:
    i=int(pair[0])
    j=int(pair[1])
    irestype=pair[2][-1:]
    jrestype=pair[3][-1:]
    if float(pair[probability_column]) > float(filter_threshold):
        if sum((distancesCACB['i']==irestype)&(distancesCACB['j']==jrestype))>0: #check if pair is in correct order
            well_centerCACB = distancesCACB[(distancesCACB['i']==irestype)&(distancesCACB['j']==jrestype)]['dist'].values[0]
            well_centerCACA = distancesCACA[(distancesCACA['i']==irestype)&(distancesCACA['j']==jrestype)]['dist'].values[0]
            well_centerCBCB = distancesCBCB[(distancesCBCB['i']==irestype)&(distancesCBCB['j']==jrestype)]['dist'].values[0]            
        else:
            well_centerCACB = distancesCACB[(distancesCACB['i']==jrestype)&(distancesCACB['j']==irestype)]['dist'].values[0]
            well_centerCACA = distancesCACA[(distancesCACA['i']==jrestype)&(distancesCACA['j']==irestype)]['dist'].values[0]
            well_centerCBCB = distancesCBCB[(distancesCBCB['i']==jrestype)&(distancesCBCB['j']==irestype)]['dist'].values[0]

        rnative_matrixCACB[i-1, j-1] = well_centerCACB
        rnative_matrixCACB[j-1, i-1] = well_centerCACB
        rnative_matrixCACA[i-1, j-1] = well_centerCACA
        rnative_matrixCACA[j-1, i-1] = well_centerCACA
        rnative_matrixCBCB[i-1, j-1] = well_centerCBCB
        rnative_matrixCBCB[j-1, i-1] = well_centerCBCB

np.savetxt('go_rnativeCACB.dat', rnative_matrixCACB, fmt='%10.5f')
np.savetxt('go_rnativeCACA.dat', rnative_matrixCACA, fmt='%10.5f')
np.savetxt('go_rnativeCBCB.dat', rnative_matrixCBCB, fmt='%10.5f')

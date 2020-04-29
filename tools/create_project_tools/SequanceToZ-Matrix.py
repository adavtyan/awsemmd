#!/usr/bin/python

# ----------------------------------------------------------------------
# Copyright (2010) Aram Davtyan and Garegin Papoian

# Papoian's Group, University of Maryland at Collage Park
# http://papoian.chem.umd.edu/

# Last Update: 03/04/2011
# ----------------------------------------------------------------------
import os
import sys

class ZAtom:
    def __init__(self, no=0, ch=1, ty='', desc='', d=0, d_no=0, a=0, a_no=0, b=0, b_no=0):
        self.no = no
        self.ch = ch
        self.ty = ty
        self.d = d
        self.d_no = d_no
        self.a = a
        self.a_no = a_no
        self.b = b
        self.b_no = b_no
        self.desc = desc

    def print_(self, wdesc=False):
        pr = [ str(self.no), str(self.ch), self.ty ]
        if wdesc: pr += [ (self.desc+"            ")[0:12] ]
        if self.d_no>0: pr += [ str(round(self.d,3)), str(self.d_no) ]
        if self.a_no>0: pr += [ str(round(self.a,1)), str(self.a_no) ]
        if self.b_no>0: pr += [ str(round(self.b,1)), str(self.b_no) ]
        print '\t'.join(pr)

    def write_(self, f, wdesc=False):
        pr = [ str(self.no), str(self.ch), self.ty ]
        if wdesc: pr += [ (self.desc+"            ")[0:12] ]
        if self.d_no>0: pr += [ str(round(self.d,3)), str(self.d_no) ]
        if self.a_no>0: pr += [ str(round(self.a,1)), str(self.a_no) ]
        if self.b_no>0: pr += [ str(round(self.b,1)), str(self.b_no) ]
        pr += ['\n']
        f.write('\t'.join(pr))

if len(sys.argv)==1:
    print "\nSequanceToZ-Matrix.py input_file [output_file] [-d]\n"
    print "\t-d\tInclude descripton string"
    exit()

incDesc = False
for arg in sys.argv:
    if arg=="-d":
        incDesc = True
        sys.argv.remove(arg)

inpfname = ""
outfname = ""
if len(sys.argv)>1: inpfname = sys.argv[1]
if len(sys.argv)>2: outfname = sys.argv[2]

myhome = os.environ.get("HOME")
fdata = open(myhome + "/opt/script/SequanceToZ-Matrix.data")
raw_data = []
iFirst = 0
iSecond = 0
iFirst_GLY = 0
iSecond_GLY = 0
Error = True
i = 0
for l in fdata:
    if l.strip()=="": break
    rd = l.strip().split()
    if rd[0][0]=='[' and rd[0][-1]==']':
        if rd[0]=='[First_Residue_Chain]':
            iFirst = i+1
        elif rd[0]=='[Second_Residue_Chain]':
            iSecond = i+1
        elif rd[0]=='[First_GLY_Residue_Chain]':
            iFirst_GLY = i+1
        elif rd[0]=='[Second_GLY_Residue_Chain]':
            iSecond_GLY = i+1
    
    raw_data.append(rd)
    i += 1
fdata.close()


inp = open(inpfname)

seq = inp.readline().rstrip()

nAtoms = 0
ZMatrix = []
nAltered = [0, 0]
standartPhi = 60.0
standartPsi = 180.0
phi = []
psi = []
for i in range(0,len(seq)):
    if seq[i]!='G': iChainStart = iFirst
    else: iChainStart = iFirst_GLY
    if i>0:
        if seq[i]!='G': iChainStart = iSecond
        else: iChainStart = iSecond_GLY

    noFirstInRes = nAtoms + 1
    nAltered[0] = nAltered[1]
    nAltered[1] = 0
    
    angles = inp.readline().split()
    phi = [float(angles[0])] + phi
    psi = [float(angles[1])] + psi

    for jrd in raw_data[iChainStart:]:
        if jrd[0][0]=='[': break

        if jrd[2]=='O-In-The-End' and i!=len(seq)-1: continue
#        if jrd[2]=='C-Beta' and seq[i]=='G':
#            nAltered[1] += 1
#            continue
        
        if len(jrd)==3:
            atom = ZAtom(int(jrd[0]), 1, jrd[1], jrd[2])
        elif len(jrd)==5:
            atom = ZAtom(int(jrd[0]), 1, jrd[1], jrd[2], float(jrd[3]), int(jrd[4]))
        elif len(jrd)==7:
            atom = ZAtom(int(jrd[0]), 1, jrd[1], jrd[2], float(jrd[3]), int(jrd[4]), float(jrd[5]), int(jrd[6]))
        elif len(jrd)==9:
            atom = ZAtom(int(jrd[0]), 1, jrd[1], jrd[2], float(jrd[3]), int(jrd[4]), float(jrd[5]), int(jrd[6]), float(jrd[7]), int(jrd[8]))

        dn = nAtoms - atom.no + 1
        atom.d_no += dn
        if atom.d_no < noFirstInRes: atom.d_no += nAltered[0]
        atom.a_no += dn
        if atom.a_no < noFirstInRes: atom.a_no += nAltered[0]
        atom.b_no += dn
        if atom.b_no < noFirstInRes: atom.b_no += nAltered[0]
        atom.no = nAtoms + 1

        if jrd[2]=='C-Prime':
            atom.b +=  phi[0] - standartPhi
        elif jrd[2]=='N':
            atom.b +=  psi[1] - standartPsi
        elif jrd[2]=='O' or jrd[2]=='O-In-The-End' or jrd[2]=='C-Beta' or jrd[2]=='H-Beta':
            atom.b +=  psi[0] - standartPsi

        ZMatrix.append(atom)
        nAtoms += 1
inp.close()

if outfname=="":
    for iZM in ZMatrix:
        iZM.print_(incDesc)
else:
    out = open(outfname, 'w')
    for iZM in ZMatrix:
        iZM.write_(out, incDesc)
    out.close()


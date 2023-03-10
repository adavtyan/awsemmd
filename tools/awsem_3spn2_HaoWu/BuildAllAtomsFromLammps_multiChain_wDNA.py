#!/usr/bin/python

# ----------------------------------------------------------------------
# Copyright (2010) Aram Davtyan and Garegin Papoian

# Papoian's Group, University of Maryland at Collage Park
# http://papoian.chem.umd.edu/

# Last Update: 03/04/2011

# Fix problem for single snapshot version
# Reformat for correct indentation & clear comments
# Hao Wu Apr 06 2018
# ----------------------------------------------------------------------

import sys
from numpy import *

#from Bio.PDB.PDBParser import PDBParser

#----------------------------------------------------------------------
# Predefined constants, tables and dictionaries
#----------------------------------------------------------------------

atom_type = {'1' : 'C', '2' : 'N', '3' : 'O', '4' : 'C', '5' : 'H', '6' : 'C', \
             '15' : 'C', '16' : 'N', '17' : 'O', '18' : 'C', '19' : 'H', '20' : 'C' \
            }
atom_desc = {'1' : 'C-Alpha', '2' : 'N', '3' : 'O', '4' : 'C-Beta', '5' : 'H-Beta', '6' : 'C-Prime',        \
             '15' : 'C-Alpha', '16' : 'N', '17' : 'O', '18' : 'C-Beta', '19' : 'H-Beta', '20' : 'C-Prime'   \
            }
PDB_type = {'1' : 'CA', '2' : 'N', '3' : 'O', '4' : 'CB', '5' : 'HB', '6' : 'C',        \
            '15' : 'CA', '16' : 'N', '17' : 'O', '18' : 'CB', '19' : 'HB', '20' : 'C' }

d_res = {"C" : "CYS", "I" : "ILE", "S" : "SER", "Q" : "GLN", "K" : "LYS",
         "N" : "ASN", "P" : "PRO", "T" : "THR", "F" : "PHE", "A" : "ALA",
         "H" : "HIS", "G" : "GLY", "D" : "ASP", "L" : "LEU", "R" : "ARG",
         "W" : "TRP", "V" : "VAL", "E" : "GLU", "Y" : "TYR", "M" : "MET"}

chainId2Letter = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', \
                  'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', \
                  'Y', 'Z']

an = 0.4831806
bn = 0.7032820
cn = -0.1864262
ap = 0.4436538
bp = 0.2352006
cp = 0.3211455

#----------------------------------------------------------------------
# Class: PDB_Atom
#----------------------------------------------------------------------

class PDB_Atom:
    no = 0
    ty = ''
    res = 'UNK'
    res_no = 0
    x = 0.0
    y = 0.0
    z = 0.0
    atm = 'C'

    def __init__(self, no, ty, res, res_no, x, y, z, cid, atm):
        self.no = no
        self.ty = ty
        self.res = res
        self.res_no = res_no
        self.x = x
        self.y = y
        self.z = z
        self.cid = cid
        self.atm = atm

    def write_(self, f):
        f.write('ATOM')
        f.write(('       ' + str(self.no))[-7:])
        f.write('  ')
        f.write((self.ty + '    ')[:4])
        f.write(self.res)
        f.write(' ')
        f.write(chainId2Letter[self.cid])
        f.write(('    '+str(self.res_no))[-4:])
        f.write(('            '+str(round(self.x,3)))[-12:])
        f.write(('        '+str(round(self.y,3)))[-8:])
        f.write(('        '+str(round(self.z,3)))[-8:])
        f.write('  1.00')
        f.write('  0.00')
        f.write(('            '+self.atm)[-12:]+'  ')
        f.write('\n')

#----------------------------------------------------------------------
# Class: Atom
#----------------------------------------------------------------------

class Atom:
    No = 0
    ty = ''
    x = 0.0
    y = 0.0
    z = 0.0
    desc = ''

    def __init__(self, No, ty, No_m, x, y, z, chainId, desc=''):
        self.No = No
        self.ty = ty
        self.No_m = No_m
        self.x = x
        self.y = y
        self.z = z
        self.cid = chainId
        self.desc = desc

    def write_(self, f):
        f.write(str(self.No))
        f.write(' ')
        f.write(PDB_type[self.No_m])
        f.write(' ')
        f.write(str(round(self.x,8)))
        f.write(' ')
        f.write(str(round(self.y,8)))
        f.write(' ')
        f.write(str(round(self.z,8)))
        f.write(' ')
        f.write(self.desc)
        f.write('\n')

def One2ThreeLetters(txt):
    seq_array = []
    for s in txt:
        seq_array.append( d_res.get(s,"ALA") )
    return seq_array

#----------------------------------------------------------------------
# Read input variables (including DNA)
#----------------------------------------------------------------------

if len(sys.argv) < 3 or len(sys.argv) > 10:
    print ("\n" + sys.argv[0] + "Input_file Output_file [snapshot] [-seq sequence_file] [-dna dnaPDBTemplateFile] [-dna dnaBondFile]\n")
    exit()

del_list = []
seq_file = ""
dnaPDBFIle  = ""
dnaBondFIle = ""
buildDNA = False
sequence = []
seq_txt = ""
chainLen = []
bseq = False

for iarg in range(3, len(sys.argv)):
    if sys.argv[iarg] == "-seq":
        bseq = True
        seq_file = sys.argv[iarg+1]
        del_list.insert(0, iarg)
        del_list.insert(0, iarg+1)
    if sys.argv[iarg] == "-dnaPdb":
        buildDNA = True
        dnaPDBFile = sys.argv[iarg+1]
        del_list.insert(0, iarg)
        del_list.insert(0, iarg+1)
    if sys.argv[iarg] == "-dnaBond":
        dnaBondFile = sys.argv[iarg+1]
        del_list.insert(0, iarg)
        del_list.insert(0, iarg+1)
for idel in del_list:
    sys.argv.pop(idel)

lammps_file = sys.argv[1]

output_file = ""
if len(sys.argv) > 2: output_file = sys.argv[2]
psf_file = output_file

if output_file[-4:] != ".pdb": output_file = output_file + ".pdb"
if psf_file[-4:] == ".pdb": psf_file = psf_file[:-3] + "psf"
if psf_file[-4:] != ".psf": psf_file = psf_file + ".psf"

snapshot = -1
if len(sys.argv) > 3: snapshot = int(sys.argv[3])

# protein sequence
if seq_file != "":
    seq_txt = ''
    fh = open(seq_file, 'r')
    for line in fh.readlines():
        seq = line.strip()
        chainLen.append(len(seq))
        seq_txt=seq_txt+seq
    fh.close()
    sequence = One2ThreeLetters(seq_txt)

#----------------------------------------------------------------------
# Read DNA atoms from PDB and bond
#----------------------------------------------------------------------

# ---               build dna Atoms             --- #
dnaAtomList = []
dnaBondList = []
if buildDNA is True:

    # read DNA PDB file
    fh = open(dnaPDBFile, 'r')
    for line in fh.readlines()[1:-1]:
        # items   = line.split()
        # No      = int(items[1])
        # type    = items[2]
        # res_ty  = items[3]
        # chainId = chainId2Letter.index(items[4])
        # ires    = int(items[5])
        # x       = float(items[6])
        # y       = float(items[7])
        # z       = float(items[8])
        # atomType = items[-1]

        # Use PDB format to split each item instead of splitting by whitespaces
        # based on http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html
        # avoid long coordinates
        # Hao Wu Apr 05 2018
        No = int(line[6:11].replace(" ", ""))
        type = line[12:16].replace(" ", "")
        res_ty = line[17:20].replace(" ", "")
        chainId = chainId2Letter.index(line[21])
        ires = int(line[22:26].replace(" ", ""))
        x = float(line[30:38].replace(" ", ""))
        y = float(line[38:46].replace(" ", ""))
        z = float(line[46:54].replace(" ", ""))
        atomType = line[77]

        atom = PDB_Atom(No, type, res_ty, ires, x, y, z, chainId, atomType)
        dnaAtomList.append(atom)
    fh.close()

    # read dna bonds
    fh = open(dnaBondFile, 'r')
    for line in fh.readlines()[1:-1]:
        items   = line.split()
        if len(items) > 0:
            dnaBondList.append([int(items[0]), int(items[1])])
            if len(items) > 2:
                dnaBondList.append([int(items[2]), int(items[3])])
                if len(items) > 4:
                    dnaBondList.append([int(items[4]), int(items[5])])
                    if len(items) > 6:
                        dnaBondList.append([int(items[6]), int(items[7])])

    fh.close()

n_atoms = 0
i_atom = 0
item = ''
step = 0
atoms = []
atoms2 = []
atoms3 = []
bonds = []
box = []
A = []

out = open(output_file, 'w')

#----------------------------------------------------------------------
# calculate chain ID
#----------------------------------------------------------------------
atomInd2ChainId = []
def calcChainId(n_atoms):
    chLenList = array(chainLen)
    chLenListCS = cumsum(chLenList)
    for atomInd in range(n_atoms):
        resInd = floor(atomInd/3) + 1
        chainId = searchsorted(chLenListCS, resInd)
        atomInd2ChainId.append(chainId)

#----------------------------------------------------------------------
# convert protein atoms2 to atoms3 based on PDB_Atom class format
#----------------------------------------------------------------------
def convertToPDB():
    ires = 1
    ires4seq = 1
    # newChain = False
    for i, ia in enumerate(atoms2):
        # new chain
        if (i > 0) and (ia.cid != atoms2[i-1].cid):
            ires=1
            ires4seq += 1
        # new residue (N-terminal)
        if ia.desc == 'N':
            ires = ires + 1
            ires4seq +=  1
        # residue type based on input sequence, "ALA" if not specified
        res_ty = "ALA"
        if bseq: res_ty = sequence[ires4seq-1]
        # write atom
        atom = PDB_Atom(ia.No, PDB_type[ia.No_m], res_ty, ires, ia.x, ia.y, ia.z, ia.cid, ia.ty)
        atoms3.append(atom)

#----------------------------------------------------------------------
# Build protein atoms2Tmp based on atomsTmp
#----------------------------------------------------------------------
def buildAllAtoms(atomsTmp):

    atoms2Tmp = []
    index = 0
    last_Ca_index = -1
    last_O_index = -1
    Cp_index = -1
    NullVal = Atom(0, '', '6', 0.0, 0.0, 0.0, '')

    for i in range(0, len(atomsTmp)):
        ia = atomsTmp[i]
        index = index + 1
        if ia.desc == 'O': last_O_index = i
        if ia.desc == 'C-Alpha':
            if last_Ca_index != -1:
                Cai = atomsTmp[last_Ca_index]
                Cai1 = ia
                Oi = atomsTmp[last_O_index]

                nx = an*Cai.x + bn*Cai1.x + cn*Oi.x
                ny = an*Cai.y + bn*Cai1.y + cn*Oi.y
                nz = an*Cai.z + bn*Cai1.z + cn*Oi.z

                px = ap*Cai.x + bp*Cai1.x + cp*Oi.x
                py = ap*Cai.y + bp*Cai1.y + cp*Oi.y
                pz = ap*Cai.z + bp*Cai1.z + cp*Oi.z

                N = Atom(index, 'N', '2', nx, ny, nz, ia.cid, 'N')
                index = index + 1
                Cp = Atom(int(Cai.No) + 1, 'C', '6', px, py, pz, ia.cid, 'C-Prime')

                atoms2Tmp.append(N)
                atoms2Tmp.pop(Cp_index)
                atoms2Tmp.insert(Cp_index, Cp)

            last_Ca_index = i
        ia.No = index
        atoms2Tmp.append(ia)
        if ia.desc == 'C-Alpha':
            atoms2Tmp.append(NullVal)
            Cp_index = index
            index = index + 1
    if atoms2Tmp[Cp_index].No == 0: atoms2Tmp.pop(Cp_index)
    for i in range(Cp_index, len(atoms2Tmp)):
        atoms2Tmp[i].No = atoms2Tmp[i].No - 1

    return atoms2Tmp

#----------------------------------------------------------------------
# Build protein bonds
#----------------------------------------------------------------------
def buildBonds(atoms2Tmp):
    N_index = -1
    Ca_index = -1
    Cp_index = -1
    O_index = -1
    Cb_index = -1
    Hb_index = -1
    bondsTmp = [];
    startingIndex = atoms2Tmp[0].No

    for i in range(0, len(atoms2Tmp)):
        ia = atoms2Tmp[i]
        if ia.desc == 'N':
            if N_index!=-1 and Ca_index!=-1:
                bondsTmp.append([N_index, Ca_index])
            if Ca_index!=-1 and Cp_index!=-1:
                bondsTmp.append([Ca_index, Cp_index])
            if Cp_index!=-1 and O_index!=-1:
                bondsTmp.append([Cp_index, O_index])
            if Ca_index!=-1 and Cb_index!=-1:
                bondsTmp.append([Ca_index, Cb_index])
            if Ca_index!=-1 and Hb_index!=-1:
                bondsTmp.append([Ca_index, Hb_index])
            N_index = i + startingIndex
            if Cp_index!=-1:
                bondsTmp.append([Cp_index, N_index])
            Ca_index = -1
            Cp_index = -1
            O_index = -1
            Cb_index = -1
            Hb_index = -1
        if ia.desc == 'C-Alpha': Ca_index = i + startingIndex
        if ia.desc == 'C-Beta': Cb_index = i + startingIndex
        if ia.desc == 'H-Beta': Hb_index = i + startingIndex
        if ia.desc == 'C-Prime': Cp_index = i + startingIndex
        if ia.desc == 'O': O_index = i + startingIndex
    if N_index!=-1 and Ca_index!=-1:
        bondsTmp.append([N_index, Ca_index])
    if Ca_index!=-1 and Cb_index!=-1:
        bondsTmp.append([Ca_index, Cb_index])
    if Ca_index!=-1 and Hb_index!=-1:
        bondsTmp.append([Ca_index, Hb_index])

    return bondsTmp

#----------------------------------------------------------------------
# Print atom array from LAMMPS dump file
#----------------------------------------------------------------------
def print_atom_array():
    out.write("ITEM: TIMESTEP\n")
    out.write(str(step))
    out.write("\n")
    out.write("ITEM: NUMBER OF ATOMS\n")
    out.write(str(n_atoms))
    out.write("\n")
    out.write("ITEM: BOX BOUNDS\n")
    for ib in box:
        out.write(ib)
        out.write("\n")
    out.write("ITEM: ATOMS\n")
    for ia in atoms2:
        ia.write_(out)

#----------------------------------------------------------------------
# print (PDB file with DNA)
#----------------------------------------------------------------------
def print_pdb():
    for ia in atoms3:
        ia.write_(out)
    if buildDNA is True:
        for ia in dnaAtomList:
            ia.write_(out)
    out.write("END\n");

#----------------------------------------------------------------------
# print (PSF file)
#----------------------------------------------------------------------
def print_psf():
    space8 = "        "
    psfout = open(psf_file,'w')
    psfout.write("PDF\n\n\t2 !NTITLE\n\n")
    if buildDNA is True:
        psfout.write((space8+str(len(atoms3)+len(dnaAtomList)))[-8:]+" !NATOM\n")
    else:
        psfout.write((space8+str(len(atoms3)))[-8:]+" !NATOM\n")
    for ia in atoms2:
        psfout.write((space8+str(ia.No))[-8:]+" PROT 1")
        psfout.write("    R00")
        psfout.write("  "+ia.ty)
        psfout.write("       1")
        psfout.write("          0            1           0\n")

    if buildDNA is True:
        for ia in dnaAtomList:
            psfout.write((space8+str(ia.no))[-8:]+" DNA  1")
            psfout.write("    R00")
            psfout.write("  "+ia.ty)
            psfout.write("       1")
            psfout.write("          0            1           0\n")
    psfout.write("\n")

    if buildDNA is True:
        # include protein atoms in the index
        for i in range(0, len(dnaBondList)):
            ib = dnaBondList[i]
            ib[0] += len(atoms2)
            ib[1] += len(atoms2)
            bonds.append(ib)
    psfout.write((space8+str(len(bonds)))[-8:]+" !NBOND")

    for i in range(0, len(bonds)):
        ib = bonds[i]
        if i%4==0: psfout.write("\n")
        psfout.write((space8+str(ib[0]))[-8:])
        psfout.write((space8+str(ib[1]))[-8:])

    psfout.close()

#----------------------------------------------------------------------
# main process
#----------------------------------------------------------------------
nFrame = 0
found = False
lfile = open(lammps_file)

setChainIdList = False
setDNAAtomIndex = False

#----------------------------------------------------------------------
# snapshot not selected: entire trajectory
#----------------------------------------------------------------------
if snapshot < 0:
    for l in lfile:
        l = l.strip()
        if l[:5] == "ITEM:":
            item = l[6:]
        else:
            if item == "TIMESTEP":
                if len(atoms) > 0:
                    ista = 0; iend = 0;
                    for ic in range(len(chainLen)):
                        ista = iend
                        iend = ista + chainLen[ic]*3
                        atoms2Tmp = buildAllAtoms(atoms[ista:iend])
                        startingIndex = len(atoms2)
                        for ia in atoms2Tmp:
                            ia.No += startingIndex
                        atoms2 += atoms2Tmp
                    convertToPDB()
                    n_atoms = len(atoms2)

                    if buildDNA is True:
                        for i, ia in enumerate(atoms[iend::]):
                            if setDNAAtomIndex is False:
                                dnaAtomList[i].no += n_atoms
                            dnaAtomList[i].x = ia.x
                            dnaAtomList[i].y = ia.y
                            dnaAtomList[i].z = ia.z
                        setDNAAtomIndex = True

                    print_pdb()
                step = int(l)
                atoms = []
                atoms2 = []
                atoms3 = []
                box = []
                A = []
                nFrame = nFrame + 1
            elif item == "NUMBER OF ATOMS":
                n_atoms = int(l)
                if setChainIdList is not True:
                    calcChainId(n_atoms)
                    setChainIdList = True
            elif item[:10] == "BOX BOUNDS":
                box.append(l)
                l = l.split()
                A.append([float(l[0]), float(l[1])])
            elif item[:5] == "ATOMS":
                l = l.split()
                i_atom = l[0]
                chainId = atomInd2ChainId[int(i_atom)-1]
                x = float(l[2])
                y = float(l[3])
                z = float(l[4])
                x = (A[0][1] - A[0][0])*x + A[0][0]
                y = (A[1][1] - A[1][0])*y + A[1][0]
                z = (A[2][1] - A[2][0])*z + A[2][0]
                if (buildDNA is True) and (int(l[1]) < 15):
                    desc = "DNA"
                    atom = Atom(i_atom, 'P', l[1], x, y, z, chainId, desc)
                else:
                    desc = atom_desc[l[1]]
                    atom = Atom(i_atom, atom_type[l[1]], l[1], x, y, z, chainId, desc)
                atoms.append(atom)

    if len(atoms) > 0:
        ista = 0; iend = 0;
        for ic in range(len(chainLen)):
            ista = iend
            iend = ista + chainLen[ic]*3
            atoms2Tmp = buildAllAtoms(atoms[ista:iend])
            startingIndex = len(atoms2)
            for ia in atoms2Tmp:
                ia.No += startingIndex
            atoms2 += atoms2Tmp
            bonds += buildBonds(atoms2Tmp)
        convertToPDB()
        n_atoms = len(atoms2)

        if buildDNA is True:
            for i, ia in enumerate(atoms[iend::]):
                if setDNAAtomIndex is False:
                    dnaAtomList[i].no += n_atoms
                dnaAtomList[i].x = ia.x
                dnaAtomList[i].y = ia.y
                dnaAtomList[i].z = ia.z
            setDNAAtomIndex = True

        print_pdb()
        print_psf()

#----------------------------------------------------------------------
# snapshot selected: single snapshot
#----------------------------------------------------------------------
else:
    for l in lfile:
        l = l.strip()
        if l[:5] == "ITEM:":
            item = l[6:]
            if item == "TIMESTEP":
                if found: break
                elif nFrame == snapshot: found = True
                nFrame = nFrame + 1
        elif found:
            if item == "TIMESTEP":
                if len(atoms) > 0:
                    ista = 0; iend = 0;
                    for ic in range(len(chainLen)):
                        ista = iend
                        iend = ista + chainLen[ic]*3
                        atoms2Tmp = buildAllAtoms(atoms[ista:iend])
                        startingIndex = len(atoms2)
                        for ia in atoms2Tmp:
                            ia.No += startingIndex
                        atoms2 += atoms2Tmp
                    #buildAllAtoms()
                    convertToPDB()
                    n_atoms = len(atoms2)

                    if buildDNA is True:
                        for i, ia in enumerate(atoms[iend::]):
                            if setDNAAtomIndex is False:
                                dnaAtomList[i].no += n_atoms
                            dnaAtomList[i].x = ia.x
                            dnaAtomList[i].y = ia.y
                            dnaAtomList[i].z = ia.z
                        setDNAAtomIndex = True

                    print_pdb()
                step = int(l)
            elif item == "NUMBER OF ATOMS":
                n_atoms = int(l)
                if setChainIdList is not True:
                    calcChainId(n_atoms)
                    setChainIdList = True
            elif item[:10] == "BOX BOUNDS":
                box.append(l)
                l = l.split()
                A.append([float(l[0]), float(l[1])])
            elif item[:5] == "ATOMS":
                l = l.split()
                i_atom = l[0]
                chainId = atomInd2ChainId[int(i_atom)-1]
                x = float(l[2])
                y = float(l[3])
                z = float(l[4])
                x = (A[0][1] - A[0][0])*x + A[0][0]
                y = (A[1][1] - A[1][0])*y + A[1][0]
                z = (A[2][1] - A[2][0])*z + A[2][0]
                if (buildDNA is True) and (int(l[1]) < 15):
                    desc = "DNA"
                    atom = Atom(i_atom, 'P', l[1], x, y, z, chainId, desc)
                else:
                    desc = atom_desc[l[1]]
                    atom = Atom(i_atom, atom_type[l[1]], l[1], x, y, z, chainId, desc)
                atoms.append(atom)

    if len(atoms) > 0:
        ista = 0; iend = 0;
        for ic in range(len(chainLen)):
            ista = iend
            iend = ista + chainLen[ic]*3
            atoms2Tmp = buildAllAtoms(atoms[ista:iend])
            startingIndex = len(atoms2)
            for ia in atoms2Tmp:
                ia.No += startingIndex
            atoms2 += atoms2Tmp
            bonds += buildBonds(atoms2Tmp)
        convertToPDB()
        n_atoms = len(atoms2)

        if buildDNA is True:
            for i, ia in enumerate(atoms[iend::]):
                if setDNAAtomIndex is False:
                    dnaAtomList[i].no += n_atoms
                dnaAtomList[i].x = ia.x
                dnaAtomList[i].y = ia.y
                dnaAtomList[i].z = ia.z
            setDNAAtomIndex = True

        print_pdb()
        print_psf()

lfile.close()
out.close()

#!/usr/bin/python
#rmbh11
#v0.72
#corects and builds on release version - v0.5

#USAGE: prepFragsLAMW.py database-prefix file.fasta > logfile

import sys,os,re
from Pdb2GroLib import *

if len(sys.argv)!=3:
	print "\n prepFragsLAMW.py database file.fasta > logfile \n"
	exit()

database=sys.argv[1]
inFASTA=sys.argv[2]
weight=1 #feature in match file

from Bio import SeqIO
inseq=SeqIO.read(inFASTA,'fasta')
print "processing: ",inseq.name
query=str(inseq.name)[0:4]

pdbDir = "./PDBs/"
fLibDir = "./fraglib/"

fragmentLength=9 #needs to be an odd number
memoriesPerPosition=10 #can be any integer > 0
EvalueThreshold=10000 #needs to be large enough that PSI-BLAST returns at least memoriesPerPosition
#
##
##SANITY CHECKING
#is length greater than 9 residues?
if(len(inseq.seq) < fragmentLength):
    print "Exception::query sequence is smaller than "+str(fragmentLength)+" residues"
    print "This version has no means to handle smaller queries"
    sys.exit()

##Create necessary directories
if not os.path.exists(pdbDir): os.makedirs(pdbDir)
if not os.path.exists(fLibDir): os.makedirs(fLibDir)

if not os.path.exists(pdbDir) or not os.path.exists(fLibDir):
	print "Can't create necessary directories"
	sys.exit()

##open match file
match=open('prepFrags.match','w')
match.write(query+"\n")
##FRAGMENT GENERATION LOOP
iterations=len(inseq.seq)-fragmentLength+1 #number of sliding window positions
for i in range(1,iterations+1): 
##select subrange
    print "window position:::"+str(i)
    rangeStart=i-1
    rangeEnd=i+fragmentLength-1
    subrange=str(inseq[rangeStart:rangeEnd].seq)
    fragment=open('fragment.fasta','w')
    print "fragment subrange:::"+subrange
    fragment.write(subrange)
    fragment.close()
##submit PSI-BLAST -- run psiblast -help for explanation of format 6 output
#    exeline="psiblast -query_loc "+str(rangeStart)+"-"+str(rangeEnd)
    exeline="psiblast -num_iterations 1 -word_size 2 -evalue "+str(EvalueThreshold)
    exeline+=" -outfmt 6 -matrix BLOSUM62 -db "
    exeline+=database+" -query fragment.fasta"
    print "executing:::"+exeline
    psiblastOut=os.popen(exeline).read()
    psiblastOut=psiblastOut.splitlines() #now an array
#column 7,8,9,10 (starting at 1) are the indices for the aligned ranges
    print "PDB INSEQ-START INSEQ-END MATCH-START MATCH-END EVALUE"
    for line in psiblastOut[0:memoriesPerPosition]:
        this=line.split()
#        print this[1],this[6],this[7],this[8],this[9],this[11]
        queryStart=rangeStart+int(this[6])
        queryEnd=queryStart+int(this[7])-1
        print this[1],str(queryStart),str(queryEnd),this[8],this[9],this[11]
        out=this[1]+' '+str(queryStart)+' '+str(queryEnd)+' '
        out+=this[8]+' '+this[9]+' '+str(weight)+"\n"
        delQuery=queryEnd-queryStart
        delAlign=int(this[9])-int(this[8])
        #if residue ranges do not match, this alignment was gapped
        #skip gapped alignments:                                  ################################
        if ((delQuery-delAlign)==0):
            match.write(out)
##
match.close()
match=open('prepFrags.match','r') #match is read-only now
LAMWmatch=open('fragsLAMW.mem','w')
LAMWmatch.write('[Target]'+"\n")
LAMWmatch.write(query+"\n\n"+'[Memories]'+"\n")
##get pdbs
matchlines=list()
keys = {}
for line in match.readlines():
    matchlines.append(line)
    entries=line.split()
    pdbfull=str(entries[0])
    keys[pdbfull]=1
    unique=keys.keys()
#
from Bio.PDB.PDBParser import PDBParser
pdbparse=PDBParser(PERMISSIVE=1)
atomLine=re.compile('\AATOM')
#from Bio.PDB import PDBIO
for pdbfull in unique:
    pdbID=pdbfull[0:4].lower()
    pdbIDsecond=pdbfull[1:2].lower()
    pdbIDthird=pdbfull[2:3].lower()
    chainID=pdbfull[4:5].lower()
    if not os.path.isfile(pdbDir+pdbID.upper()+".pdb"):
##from script 'pdbget' (original author unknown)
        exeline="wget ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb/"
        exeline+=pdbIDsecond+pdbIDthird+"/pdb"+pdbID+".ent.gz"
        os.system(exeline);
        os.system("nice gunzip pdb"+pdbID+".ent.gz; mv pdb"+pdbID+".ent "+pdbDir+pdbID.upper()+".pdb");

    if not os.path.isfile(pdbDir+pdbID.upper()+".pdb"):
        print ":::Cannot build PDB for PDB ID, failed to download:"+pdbID.upper()
        
iter=0
for line in matchlines:
    iter+=1
    if not(iter==1):
        print ":::here: match line:"+line
        entries=line.split()
        pdbfull=str(entries[0])
        pdbID=pdbfull[0:4].lower()
        pdbIDsecond=pdbfull[1:2].lower()
        pdbIDthird=pdbfull[2:3].lower()
        chainID=pdbfull[4:5].lower()
	groFile=fLibDir+pdbID+chainID+".gro"
	pdbFile=pdbDir+pdbID.upper()+".pdb"
        if os.path.isfile(pdbFile):
            print ":::convert: "+pdbFile+" --> "+groFile
            Pdb2Gro(pdbFile, groFile, chainID.upper())
            
            print ":::here2: writing line to LAMWmatch"
            length=int(entries[2])-int(entries[1])+1
            out=groFile+' '+entries[1]+' '
            out+=entries[3]+' '+str(length)+' '+str(weight)+"\n"
            LAMWmatch.write(out)

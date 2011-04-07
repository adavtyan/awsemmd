#!/usr/bin/python
#rmbh11
#release version - v0.6 (corrects 0.5)

#USAGE: prepFrags.py database-prefix file.fasta > logfile

import sys,os

database=sys.argv[1]
inFASTA=sys.argv[2]
weight=1 #feature in match file

from Bio import SeqIO
inseq=SeqIO.read(inFASTA,'fasta')
print "processing: ",inseq.name
#print inseq.seq

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

##open match file
match=open('prepFrags.match','w')
match.write(inseq.id+"\n")
##FRAGMENT GENERATION LOOP
iterations=len(inseq.seq)-fragmentLength+1 #number of sliding window positions
for i in range(1,iterations+1): 
##select subrange
    rangeStart=i-1
    rangeEnd=i+fragmentLength-1
    subrange=str(inseq[rangeStart:rangeEnd].seq)
    fragment=open('/usr/tmp/fragment.fasta','w')
    print subrange
    fragment.write(subrange)
    fragment.close()
##submit PSI-BLAST -- run psiblast -help for explanation of format 6 output
#    exeline="psiblast -query_loc "+str(rangeStart)+"-"+str(rangeEnd)
    exeline="psiblast -num_iterations 1 -word_size 2 -evalue "+str(EvalueThreshold)
    exeline+=" -outfmt 6 -matrix BLOSUM62 -db "
    exeline+=database+" -query /usr/tmp/fragment.fasta"
    print exeline
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
        #skip gapped alignments:
        if ((delQuery-delAlign)==0):
            match.write(out)
##get pdbcodes and alignment ranges
##from script 'pdbget' (original author unknown)
#system ("wget ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb/$two$three/pdb$pdb.ent.gz");
#system ("gunzip pdb$pdb.ent.gz");
#system ("mv pdb$pdb.ent $pdb.pdb");
#foreach line in the output
#if this line introduced a new memory fragment
##write to match
#line 1: target name target length
#line 1+n: memory-name target-start memory-start fragment-length weight
##over match


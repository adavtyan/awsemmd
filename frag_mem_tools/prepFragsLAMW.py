#!/usr/bin/python
#rmbh11
#v0.72
#corects and builds on release version - v0.5

######################################
#USAGE: prepFragsLAMW.py database-prefix file.fasta > logfile
#Originally written by Ryan Hoffmann
#Modified by Weihua Zheng 06/01/2011

#######################################################################
#NOTE: Before running this script, please make sure the fasta 
#file contains only the sequences that have coordinates in the PDB file
######################################################################

import sys,os,re
from Pdb2GroLib import *
from Bio.PDB.Polypeptide import * #func three_to_one()

if len(sys.argv)!=5:
	print "\n prepFragsLAMW.py database-prefix file.fasta N_mem brain_damage_flag (1/0 for yes/no) > logfile \n\n"
	print "#######################################################################"
	print "#NOTE: Before running this script, please make sure the fasta file "
	print "#contains only the sequences that have coordinates in the PDB file"
	print "######################################################################"
	exit()

################################################################
## NoMissingAtoms function

def NoMissingAtoms(atom_list, residue_list, res_Start, pdbID, pdbFile):
	res_End = res_Start + len(residue_list) - 1
	#print "res:", res_Start, res_End
	p = PDBParser(PERMISSIVE=1)
	s = p.get_structure(pdbID, pdbFile)
	chains = s[0].get_list()
	chain_name = "A" 
	count_chain = 1 #do the following for just one chain

	keys_res = {}
	keys = {}

	for chain in chains:
		name = chain.get_id()
		if (name == '' or name == chain_name) and count_chain == 1:
	        	count_chain = 0     
			i = 0
			for res in chain:
		                res_index = res.get_id()[1]

		                if (res_index < res_Start ):
		                    continue
		                if (res_index > res_End   ):
		                    break

		                is_regular_res = res.has_id('N') and res.has_id('CA') and res.has_id('C')
				res_id = res.get_id()[0]
		                if not ((res_id ==' ' or res_id =='H_MSE') and is_regular_res):
					print 'Non-regular residue:', res.get_id()[0], 'at position', res_index,  'in pdb:', pdbID
					return False
				res_name = res.get_resname()
				#convert to 1-letter code
				if res_name == 'MSE':
					res_code = 'M'
				else:
					res_code = three_to_one(res_name)

				#Add sanity check, residues have to match the blast-out seq
				if ( res_code != residue_list[i] ):
					print "Mismatching residue in the PDB file:", pdbID, "letter code:", res_code
					return False

				i += 1

			        
				keys = {}
				if res_name == 'GLY':  #GLY has no CB atoms  
					keys['CB']=1
			        for atom in res:
			                atom_name = atom.get_name()
			                for target_atom_name in atom_list:
				               	if atom_name == target_atom_name:
			                 	        keys[target_atom_name]=1
							#print "matching:", atom_name
						if len(keys) == len(atom_list):
							break;
		                if len(keys) == len(atom_list):
				#	print "matching res:", res_index
	        	        	keys_res[res_index] = 1 

	if len(keys_res) == res_End - res_Start + 1:
		return True
	else:
		return False
## NoMissingAtoms function
################################################################

database=sys.argv[1]
fasta =sys.argv[2]
N_mem = int(sys.argv[3])
brain_damage = int(sys.argv[4])

inFASTA=open(fasta, 'r')
weight=1 #feature in match file

from Bio import SeqIO
inseq=SeqIO.read(inFASTA,'fasta')
print "processing: ",inseq.name
query=str(inseq.name)[0:4]

myhome = os.environ.get("HOME")
pdbDir = myhome + "/opt/script/PDBs/"
fLibDir = "./fraglib/"

fragmentLength=9 #needs to be an odd number
memoriesPerPosition=N_mem  #can be any integer > 0
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
iterations=len(inseq.seq)-fragmentLength+1 #number of sliding windows

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
    ##submit PSI-BLAST
    ##run "psiblast -help" for more details of output format (outfmt)
    exeline="psiblast -num_iterations 1 -word_size 2 -evalue "+str(EvalueThreshold)
    exeline+=" -outfmt '6 sseqid qlen slen qstart qend sstart send qseq sseq length gaps bitscore evalue' -matrix BLOSUM62 -db "
    exeline+=database+" -query fragment.fasta"
    print "executing:::"+exeline
    psiblastOut=os.popen(exeline).read()
    psiblastOut=psiblastOut.splitlines() #now an array
    print "Number of searched PDBs:  ", len(psiblastOut)
#    print psiblastOut
#   exit()
#    print "PDB INSEQ-START INSEQ-END MATCH-START MATCH-END EVALUE"
    for line in psiblastOut:#[0:memoriesPerPosition]:
        this=line.split()
	this.append(str(i))
	# 0:sseqid 1:qlen 2:slen 3:qstart 4:qend 5:sstart 6:send 7:qseq 8:sseq 9:length 10:gaps 11:bitscore 12:evalue 13:window_index
        queryStart=int(this[3])+rangeStart  #+int(this[6])
        queryEnd  =queryStart+int(this[4])-1
        #print this #[1],str(queryStart),str(queryEnd),this[8],this[9],this[11]
	this[3] = str(queryStart)
	this[4] = str(queryEnd)

	out = ' '.join(this)
	out+='\n'
	gaps = this[10]
	if(gaps == '0'): #skip gapped alignments
		match.write(out)
	
        #out=this[1]+' '+str(queryStart)+' '+str(queryEnd)+' '
#        out+=this[8]+' '+this[9]+' '+str(weight)+"\n"
#        delQuery=queryEnd-queryStart
#        delAlign=int(this[9])-int(this[8])
        #if residue ranges do not match, this alignment was gapped
        #skip gapped alignments:       ################################
#        if ((delQuery-delAlign)==0):
#            match.write(out)

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

from Bio.PDB.PDBParser import PDBParser
pdbparse=PDBParser(PERMISSIVE=1)

#atomLine=re.compile('\AATOM')

#Finding homologs
print inseq.seq
fragment=open('fragment.fasta','w')
fragment.write(str(inseq.seq))
fragment.close()
homo={} 
failed_pdb = {}
for pdbfull in unique:
    pdbID=pdbfull[0:4].lower()
    pdbIDsecond=pdbfull[1:2].lower()
    pdbIDthird=pdbfull[2:3].lower()
    chainID=pdbfull[4:5].lower()
    failed_pdb[pdbID] = 0
    homo[pdbID] = 0
    if not os.path.isfile(pdbDir+pdbID.upper()+".pdb"):
	##from script 'pdbget' (original author unknown)
        exeline="wget ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb/"
        exeline+=pdbIDsecond+pdbIDthird+"/pdb"+pdbID+".ent.gz"
        os.system(exeline);
        os.system("nice gunzip pdb"+pdbID+".ent.gz; mv pdb"+pdbID+".ent "+pdbDir+pdbID.upper()+".pdb");

    if not os.path.isfile(pdbDir+pdbID.upper()+".pdb"):
        print ":::Cannot build PDB for PDB ID, failed to download:"+pdbID.upper()
	failed_pdb[pdbID] = 1
	
        
if brain_damage == 1:
# blast the whole sequence to identify homologs
# Evalue 0.01 
    exeline="psiblast -num_iterations 1 -word_size 2 -evalue 0.005"
    exeline+=" -outfmt '6 sseqid slen bitscore score evalue' -matrix BLOSUM62 -db "
    exeline+=database+" -query fragment.fasta"
    print "brain damamge, finding homologs"
    print "executing::: "+exeline
    psiblastOut=os.popen(exeline).read()
    psiblastOut=psiblastOut.splitlines() #now an array
    #print psiblastOut
    for line in psiblastOut:
    	entries=line.split()
	print entries
	pdbfull = entries[0]
	pdbID = pdbfull[0:4].lower()
	homo[pdbID] = 1

iter=0
count = {}
#count number of mem per fragments
for i in range(1,iterations+1): 
	count[str(i)]=0

Missing_count = 0
for line in matchlines:
    iter+=1
    if not(iter==1):
        print ":::here: match line:"+line.rstrip('\n')
        entries=line.split()

	windows_index_str = entries[13]

        pdbfull=str(entries[0])
        pdbID=pdbfull[0:4].lower()
        pdbIDsecond=pdbfull[1:2].lower()
        pdbIDthird=pdbfull[2:3].lower()
        chainID=pdbfull[4:5].lower()
	groFile=fLibDir+pdbID+chainID+".gro"
	pdbFile=pdbDir+pdbID.upper()+".pdb"

	if failed_pdb[pdbID]:#failed-downloaded ones are still in matchlines, need to be ignored
		continue
	#ignore homologs
	if brain_damage and  homo[pdbID]:
		continue
	atoms_list = ('CA', 'CB')
	residue_list = entries[8]  ##sseq
	res_Start = int(entries[5])
	res_End   = int(entries[6])
	print "start: ", res_Start, "end: ", res_End

	#check missing atoms
	##have to check residue list, not residue index.

	if count[windows_index_str] >= N_mem:
		continue

	if NoMissingAtoms(atoms_list, residue_list, res_Start, pdbID, pdbFile):
	        if os.path.isfile(pdbFile):
        	    print ":::convert: "+pdbFile+" --> "+groFile
        	    Pdb2Gro(pdbFile, groFile, chainID.upper())
		    count[windows_index_str] += 1
            
        	    print ":::here2: writing line to LAMWmatch\n"
        	    length=res_End - res_Start + 1  #int(entries[2])-int(entries[1])+1
	            out=groFile+' '+entries[3]+' '
        	    out+=entries[5]+' '+str(length)+' '+str(weight)+"\n"
        	    LAMWmatch.write(out)
		else:
		    print pdbFile, "does not exist! Go figure..."
	else:
		Missing_count += 1
print "MemPerPosition: ", count
print "Number of failed downloaded PDB: ", len(failed_pdb)
print "Number of PDB with Missing atoms: ", Missing_count

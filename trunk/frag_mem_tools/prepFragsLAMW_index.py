#!/usr/bin/python
#rmbh11
#v0.72
#corects and builds on release version - v0.5

######################################
#USAGE: prepFragsLAMW.py database-prefix file.fasta > logfile
#Originally written by Ryan Hoffmann from Wolynes group

#Modified by Weihua Zheng 06/01/2011

#######################################################################
#NOTE: Before running this script, please make sure the fasta 
#file contains only the sequences that have coordinates in the PDB file
######################################################################

import sys,os,re
from IndexPdb import *
from Pdb2GroLib import *
from Bio.PDB.Polypeptide import * #func three_to_one()
from Bio import SeqIO

if len(sys.argv)!=5:
	print "\n prepFragsLAMW.py database-prefix file.fasta N_mem brain_damage_flag (1/0 for yes/no) > logfile \n\n"
	print "#######################################################################"
	print "#NOTE: Before running this script, please make sure the fasta file "
	print "#contains only the sequences that have coordinates in the PDB file"
	print "######################################################################"
	exit()

################################################################
## NoMissingAtoms function

def NoMissingAtoms(atom_list, residue_list, res_Start, pdbID, ch_name, pdbFile):
	res_End = res_Start + len(residue_list) - 1
	p = PDBParser(PERMISSIVE=1)
	s = p.get_structure(pdbID, pdbFile)
	chains = s[0].get_list()
	if ch_name == '':
		ch_name = "A" 

	keys_res = {}
	keys = {}

	for chain in chains:
		if chain.get_id()==ch_name:
			i = 0
			for res in chain:
		                res_index = res.get_id()[1]
		                if (res_index < res_Start ):
		                    continue
		                if (res_index > res_End and i == 0 ):
				   print "Residue index shifted: ", res_index, "mismatch: ", res_Start
		                   return False
		                if (res_index > res_End   ):
		                    break

		                is_regular_res = res.has_id('N') and res.has_id('CA') and res.has_id('C')
				res_id = res.get_id()[0]
		                if not (res_id ==' ' or res_id =='H_MSE' or res_id =='H_M3L' or res_id=='H_CAS') and is_regular_res :
					print 'Discard Fragment: Non-regular residue:', res.get_id()[0], 'at position', res_index,  'in pdb:', pdbID
					return False
				res_name = res.get_resname()
				#convert to 1-letter code
				if res_name == 'MSE':
					res_code = 'M'
				elif res_name == 'M3L':
					res_code = 'K'
				elif res_name == 'CAS':
					res_code = 'C'

				else:
					res_code = three_to_one(res_name)

				#Add sanity check, residues have to match the blast-out seq
				if ( res_code != residue_list[i] ):
					print "Mismatching residue in the PDB file:", pdbID, "residue :", res_code
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
		print "Missing CA or CB in the residues for PDB ", pdbID, ch_name
		print "Good residues: "
		for j in keys_res:
			print j
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

myhome  = os.environ.get("HOME")
pdbDir  = myhome + "/opt/script/PDBs/"
fLibDir = "./fraglib/"
indexDir = myhome + "/opt/script/Indices/"
#fLibDir = myhome + "/opt/script/fraglib/"
pdbSeqres= myhome + "/opt/script/pdb_seqres.txt"
fasta_database = database+".fasta"

# Index database fasta file
if not os.path.isfile(fasta_database):
	print "Can't find database fasta file"
        exit()
seq_records = SeqIO.index(fasta_database, "fasta")

fragmentLength=9 #needs to be an odd number
memoriesPerPosition=N_mem  #can be any integer > 0
EvalueThreshold=10000 #needs to be large enough that PSI-BLAST returns at least memoriesPerPosition

##SANITY CHECKING
#is length greater than fragmentLength?
if(len(inseq.seq) < fragmentLength):
    print "Exception::query sequence is smaller than "+str(fragmentLength)+" residues"
    print "This version has no means to handle smaller queries"
    sys.exit()

##Create necessary directories

if not os.path.exists(indexDir): os.makedirs(indexDir)

if not os.path.exists(pdbDir): os.makedirs(pdbDir)
if not os.path.exists(fLibDir): os.makedirs(fLibDir)
if not os.path.exists(pdbDir) or not os.path.exists(fLibDir) or not os.path.exists(indexDir) :

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
    exeline+=" -outfmt '6 sseqid qstart qend sstart send qseq sseq length gaps bitscore evalue' -matrix BLOSUM62 -db "
    exeline+=database+" -query fragment.fasta"
    print "executing:::"+exeline
    psiblastOut=os.popen(exeline).read()
    psiblastOut=psiblastOut.splitlines() #now an array
    print "Number of searched PDBs:  ", len(psiblastOut)
#   print psiblastOut
#   exit()
#    print "PDB INSEQ-START INSEQ-END MATCH-START MATCH-END EVALUE"
    for line in psiblastOut:#[0:memoriesPerPosition]:
        this=line.split()
	this.append(str(i))
	print this
	# 0:sseqid 1:qlen 2:slen 3:qstart 4:qend 5:sstart 6:send 7:qseq 8:sseq 9:length 10:gaps 11:bitscore 12:evalue 13:window_index
        queryStart=int(this[1])+rangeStart  #+int(this[6])
        queryEnd  = rangeStart + int(this[2])
        #print this #[1],str(queryStart),str(queryEnd),this[8],this[9],this[11]
	this[1] = str(queryStart)
	this[2] = str(queryEnd)

	out = ' '.join(this)
	out+='\n'
	gaps = this[8]
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

log_match=open('log.mem','w')

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
# blast the whole sequence to identify homologs Evalue 0.005 
    exeline="psiblast -num_iterations 1 -word_size 2 -evalue 0.005"
    exeline+=" -outfmt '6 sseqid slen bitscore score evalue' -matrix BLOSUM62 -db "
    exeline+=database+" -query fragment.fasta"
    print "brain damamge, finding homologs"
    print "executing::: "+exeline
    homoOut=os.popen(exeline).read()
    homoOut=homoOut.splitlines() #now an array
    for line in homoOut:
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
Missing_pdb = {}
fastFile="./tmp.fasta"
for line in matchlines:
    iter+=1
    if not(iter==1):
        #print ":::here: match line:"+line.rstrip('\n')
        entries=line.split()

	windows_index_str = entries[11]
	if count[windows_index_str] >= N_mem:
		continue

        pdbfull=str(entries[0])
        pdbID=pdbfull[0:4].lower()
        pdbIDsecond=pdbfull[1:2].lower()
        pdbIDthird=pdbfull[2:3].lower()
        chainID=pdbfull[4:5].lower()
	groFile=fLibDir+pdbID+chainID+".gro"
	groName=pdbID+chainID+".gro"
	pdbFile=pdbDir+pdbID.upper()+".pdb"

	indexFile=indexDir+pdbID+chainID+".index"


	if failed_pdb[pdbID]:#failed-downloaded ones are still in matchlines, need to be ignored
		continue
	#ignore homologs
	if brain_damage and  homo[pdbID]:
		print pdbID, " is a homolog, discard"
		continue
	atoms_list = ('CA', 'CB')
	residue_list = entries[6]  ##sseq
	res_Start = int(entries[3])
	res_End   = int(entries[4])
	print "start: ", res_Start, "end: ", res_End

	#check missing atoms
	##have to check residue list, not residue index.
	#if NoMissingAtoms(atoms_list, residue_list, res_Start, pdbID, chainID.upper(), pdbFile):

	#Do I have the index file?
	#No, write it
	if not os.path.isfile(indexFile):
		#generate fasta file
		seq_id = pdbID.upper() + chainID.upper()
		handle = open(fastFile, "w")
		SeqIO.write(seq_records[seq_id], handle, "fasta")
		handle.close()
#		print str(seq_records[seq_id].seq)

#		if not os.path.isfile(pdbSeqres):
#			print "Need to download pdb_seqres.txt from PDB!"
#			print "ftp://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt"
#			print "Copy to $HOME/opt/script/"
#			exit()
#		fastaFile=pdbID+'_'+chainID.upper()
#		exeline="grep -A1 "+fastaFile+" "+pdbSeqres+" > ./tmp.fasta"
#		os.popen(exeline)

		#write index file
		if os.path.getsize('tmp.fasta') > 0 :
			print "Writing indexFile: ", indexFile
			writeIndexFile(fastFile, pdbFile, indexFile, chainID.upper())

	#Read index file
	if not os.path.isfile(indexFile):
		print "Can't create index file, ignore and go on!"
		continue

	index=open(indexFile,'r')
	#create new_index for frag_seq starting position
	line_count = 0
	flag=' '
	index_shift=0
	#read and get the flag
	indexlines=list()
	for index_line in index.readlines():
		#print index_line
		tmp_line=index_line.split()
		line_count += 1
		indexlines.append(index_line)
		if line_count == 1 : #first line is the flag
			flag = tmp_line[0] #index_line
		if flag == "SHIFT" and line_count == 2 :
			print "shift: ", tmp_line[0]
			index_shift = int(tmp_line[0])

	r_list = '' #list()
	if flag == "SKIP":
		Missing_pdb[pdbID] = 1
		Missing_count += 1
		print "***********", flag
		print "SKIP pdb:", pdbID+chainID
		continue
	elif flag == "FULLMATCH":
		new_index = int(entries[3])
		r_list = residue_list
		print "***********", flag
	elif flag == "SHIFT":
		new_index = int(entries[3]) + index_shift
		r_list = residue_list
		print "***********", flag
	elif flag == "INDEXED":
		print "***********", flag
		#check if there is gaps
		count_flag = 0
		line_count1=0
		for index_line in indexlines:
		  line_count1 += 1
		  if not line_count1==1:
			index_entries = index_line.split()
			seq_id = int(index_entries[0])
			res_id = int(index_entries[1])
			#print "seq_id:", seq_id, "res_id:", res_id
			if seq_id < res_Start:
				continue
			if seq_id > res_End:
				break
			if res_id == -1:
				print "Missing residues in PDB: ", pdbID+chainID
				break
			if count_flag == 0 :
				new_index = res_id
				count_flag += 1
			res_nm = index_entries[2]
			#print "res_name: ", res_nm
			#r_list.append(res_nm)
			r_list+=res_nm
			#print r_list

	if r_list != residue_list:
		print "Missing residues: ", pdbID+chainID, residue_list, " incomplete: ", r_list
		#print
		Missing_pdb[pdbID] = 1
		Missing_count += 1
		continue

        if os.path.isfile(pdbFile): 
		    if not os.path.isfile(groFile) :
        	        Pdb2Gro(pdbFile, groFile, chainID.upper())
        	    print ":::convert: "+pdbFile+" --> "+groFile
		    count[windows_index_str] += 1
            
        	    print ":::here2: writing line to LAMWmatch\n"
        	    length=res_End - res_Start + 1  
	            out=groFile+' '+entries[1]+' ' #queue start
        	    #out+=entries[3]+' '+str(length)+' '+str(weight)+"\n" #frag_seq start
        	    out+= str(new_index) +' '+str(length)+' '+str(weight)+"\n" #frag_seq start
        	    LAMWmatch.write(out)
		    #out1 = out
		    out1 = windows_index_str 
		    out1 += ' '+str(count[windows_index_str])
		    out1 += ' '+entries[9]+' '+entries[10]+' '+groName
		    out1 += ' '+entries[1]+' '+ str(new_index) +' '+str(length)+' '+str(weight)+"\n"
		    log_match.write(out1)
	else:
		    print pdbFile, "does not exist! Go figure..."

if brain_damage == 1:
  for line in homoOut:
    entries=line.split()
    print "HOMOLOGS:::"
    print entries
print "memories per position that is fewer than expected:"  
for i in count:
  if count[i] < N_mem:
    print i, count[i]

#print "MemPerPosition: ", count
print "Number of blasted PDB: ", len(failed_pdb)
print "Number of failed downloaded PDB: ", sum(failed_pdb.values())
print "Number of PDB with Missing atoms: ", len(Missing_pdb)
print "Discarded fragments with Missing atoms: ", Missing_count 

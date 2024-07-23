#!/shared/local/bin/python

import sys
from VectorAlgebra import *

def calc_dihedral_angle(p1, p2, p3, p4):
    v1 = vector(p1, p2)
    v2 = vector(p2, p3)
    v3 = vector(p3, p4)
    return dihedral_angle(v1, v2, v3)

def calc_angle(p1, p2, p3):
    v1 = vector(p1, p2)
    v2 = vector(p3, p2)
    return vangle(v1, v2)

def calc_bond(p1, p2):
    v = vector(p1, p2)
    return vabs(v)

def three2one(prot):
    """ translate a protein sequence from 3 to 1 letter code"""
    
    code = {"GLY" : "G", "ALA" : "A", "LEU" : "L", "ILE" : "I",
            "ARG" : "R", "LYS" : "K", "MET" : "M", "CYS" : "C",
            "TYR" : "Y", "THR" : "T", "PRO" : "P", "SER" : "S",
            "TRP" : "W", "ASP" : "D", "GLU" : "E", "ASN" : "N",
	    "GLN" : "Q", "PHE" : "F", "HIS" : "H", "VAL" : "V"}
    
    newprot = ""
    for a in prot:
        newprot += code.get(a, "?")
    
    return newprot


def checkIfNative(ires, jres):
    xyz_CAi = ires['CA'].get_coord()
    xyz_CAj = jres['CA'].get_coord()
    v = vector(xyz_CAi, xyz_CAj)
    r = vabs(v)
    if r<12.0: return True
    else: return False

if len(sys.argv)==1:
    print "\nReadingPDBFile.py PDB_Id [Output_file [-s]]\n"
    print "-s\tSplit into files for each chain"
#    sys.argv.append("1BG8")
    exit()

from Bio.PDB.PDBParser import PDBParser

p = PDBParser(PERMISSIVE=1)

struct_id = sys.argv[1]
filename = struct_id + ".pdb"

######## Weihua ###########
#shadow_list is the modified output from SMOG script SCM.jar
#it's a 2 col data contains all the native contacts
# res_i res_j
def IfNative(res_i, res_j, shadow_list):
	for sublist in shadow_list:
		if sublist[0] == str(res_i) and sublist[1] == str(res_j):
			return True
	return False

shadowfile = struct_id + ".shadow"
fh_shadow = open(shadowfile, "r")
lines = fh_shadow.readlines()
ls_shadow = []
for line in lines:
	ls = line.strip().split()
	ls_shadow.append(ls)
fh_shadow.close()

############################

splite = False
for av in sys.argv:
    if av=="-s":
        splite = True
        sys.argv.remove(av)
        break

output_fn = ""
if len(sys.argv)>2: output_fn = sys.argv[2]
if output_fn[-5:]==".data": output_fn = output_fn[:-5]

if output_fn!="" and not splite:
    out = open( (output_fn+".data"), 'w' )

k_bond = 100
k_angle = 20
k_dihedral = [1, 0.5]
epsilon = 1
epsilon2 = 1
sigma0 = 4

xyz_CA1 = []
xyz_CA2 = []
xyz_CA3 = []
xyz_CA4 = []

s = p.get_structure(struct_id, filename)
chains = s[0].get_list()
for ch in chains:
    sequance = []
    bonds = []
    angles = []
    dihedrals = []
    if output_fn!="":
	pass
    else:
        print "Chain:", ch.get_id()
    four_res = [None, None, None, None]
    all_res = []
    for res in ch:
#        is_regular_res = res.has_id('N') and res.has_id('CA') and res.has_id('C')
	is_regular_res = res.has_id('CA') and res.has_id('O')
        if res.get_id()[0]==' ' and is_regular_res:
            all_res.append(res)
            four_res.append(res)
            p_res = four_res.pop(0)
            sequance.append(res.get_resname())
            if four_res[2] and four_res[3]:
                if four_res[2].get_id()[1]+1!=four_res[3].get_id()[1]:
                    print "Error: Wrong residue order"
                xyz_CA3 = four_res[2]['CA'].get_coord()
                xyz_CA4 = four_res[3]['CA'].get_coord()
                r = calc_bond(xyz_CA3, xyz_CA4)
                bonds.append(r)
            if four_res[1] and four_res[2] and four_res[3]:
                if four_res[1].get_id()[1]+1!=four_res[2].get_id()[1]:
                    print "Error: Wrong residue order"
                xyz_CA2 = four_res[1]['CA'].get_coord()
                theta = calc_angle(xyz_CA2, xyz_CA3, xyz_CA4)
                angles.append(theta)
            if four_res[0] and four_res[1] and four_res[2] and four_res[3]:
                if four_res[0].get_id()[1]+1!=four_res[1].get_id()[1]:
                    print "Error: Wrong residue order"
                xyz_CA1 = four_res[0]['CA'].get_coord()
                phi = calc_dihedral_angle(xyz_CA1, xyz_CA2, xyz_CA3, xyz_CA4)
                dihedrals.append(phi)

    isNative = []
    sigma = []
    for i in range( 0, len(all_res) ):
	zero = []
	zerob = []
	for ii in range( 0, len(all_res) ):
	    zero.append(0)
	    zerob.append(False)
	isNative.append(zerob)
	sigma.append(zero)
        for j in range( i+4, len(all_res) ):
            ires = all_res[i]
	    jres = all_res[j]
            #isNative[i][j] = checkIfNative(ires, jres)

	    isNative[i][j] = IfNative(i+1, j+1, ls_shadow)
	    
	    if isNative[i][j]:
                xyz_CAi = ires['CA'].get_coord()
                xyz_CAj = jres['CA'].get_coord()
		v = vector(xyz_CAi, xyz_CAj)
                sigma[i][j] = vabs(v)
	    else:
                sigma[i][j] = sigma0

    if output_fn!="":
        if splite and len(sequance)==0: continue
        if splite:
	    if len(chains)==1:
		file_name = output_fn+".data"
	    else:
		file_name = output_fn+"_"+ch.get_id()+".data"
            out = open( file_name, 'w' )
	
	out.write('[Go-Model_LJ]\n')
	out.write(str(epsilon))
	out.write(' ')
	out.write(str(epsilon2))
	out.write('\n\n')
		
	out.write('[Bonds]\n')
	out.write(str(k_bond))
	out.write('\n')
	for ir in bonds:
		out.write( str(round(ir, 5)) )
		if bonds.index(ir)!=len(bonds)-1:
			out.write(' ')
	out.write('\n\n')
	
	out.write('[Angles]\n')
	out.write(str(k_angle))
	out.write('\n')
	for ia in angles:
		out.write( str(round(ia, 5)) )
		if angles.index(ia)!=len(angles)-1:
			out.write(' ')
	out.write('\n\n')
	
	out.write('[Dihedrals]\n')
	out.write(str(k_dihedral[0]))
	out.write(' ')
	out.write(str(k_dihedral[1]))
	out.write('\n')
	for idha in dihedrals:
		out.write( str(round(idha, 4)) )
		if dihedrals.index(idha)!=len(dihedrals)-1:
			out.write(' ')
	out.write('\n\n')
	
	out.write('[Contacts]\n')
	for i in range( 0, len(all_res)-4 ):
        	for j in range( i+4, len(all_res) ):
			if isNative[i][j]:
				out.write( '1' )
			else:
				out.write( '0' )
			if j!=len(all_res)-1:
				out.write( ' ' )
		out.write('\n')
	out.write('\n\n')
	for i in range( 0, len(all_res)-4 ):
        	for j in range( i+4, len(all_res) ):
			out.write( str(round(sigma[i][j], 4)) )
			if j!=len(all_res)-1:
				out.write(' ')
		out.write('\n')
	out.write('\n\n')
	
        if not splite:
            out.write('\n')
        if splite:
            out.close()
    else:
        print three2one(sequance)
        print "[Go-Model_LJ]"
	print str(epsilon), str(epsilon2)
	print '\n'
		
	print "[Bonds]"
	print str(k_bond)
	for ir in bonds:
		print str(round(ir, 5)), ' ',
	print '\n'
	
	print "[Angles]"
	print str(k_angle)
	for ia in angles:
		print str(round(ia, 4)), ' ',
	print '\n'
	
	print "[Dihedrals]"
	print str(k_dihedral[0]), str(k_dihedral[1])
	for idha in dihedrals:
		print str(round(idha, 4)), ' ',
	print '\n'
	
	print "[Contacts]"
	for i in range( 0, len(all_res)-4 ):
        	for j in range( i+4, len(all_res) ):
			if isNative[i][j]:
				print 1,
			else:
				print 0,
		print '\n',
	print '\n',
	for i in range( 0, len(all_res)-4 ):
        	for j in range( i+4, len(all_res) ):
			print str(round(sigma[i][j], 4)),
		print '\n',
        print '\n'

if output_fn!="" and not splite:
    out.close()

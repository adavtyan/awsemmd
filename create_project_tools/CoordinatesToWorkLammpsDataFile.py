#!/usr/bin/python

# ----------------------------------------------------------------------
# Copyright (2010) Aram Davtyan and Garegin Papoian

# Papoian's Group, University of Maryland at Collage Park
# http://papoian.chem.umd.edu/

# Last Update: 03/04/2011
# ----------------------------------------------------------------------
import os
import sys

class Atom:
    def __init__(self, no, ch, res, ty, q, x, y, z):
        self.no = no
        self.ch = ch
        self.res = res
        self.ty = ty
        self.q = q
        self.x = x
        self.y = y
        self.z = z

    def write_(self, f):
        space11 = "           "
        f.write( (space11+str(self.no))[-12:] + "\t" )
        f.write( "\t".join([ str(self.ch), str(self.res), str(self.ty), str(self.q), str(self.x), str(self.y), str(self.z) ]) )
        f.write( "\n" )

class Bond:
    def __init__(self, no, ty, I, J):
        self.no = no
        self.ty = ty
        self.I = I
        self.J = J

    def write_(self, f):
        f.write( (space11+str(self.no))[-12:] + "\t" )
        f.write( "\t".join([ str(self.ty), str(self.I), str(self.J) ]) )
        f.write( "\n" )

inp_file = ""
out_file = ""
if len(sys.argv)>1: inp_file = sys.argv[1]
if len(sys.argv)>2: out_file = sys.argv[2]

if inp_file=="":
    print "\nCoordinatesToLammpsDataFile.py input_file [output_file] [-b] [-go]\n\n"
    print "\t-b\tadd bonds between CA & CA, CA & O and CA & CB in the case of coarse graining\n"
    print "\t-go\tcoarse-grained setup\n\n"
    exit()

cg_bonds = False
go = False
for cl in sys.argv[3:]:
	if cl == '-b': cg_bonds = True
	if cl == '-go': go = True
		

seq_file = "sequance.seq"    
lammps_out_file = "file.in"
if out_file[:5]=="data.":
    lammps_out_file = out_file[5:] + ".in"
    seq_file = out_file[5:] + ".seq"
elif out_file[-5:]==".data":
    lammps_out_file = out_file[:-5] + ".in"
    seq_file = out_file[:-5] + ".seq"
else:
    lammps_out_file = out_file + ".in"
    out_file = "data." + out_file
    seq_file = out_file + ".seq"

cg = True

xlo = -200.0
xhi = 200.0
ylo = -200.0
yhi = 200.0
zlo = -200.0
zhi = 200.0
masses = [12.0, 14.0, 16.0, 12.0, 1.0]
if cg and not go:
	masses = [27.0, 14.0, 28.0, 60.0, 2.0]
n_atom_types = 5
if cg:
	if cg_bonds: n_bond_types = 5
	else: n_bond_types = 0
else: n_bond_types = 7

last_nos = { 'N' : 0, 'C-Alpha' : 0, 'C-Prime' : 0, 'O' : 0 }
last_chno = { 'N' : 0, 'C-Alpha' : 0, 'C-Prime' : 0, 'O' : 0 }

n_atoms = 0
n_bonds = 0
n_res = 0
group_id = 0
atoms = []
bonds = []
groups = []

fix_string = "2 alpha_carbons backbone beta_atoms oxygens fix_backbone_coeff.data " + seq_file
if go:
	fix_string = "2 alpha_carbons gomodel fix_gomodel_coeff.data"

groups.append(["alpha_carbons", "id"])
if not go:
	groups.append(["beta_atoms", "id"])
	groups.append(["oxygens", "id"])

inp = open(inp_file)
atom_type = 0
for l in inp:
    l = l.strip().split()
    if len(l)==6:
        print "Input file lacks description field!"
        exit()

    desc = l[6]
    chain_no = l[1]
    if not go:
	if desc == 'C-Beta' or desc == 'H-Beta' or desc == 'C-Alpha' or desc == 'O':
		n_atoms += 1
    else:
	if desc == 'C-Alpha':
		n_atoms += 1

    if not go:
        if desc == 'N0' or desc == 'N':
            atom_type = 2
            if last_nos['C-Prime']!=0 and last_chno['C-Prime']==chain_no and not cg:
                n_bonds += 1
                bonds.append( Bond(n_bonds, 3, last_nos['C-Prime'], n_atoms) )
            desc = 'N'
            last_nos[desc] = n_atoms
            last_chno[desc] = chain_no
            n_res += 1
        elif desc == 'C-Alpha':
            if last_nos['N']!=0 and last_chno['N']==chain_no and not cg:
                n_bonds += 1
                bonds.append( Bond(n_bonds, 1, last_nos['N'], n_atoms) )
            if cg and cg_bonds:
                if last_nos['C-Alpha']!=0 and last_chno['C-Alpha']==chain_no:
                    n_bonds += 1
                    bonds.append( Bond(n_bonds, 1, last_nos['C-Alpha'], n_atoms) )
                if last_nos['O']!=0 and last_chno['O']==chain_no:
                    n_bonds += 1
                    bonds.append( Bond(n_bonds, 3, last_nos['O'], n_atoms) )
            atom_type = 1
            last_nos[desc] = n_atoms
            last_chno[desc] = chain_no
            group_id = 1
        elif desc == 'C-Prime':
            if last_nos['C-Alpha']!=0 and last_chno['C-Alpha']==chain_no and not cg:
                n_bonds += 1
                bonds.append( Bond(n_bonds, 2, last_nos['C-Alpha'], n_atoms) )
            atom_type = 1
            last_nos[desc] = n_atoms
            last_chno[desc] = chain_no
        elif desc == 'O':
            if last_nos['C-Prime']!=0 and last_chno['C-Prime']==chain_no and not cg:
                n_bonds += 1
                bonds.append( Bond(n_bonds, 6, last_nos['C-Prime'], n_atoms) )
            if cg and cg_bonds:
                    if last_nos['C-Alpha']!=0 and last_chno['C-Alpha']==chain_no:
                        n_bonds += 1
                        bonds.append( Bond(n_bonds, 2, last_nos['C-Alpha'], n_atoms) )
            atom_type = 3
            last_nos[desc] = n_atoms
            last_chno[desc] = chain_no
            group_id = 3
        elif desc == 'C-Beta':
            if last_nos['C-Alpha']!=0 and (not cg or cg_bonds):
                n_bonds += 1
                bonds.append( Bond(n_bonds, 4, last_nos['C-Alpha'], n_atoms) )
            atom_type = 4
            group_id = 2
        elif desc == 'H-Beta':
            if last_nos['C-Alpha']!=0 and (not cg or cg_bonds):
                n_bonds += 1
                bonds.append( Bond(n_bonds, 5, last_nos['C-Alpha'], n_atoms) )
            atom_type = 5
            group_id = 2
        elif desc == 'O-In-The-End':
            if last_nos['C-Prime']!=0 and not cg:
                n_bonds += 1
                bonds.append( Bond(n_bonds, 7, last_nos['C-Prime'], n_atoms) )
            atom_type = 3

    if not go:
        if desc == 'C-Beta' or desc == 'H-Beta' or desc == 'C-Alpha' or desc == 'O':
#            n_atoms += 1
            atoms.append( Atom(n_atoms, chain_no, n_res, atom_type, 0.0, float(l[3]), float(l[4]), float(l[5])) )
            groups[group_id - 1].append(str(n_atoms))
    else:
        if desc == 'C-Alpha':
	    atom_type = 1
	    n_res += 1
            atoms.append( Atom(n_atoms, chain_no, n_res, atom_type, 0.0, float(l[3]), float(l[4]), float(l[5])) )
            groups[group_id - 1].append(str(n_atoms))
inp.close()

if go:
	n_atoms = len(atoms)
	n_bonds = 0
	n_bond_types = 0
	n_atom_types = 1
	masses = [118.0]

space11 = "           "
out = open(out_file,'w')
out.write("LAMMPS protein data file\n\n")

out.write( (space11+str(n_atoms))[-12:] + "  atoms\n" )
out.write( (space11+str(n_bonds))[-12:] + "  bonds\n" )
out.write( space11 + "0  angles\n" )
out.write( space11 + "0  dihedrals\n" )
out.write( space11 + "0  impropers\n\n" )

out.write( (space11+str(n_atom_types))[-12:] + "  atom types\n" )
out.write( (space11+str(n_bond_types))[-12:] + "  bond types\n" )
out.write( space11 + "0  angle types\n" )
out.write( space11 + "0  dihedral types\n" )
out.write( space11 + "0  improper types\n\n" )

out.write ( "\t".join([ str(xlo), str(xhi), "xlo xhi\n" ]) )
out.write ( "\t".join([ str(ylo), str(yhi), "ylo yhi\n" ]) )
out.write ( "\t".join([ str(zlo), str(zhi), "zlo zhi\n\n" ]) )

out.write( "Masses\n\n" )
for i in range(0, len(masses)):
    out.write( (space11+str(i+1))[-12:] + "\t" + str(masses[i]) + "\n" )
out.write( "\n" )

out.write( "Atoms\n\n" )
for iAtom in atoms:
    iAtom.write_(out)
out.write( "\n" )

if cg and cg_bonds and not go:
	out.write( "Bond Coeffs\n\n" )
	out.write( space11 + "1\t200.0\t3.77\n" )
	out.write( space11 + "2\t200.0\t2.41\n" )
	out.write( space11 + "3\t200.0\t2.50\n" )
	out.write( space11 + "4\t200.0\t1.54\n" )
	out.write( space11 + "5\t200.0\t1.54\n" )

if (cg_bonds or not cg) and not go:
	out.write( "Bonds\n\n" )
	for iBond in bonds:
	    iBond.write_(out)
	out.write( "\n" )
out.close()



groups_string = ""
for igroup in groups:
    groups_string += "group\t\t" + " ".join(igroup) + "\n\n"

bonds_string = ""
if cg and cg_bonds and not go:
	bonds_string = "bond_style harmonic"

pair_string = ""
if cg and not go:
	pair_string = "pair_style vexcluded 2 3.5 3.5"

pair_coeff_string = ""
if cg and not go:
	pair_coeff_string = "pair_coeff * * 0.0\n"
	pair_coeff_string += "pair_coeff 1 1 20.0 3.5 4.5\n"
	pair_coeff_string += "pair_coeff 1 4 20.0 3.5 4.5\n"
	pair_coeff_string += "pair_coeff 4 4 20.0 3.5 4.5\n"
	pair_coeff_string += "pair_coeff 3 3 20.0 3.5 3.5\n"

replace_rules = [ ["``read_data_file",  "read_data " +  out_file],
                  ["``groups", groups_string],
		  ["``bonds", bonds_string],
		  ["``main_fix", fix_string],
		  ["``pair_interactions", pair_string],
		  ["``pair_coeff", pair_coeff_string] ]
myhome = os.environ.get("HOME")
inp = open(myhome + "/opt/script/inFilePattern.data")
inFile = inp.read()
inp.close()

for ir in replace_rules:
    inFile = inFile.replace(ir[0], ir[1])

out = open(lammps_out_file,'w')
out.write(inFile)
out.close()

#out = open(groups_out_file,'w')
#for igroup in groups:
#    out.write( "group\t\t" )
#    out.write( " ".join(igroup) )
#    out.write( "\n\n" )
#out.close()

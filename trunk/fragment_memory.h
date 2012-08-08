/* ----------------------------------------------------------------------
Copyright (2010) Aram Davtyan and Garegin Papoian

Papoian's Group, University of Maryland at Collage Park
http://papoian.chem.umd.edu/

Last Update: 03/04/2011
------------------------------------------------------------------------- */

// fragment_memory.h

/*#ifndef SE_MAP
#define SE_MAP

// {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"};
// {"A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"};
int fm_se_map[] = {0, 0, 4, 3, 6, 13, 7, 8, 9, 0, 11, 10, 12, 2, 0, 14, 5, 1, 15, 16, 0, 19, 17, 0, 18, 0};

// Four letter classes
// 1) SHL: Small Hydrophilic (ALA, GLY, PRO, SER THR) or (A, G, P, S, T) or {0, 7, 14, 15, 16}
// 2) AHL: Acidic Hydrophilic (ASN, ASP, GLN, GLU) or (N, D, Q, E) or {2, 3, 5, 6}
// 3) BAS: Basic (ARG HIS LYS) or (R, H, K) or {1, 8, 11}
// 4) HPB: Hydrophobic (CYS, ILE, LEU, MET, PHE, TRP, TYR, VAL) or (C, I, L, M, F, W, Y, V)  or {4, 9, 10, 12, 13, 17, 18, 19}
int four_letter_map[] = {1, 3, 2, 2, 4, 2, 2, 1, 3, 4, 4, 3, 4, 4, 1, 1, 1, 4, 4, 4};

#endif*/

#ifndef FRAGMENT_MEMORY_H
#define FRAGMENT_MEMORY_H

class Fragment_Memory {
public:
  Fragment_Memory(int p, int pf, int l, double w, char *fname, bool vec_fm_flag=false);
  ~Fragment_Memory();
  int pos;    // Position of the first residue in target sequance
  int mpos;   // Middle residue position (in target)
  int fpos;   // Position of the fragment in the library protein
  int len;    // Length of the fragment
  double weight;  // Weight for the particular fragment
  char *se;
  int vfm_flag; // Vector Fragment Memory flag
  char getSe(int resno); // residue name in one letter code corresponding to resno in target
  double Rf(int ires, int iatom, int jres, int jatom); // distance between atoms in residues i and j
  double VMf(int ires, int jres); // Angle between normalized CA->CB vectors of residues i and j
  int resType(int resno); // return type of fragment residue corresponding to resno in target
  char ThreeLetterToOne(char *tl_resty);
  int min(int, int);
  int max(int, int);
  int error;
  enum Errors {ERR_NONE=0, ERR_FILE, ERR_ATOM_COUNT, ERR_RES, ERR_CALL, ERR_VFM_GLY};
  enum FM_AtomTypes { FM_CA=1, FM_CB };
private:
  // rf[0][i][j] are CA-CA (i<j) and CB-CB (i>j) distances
  // rf[1][ica][jcb] are CA-CB distances
  double **rf[2];
  double **vmf;
  inline double R(double *r1, double *r2);
  inline double VM(double *r11, double *r12, double *r21, double *r22);
};

#endif

// --------------------------------------------------------------------//

#ifndef GAMMA_ARRAY_H
#define GAMMA_ARRAY_H

class Gamma_Array {
public:
  Gamma_Array(char *fname);
  ~Gamma_Array();
  double getGamma(int ires, int jres);
  double getGamma(int ires_type, int jres_type, int ires, int jres);
  double getGamma(int ires_type, int jres_type, int ifres_type, int jfres_type, int ires, int jres);
  bool fourResTypes() { return frag_resty; }
  static bool isEmptyString(char *str);
  int minSep();
  int maxSep();
  int error;
  char *se;
  enum Errors {ERR_NONE=0, ERR_FILE, ERR_CLASS_DEF, ERR_GAMMA, ERR_G_CLASS, ERR_ASSIGN, ERR_CALL};
private:
  double *gamma;
  int nseq_cl; // Number of sequance seperation classes
  int nres_cl; // Number of residue classes
  bool frag_resty; // True if gamma also depends on fragment residue types
  int i_sep[10];
  bool allocated;
  enum Four_Letter_Types {SHL=1, AHL, BAS, HPB};
  void assign(char* iresty, char* jresty, int cl, double gm);
  void assign(char* iresty, char* jresty, char* ifresty, char* jfresty, int cl, double gm);
  int get_index_array(char *resty, int *a);
};

#endif

// --------------------------------------------------------------------//

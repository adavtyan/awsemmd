/* ----------------------------------------------------------------------
Copyright (2010) Aram Davtyan and Garegin Papoian

Papoian's Group, University of Maryland at Collage Park
http://papoian.chem.umd.edu/

Last Update: 03/04/2011
------------------------------------------------------------------------- */

// fragment_memory.cpp

#include "fragment_memory.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


// {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"};
// {"A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"};
int fm_se_map[] = {0, 0, 4, 3, 6, 13, 7, 8, 9, 0, 11, 10, 12, 2, 0, 14, 5, 1, 15, 16, 0, 19, 17, 0, 18, 0};

// Four letter classes
// 1) SHL: Small Hydrophilic (ALA, GLY, PRO, SER THR) or (A, G, P, S, T) or {0, 7, 14, 15, 16}
// 2) AHL: Acidic Hydrophilic (ASN, ASP, GLN, GLU) or (N, D, Q, E) or {2, 3, 5, 6}
// 3) BAS: Basic (ARG HIS LYS) or (R, H, K) or {1, 8, 11}
// 4) HPB: Hydrophobic (CYS, ILE, LEU, MET, PHE, TRP, TYR, VAL) or (C, I, L, M, F, W, Y, V)  or {4, 9, 10, 12, 13, 17, 18, 19}
int four_letter_map[] = {1, 3, 2, 2, 4, 2, 2, 1, 3, 4, 4, 3, 4, 4, 1, 1, 1, 4, 4, 4};

Fragment_Memory::Fragment_Memory(int p, int pf, int l, char *fname)
{
  int i, j, nAtoms, ires, iatom, nca=0, ncb=0;
  double x, y, z, **xca, **xcb;
  char buff[101], resty[4], atomty[5];
  
  FILE * file;

  error = ERR_NONE;
  
  pos = p;
  len = l;
  mpos = p + len/2 + len%2;
  fpos = pf;
  
  rf[0] = new double*[len];
  rf[1] = new double*[len];
  xca = new double*[len];
  xcb = new double*[len];
  for (i=0;i<len;++i) {
    rf[0][i] = new double[len];
    rf[1][i] = new double[len];
    xca[i] = new double[3];
    xcb[i] = new double[3];
  }
  
  for (i=0;i<len;++i) {
    xca[i][0] = xca[i][1] = xca[i][2] = 0.0;
    xcb[i][0] = xcb[i][1] = xcb[i][2] = 0.0;
  }
  
  nca = ncb = 0;
  file = fopen(fname,"r");
  if (!file) { error = ERR_FILE; return; }
  fgets(buff, 100, file);
  fscanf(file, "%d",&nAtoms);
  for (i=0;i<nAtoms;++i) {
    fscanf(file, "%d %s %s %d %lf %lf %lf",&ires,resty,atomty,&iatom,&x,&y,&z);
    if (ires>fpos && ires<=fpos+len) {
      ires -= fpos + 1;
      x *= 10; y *= 10; z *= 10;
      if (strcmp(atomty,"CA")==0) {
        if (ires>=len || nca>=len) { error = ERR_ATOM_COUNT; return; }
        xca[ires][0] = x;
        xca[ires][1] = y;
        xca[ires][2] = z;
        nca++;
      }
      if (strcmp(atomty,"CB")==0) {
        if (ires>=len || ncb>=len) { error = ERR_ATOM_COUNT; return; }
        xcb[ires][0] = x;
        xcb[ires][1] = y;
        xcb[ires][2] = z;
        ncb++;
      }
    }
  }  
  fclose(file);
  
  if (nca!=len) { error = ERR_ATOM_COUNT; return; }
  
  for (i=0;i<len;++i) {
    for (j=0;j<len;++j) {
      rf[0][i][j] = ( i<j ? R(xca[i],xca[j]) : ( i>j ? R(xcb[i],xcb[j]) : 0.0 ) );
      rf[1][i][j] = R(xca[i],xcb[j]);
    }
  }
  
  for (i=0;i<len;++i) {
    delete [] xca[i];
    delete [] xcb[i];
  }
  delete [] xca;
  delete [] xcb;
}

Fragment_Memory::~Fragment_Memory()
{
  for (int i=0;i<len;++i) {
    delete rf[0][i];
    delete rf[1][i];
  }

  delete [] rf[0];
  delete [] rf[1];
}

double Fragment_Memory::Rf(int ires, int iatom, int jres, int jatom)
{
  if (iatom==FM_CA && jatom==FM_CA ) {
    return rf[0][min(ires,jres)][max(ires,jres)];
  } else if (iatom==FM_CB && jatom==FM_CB) {
    return rf[0][max(ires,jres)][min(ires,jres)];
  } else {
    return (iatom==FM_CA ? rf[1][ires][jres] : rf[1][jres][ires]);
  }
}

inline int Fragment_Memory::min(int a, int b)
{
  return (a<b ? a : b);
}

inline int Fragment_Memory::max(int a, int b)
{
  return (a>b ? a : b);
}

inline double Fragment_Memory::R(double *r1, double *r2)
{
  return sqrt((r1[0]-r2[0])*(r1[0]-r2[0]) + (r1[1]-r2[1])*(r1[1]-r2[1]) + (r1[2]-r2[2])*(r1[2]-r2[2]));
}

// --------------------------------------------------------------------//

Gamma_Array::Gamma_Array(char *fname)
{
  allocated = false;

  int i, iline, ns, nbuf, buf_len, ngamma, cl;
  char line[100] , buf[6][10], *st;
  char *iresty, *jresty, *ifresty, *jfresty;
  bool frag_resty_was_set = false;
  double gm;
  FILE *file;
  fpos_t pos;
  
  error = ERR_NONE;
  
  frag_resty = false;
  nres_cl = 0;

  file = fopen(fname,"r");
  if (!file) { error = ERR_FILE; return; }

  iline = 0;
  while ( fgets ( line, sizeof line, file ) != NULL ) {
    if (line[0]=='#') continue;
    
    if (iline==0) {
      ns = 0;
      st=strtok (line," \t");
      while ( st!=NULL ) {
        ns++;
        if (strcmp(st, "inf")==0) { i_sep[ns-1]=-1; break; }
        else i_sep[ns-1]=atoi(st);

        st=strtok (NULL," \t\n");
      }

      nseq_cl = ns-1;
      if (nseq_cl==0) { error = ERR_CLASS_DEF; return; }

      fgetpos (file,&pos);
    } else {
      nbuf=0;
      st=strtok (line," \t");
      while ( st!=NULL ) {
        strcpy(buf[nbuf], st);
        nbuf++;
        
        st=strtok (NULL," \t\n");
      }
      if (nbuf!=4 && nbuf!=6) { error = ERR_GAMMA; return; }
      if (!frag_resty_was_set) { frag_resty = (nbuf==4 ? false : true); frag_resty_was_set = true; }
      else if ( (nbuf==4 && frag_resty) || (nbuf==6 && !frag_resty) )  { error = ERR_GAMMA; return; }

      for (i=0;i<nbuf-2;i++) {
        buf_len = strlen(buf[i]);
        if (buf_len!=1 && buf_len!=3) { error = ERR_GAMMA; return; }

        if (strcmp(buf[i], "ALL")==0) {
          if (nres_cl<1) nres_cl = 1;
        } else if (buf_len==3) {
          if (nres_cl<4) nres_cl = 4;
        } else if (buf_len==1) {
          if (nres_cl<20) nres_cl = 20;
        }
      }
    }

    iline++;
  }

  if (nres_cl==0) { error = ERR_GAMMA; return; }
  
  if (!frag_resty) ngamma = nseq_cl*nres_cl*nres_cl;
  else ngamma = nseq_cl*nres_cl*nres_cl*nres_cl*nres_cl;
  
  gamma = new double[ngamma];
  allocated = true;
    
  fsetpos (file,&pos);
  while ( fgets ( line, sizeof line, file ) != NULL ) {
    if (line[0]=='#') continue;
    
    if (!frag_resty) {
      iresty = strtok(line," \t\n");
      jresty = strtok(NULL," \t\n");
      cl = atoi(strtok(NULL," \t\n"));
      gm = atof(strtok(NULL," \t\n"));
      assign(iresty, jresty, cl, gm);
    } else {
      iresty = strtok(line," \t\n");
      jresty = strtok(line," \t\n");
      ifresty = strtok(line," \t\n");
      jfresty = strtok(line," \t\n");
      cl = atoi(strtok(line," \t\n"));
      gm = atof(strtok(line," \t\n"));
      assign(iresty, jresty, ifresty, ifresty, cl, gm);
    }
    
    if (error!=ERR_NONE) return;
  }
  
  fclose(file);
}

Gamma_Array::~Gamma_Array()
{
  if (allocated) {
    delete [] gamma;
  }
}

double Gamma_Array::getGamma(int ires, int jres)
{
  if (nres_cl!=1) { error = ERR_CALL; return 0.0; }
  
  int dij = abs(ires-jres);
  if (dij<i_sep[0] || (i_sep[nseq_cl]!=-1 && dij>i_sep[nseq_cl])) return 0.0;
  
  int seq_cl;
  for (seq_cl=1;seq_cl<nseq_cl && dij>=i_sep[seq_cl];seq_cl++) {}
  
  return gamma[seq_cl-1];
}

double Gamma_Array::getGamma(int ires_type, int jres_type, int ires, int jres)
{
  if (frag_resty || ires_type>=20 || jres_type>=20) { error = ERR_CALL; return 0.0; }
 
  int dij = abs(ires-jres);
  if (dij<i_sep[0] || (i_sep[nseq_cl]!=-1 && dij>i_sep[nseq_cl])) return 0.0;

  int seq_cl, ires_cl, jres_cl, ig;
  
  for (seq_cl=1;seq_cl<nseq_cl && dij>=i_sep[seq_cl];seq_cl++) {}
  
  if (nres_cl==1) return gamma[seq_cl-1];
  
  if (nres_cl==20) {
    ires_cl = ires_type;
    jres_cl = jres_type;
  } else if (nres_cl==4) {
    ires_cl = four_letter_map[ires_type]-1;
    jres_cl = four_letter_map[jres_type]-1;
  }
  
  ig = (seq_cl==1 ? 0 : (seq_cl-1)*nres_cl*nres_cl) + ires_cl*nres_cl + jres_cl;
  
  return gamma[ig];
}

double Gamma_Array::getGamma(int ires_type, int jres_type, int ifres_type, int jfres_type, int ires, int jres)
{
  // Is not complit
}

int Gamma_Array::get_index_array(char *resty, int *a)
{
  int i, ifour_res_cl, n, res_code;

  if (strcmp(resty,"ALL")==0) {
    for (i=0;i<nres_cl;++i) a[i]=i;
    return nres_cl;
  
  } else if (strlen(resty)==1) {
    res_code = resty[0] - 'A';
    if (fm_se_map[res_code]==0 && resty[0]!='A') { error = ERR_GAMMA; return -1; }
    if (nres_cl<20) { error = ERR_ASSIGN; return -1; }
    
    a[0] = fm_se_map[res_code];
    return 1;

  } else if (strlen(resty)==3) {
    if (nres_cl<4) { error = ERR_ASSIGN; return -1; }
    
    if (strcmp(resty,"SHL")==0) ifour_res_cl = SHL;
    else if (strcmp(resty,"AHL")==0) ifour_res_cl = AHL;
    else if (strcmp(resty,"HPB")==0) ifour_res_cl = HPB;
    else if (strcmp(resty,"BAS")==0) ifour_res_cl = BAS;
    else { error = ERR_GAMMA; return -1; }
    
    if (nres_cl==4) {
      a[0] = ifour_res_cl-1;
      return 1;
    } else {
      n=0;
      for (i=0;i<20;++i)
        if (four_letter_map[i]==ifour_res_cl) { a[n]=i; n++; }
      return n;
    }
    
  } else { error = ERR_GAMMA; return -1; }
}

void Gamma_Array::assign(char* iresty, char* jresty, int cl, double gm)
{
  if (frag_resty) { error = ERR_ASSIGN; return; }
  if (cl>nseq_cl) { error = ERR_G_CLASS; return; }
  
  int i, j, in, jn, ig;
  int is[20], js[20];
  
  in = get_index_array(iresty, is);
  jn = get_index_array(jresty, js);
  if (error!=ERR_NONE) return;
  
  for (i=0;i<in;++i) {
    for (j=0;j<jn;++j) {    
      ig = (cl-1)*nres_cl*nres_cl + is[i]*nres_cl + js[j];
      gamma[ig] = gm;
      if (is[i]!=js[j]) {
        ig = (cl-1)*nres_cl*nres_cl + js[j]*nres_cl + is[i];
        gamma[ig] = gm;
      }
    }
  }
}

void Gamma_Array::assign(char* iresty, char* jresty, char* ifresty, char* jfresty, int cl, double gm)
{
  if (!frag_resty) { error = ERR_ASSIGN; return; }
  if (cl>nseq_cl) { error = ERR_ASSIGN; return; }
  
  // Is not complite
}

// --------------------------------------------------------------------//

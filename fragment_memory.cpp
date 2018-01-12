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

#include <iostream>
using namespace std;


// {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"};
// {"A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"};
int fm_se_map[] = {0, 0, 4, 3, 6, 13, 7, 8, 9, 0, 11, 10, 12, 2, 0, 14, 5, 1, 15, 16, 0, 19, 17, 0, 18, 0};

// Four letter classes
// 1) SHL: Small Hydrophilic (ALA, GLY, PRO, SER THR) or (A, G, P, S, T) or {0, 7, 14, 15, 16}
// 2) AHL: Acidic Hydrophilic (ASN, ASP, GLN, GLU) or (N, D, Q, E) or {2, 3, 5, 6}
// 3) BAS: Basic (ARG HIS LYS) or (R, H, K) or {1, 8, 11}
// 4) HPB: Hydrophobic (CYS, ILE, LEU, MET, PHE, TRP, TYR, VAL) or (C, I, L, M, F, W, Y, V)  or {4, 9, 10, 12, 13, 17, 18, 19}
int four_letter_map[] = {1, 3, 2, 2, 4, 2, 2, 1, 3, 4, 4, 3, 4, 4, 1, 1, 1, 4, 4, 4};

Fragment_Memory::Fragment_Memory(int p, int pf, int l, double w, char *fname, bool vec_fm_flag)
{
  int i, j, nAtoms, ires, iatom, nca=0, ncb=0;
  double x, y, z, **xca, **xcb;
  char buff[101], resty[4], atomty[5];

  FILE * file;

  FILE * dout;
  dout = fopen("debug.info","w");

  error = ERR_NONE;

  pos = p;
  len = l;
  mpos = p + len/2 + len%2;
  fpos = pf;
  weight = w;
  vfm_flag = vec_fm_flag;

  se = new char[len];
  rf[0] = new double*[len];
  rf[1] = new double*[len];
  if (vfm_flag) vmf = new double*[len];
  xca = new double*[len];
  xcb = new double*[len];
  for (i=0;i<len;++i) {
    rf[0][i] = new double[len];
    rf[1][i] = new double[len];
    if (vfm_flag) vmf[i] = new double[len];
    xca[i] = new double[3];
    xcb[i] = new double[3];
  }

  for (i=0;i<len;++i) {
    xca[i][0] = xca[i][1] = xca[i][2] = 0.0;
    xcb[i][0] = xcb[i][1] = xcb[i][2] = 0.0;
  }

  nca = ncb = 0;
//  fprintf(dout, "%s\n", fname);
  file = fopen(fname,"r");
  if (!file) { error = ERR_FILE; return; }
  fgets(buff, 100, file);
  fscanf(file, "%d",&nAtoms);

  // -------------------------------------------------------------------------
  // Commented by HaoWu
  // This will assign coordinates of wrong chain when reading multi-chain gro file
  // Not applicable with AMH-Go model

  // for (i=0;i<nAtoms;++i) {
  //   fscanf(file, "%d %s %s %d %lf %lf %lf",&ires,resty,atomty,&iatom,&x,&y,&z);
  //   if (ires>fpos && ires<=fpos+len) {
  //     ires -= fpos + 1;
  //     x *= 10; y *= 10; z *= 10;
  //     if (strcmp(atomty,"CA")==0) {
  //       if (ires>=len || nca>=len) { error = ERR_ATOM_COUNT; return; }
  //       se[ires] = ThreeLetterToOne(resty);
  //       if (se[ires]=='-') { error = ERR_RES; return; }
  //       xca[ires][0] = x;
  //       xca[ires][1] = y;
  //       xca[ires][2] = z;
  //       nca++;
  //     }
  //     if (strcmp(atomty,"CB")==0) {
  //       if (ires>=len || ncb>=len) { error = ERR_ATOM_COUNT; return; }
  //       xcb[ires][0] = x;
  //       xcb[ires][1] = y;
  //       xcb[ires][2] = z;
  //       ncb++;
  //     }
  //   }
  // }
  // -------------------------------------------------------------------------

  // HaoWu copy from BinZhang's version, reads multi-chain gro file correctly
  int count=-1;
  int old_ires=-9999;
  for (i=0;i<nAtoms;++i) {
    fscanf(file, "%d %s %s %d %lf %lf %lf",&ires,resty,atomty,&iatom,&x,&y,&z);
    //if (ires>fpos && ires<=fpos+len) {
    if(count>=len){fprintf(stderr, "nca=%d, count=%d, ires=%d, fpos=%d\n", nca, count, ires, fpos); break;} //faster than reading the whole file
    if (ires>fpos) {
      if(ires!=old_ires){
    ++count;
    old_ires=ires;
    if(count>=len){
        //fprintf(stderr, "nca=%d, count=%d, ires=%d, fpos=%d\n", nca, count, ires, fpos);
        break; //faster than reading the whole file
    }
      }
      //ires -= fpos + 1; //use count instead of ires
      x *= 10; y *= 10; z *= 10;
      if (strcmp(atomty,"CA")==0) {
        if (count>=len || nca>=len) { error = ERR_ATOM_COUNT; fprintf(stderr, "CA atoms are problematic!\n"); return; }
        se[count] = ThreeLetterToOne(resty);
        if (se[count]=='-') { error = ERR_RES; return; }
        xca[count][0] = x;
        xca[count][1] = y;
        xca[count][2] = z;
        nca++;
      }
      if (strcmp(atomty,"CB")==0) {
        if (count>=len || ncb>=len) { error = ERR_ATOM_COUNT; fprintf(stderr, "CB atoms are problematic!\n"); return; }
        xcb[count][0] = x;
        xcb[count][1] = y;
        xcb[count][2] = z;
        ncb++;
      }
    }
  }

  fclose(file);

  if (nca!=len) {
      fprintf(dout, "File: %s\n", fname);
      fprintf(dout, "nca=%d\n", nca);
      fprintf(dout, "len=%d\n", len);
  }

  if (nca!=len) { error = ERR_ATOM_COUNT; return; }

  for (i=0;i<len;++i) {

    for (j=0;j<len;++j) {

      rf[0][i][j] = ( i<j ? R(xca[i],xca[j]) : ( i>j ? R(xcb[i],xcb[j]) : 0.0 ) );
      rf[1][i][j] = R(xca[i],xcb[j]);

      if (vfm_flag)
          vmf[i][j] = VM(xca[i],xcb[i],xca[j],xcb[j]); // Multiplication (normalized) of CA-CB vectors between residue i and j
    }
  }

  for (i=0;i<len;++i) {
    delete [] xca[i];
    delete [] xcb[i];
  }
  delete [] xca;
  delete [] xcb;

  fclose(dout);

}

Fragment_Memory::~Fragment_Memory()
{
  for (int i=0;i<len;++i) {
    delete rf[0][i];
    delete rf[1][i];
  }

  delete [] rf[0];
  delete [] rf[1];
  delete [] se;
}

double Fragment_Memory::Rf(int ires, int iatom, int jres, int jatom)
{
  ires -= pos;
  jres -= pos;
  if (ires<0 || ires>=len || jres<0 || jres>=len) { error = ERR_CALL; return 0.0; }
  if (iatom==FM_CA && jatom==FM_CA ) {
    return rf[0][min(ires,jres)][max(ires,jres)];
  } else if (iatom==FM_CB && jatom==FM_CB) {
    return rf[0][max(ires,jres)][min(ires,jres)];
  } else {
    return (iatom==FM_CA ? rf[1][ires][jres] : rf[1][jres][ires]);
  }
}

double Fragment_Memory::VMf(int ires, int jres)
{
  ires -= pos;
  jres -= pos;
  if (!vfm_flag || ires<0 || ires>=len || jres<0 || jres>=len) { error = ERR_CALL; return 0.0; }
  if (se[ires]=='G' || se[jres]=='G') { error = ERR_VFM_GLY; return 0.0; }

  return vmf[ires][jres];
}

char Fragment_Memory::getSe(int resno)
{
  return se[resno - pos];
}

int Fragment_Memory::resType(int resno)
{
  return fm_se_map[se[resno - pos]-'A'];
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

// Calculates theta angle between vectors (r11, r12) and (r21, r22)
inline double Fragment_Memory::VM(double *r11, double *r12, double *r21, double *r22)
{
  return acos(((r12[0]-r11[0])*(r22[0]-r21[0]) + (r12[1]-r11[1])*(r22[1]-r21[1]) + (r12[2]-r11[2])*(r22[2]-r21[2]))/(R(r11,r12)*R(r21,r22)));
}

char Fragment_Memory::ThreeLetterToOne(char *tl_resty)
{
  if (strlen(tl_resty)==3) {
    if (strcmp(tl_resty,"ALA")==0) return 'A';
    else if (strcmp(tl_resty,"ARG")==0) return 'R';
    else if (strcmp(tl_resty,"ASN")==0) return 'N';
    else if (strcmp(tl_resty,"ASP")==0) return 'D';
    else if (strcmp(tl_resty,"CYS")==0) return 'C';
    else if (strcmp(tl_resty,"GLN")==0) return 'Q';
    else if (strcmp(tl_resty,"GLU")==0) return 'E';
    else if (strcmp(tl_resty,"GLY")==0) return 'G';
    else if (strcmp(tl_resty,"HIS")==0) return 'H';
    else if (strcmp(tl_resty,"ILE")==0) return 'I';
    else if (strcmp(tl_resty,"LEU")==0) return 'L';
    else if (strcmp(tl_resty,"LYS")==0) return 'K';
    else if (strcmp(tl_resty,"MET")==0) return 'M';
    else if (strcmp(tl_resty,"PHE")==0) return 'F';
    else if (strcmp(tl_resty,"PRO")==0) return 'P';
    else if (strcmp(tl_resty,"SER")==0) return 'S';
    else if (strcmp(tl_resty,"THR")==0) return 'T';
    else if (strcmp(tl_resty,"TRP")==0) return 'W';
    else if (strcmp(tl_resty,"TYR")==0) return 'Y';
    else if (strcmp(tl_resty,"VAL")==0) return 'V';
  }

  return '?';
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
    if (isEmptyString(line)) continue;

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
      st=strtok (line," \t\n");
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
    if (isEmptyString(line)) continue;

    if (!frag_resty) {
      iresty = strtok(line," \t\n");
      jresty = strtok(NULL," \t\n");
      cl = atoi(strtok(NULL," \t\n"));
      gm = atof(strtok(NULL," \t\n"));
      assign(iresty, jresty, cl, gm);
    } else {
      iresty = strtok(line," \t\n");
      jresty = strtok(NULL," \t\n");
      ifresty = strtok(NULL," \t\n");
      jfresty = strtok(NULL," \t\n");
      cl = atoi(strtok(NULL," \t\n"));
      gm = atof(strtok(NULL," \t\n"));
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
  if (!frag_resty || ires_type>=20 || jres_type>=20 || ifres_type>=20 || jfres_type>=20) { error = ERR_CALL; return 0.0; }

  int dij = abs(ires-jres);
  if (dij<i_sep[0] || (i_sep[nseq_cl]!=-1 && dij>i_sep[nseq_cl])) return 0.0;

  int seq_cl, ires_cl, jres_cl, ifres_cl, jfres_cl, ig;

  for (seq_cl=1;seq_cl<nseq_cl && dij>=i_sep[seq_cl];seq_cl++) {}

  if (nres_cl==1) return gamma[seq_cl-1];

  if (nres_cl==4) {
    ires_cl = four_letter_map[ires_type]-1;
    jres_cl = four_letter_map[jres_type]-1;
    ifres_cl = four_letter_map[ifres_type]-1;
    jfres_cl = four_letter_map[jfres_type]-1;
  } else if (nres_cl==20) {
    ires_cl = ires_type;
    jres_cl = jres_type;
    ifres_cl = ifres_type;
    jfres_cl = jfres_type;
  }

  ig = (seq_cl==1 ? 0 : (seq_cl-1)*nres_cl*nres_cl*nres_cl*nres_cl) + ires_cl*nres_cl*nres_cl*nres_cl + jres_cl*nres_cl*nres_cl + ifres_cl*nres_cl + jfres_cl;

  return gamma[ig];
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
      // gamma[iT][jT]
      ig = (cl-1)*nres_cl*nres_cl + is[i]*nres_cl + js[j];
      gamma[ig] = gm;

      // gamma[jT][iT]
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

  int i, j, k, l, in, jn, kn, ln, ig;
  int is[20], js[20], ks[20], ls[20];

  in = get_index_array(iresty, is);
  jn = get_index_array(jresty, js);
  kn = get_index_array(ifresty, ks);
  ln = get_index_array(jfresty, ls);
  if (error!=ERR_NONE) return;

  for (i=0;i<in;++i) {
    for (j=0;j<jn;++j) {
      for (k=0;k<kn;++k) {
        for (l=0;l<ln;++l) {
          // gamma[iT][jT][iF][jF]
          ig = (cl-1)*nres_cl*nres_cl*nres_cl*nres_cl + is[i]*nres_cl*nres_cl*nres_cl + js[j]*nres_cl*nres_cl + ks[k]*nres_cl + ls[l];
          gamma[ig] = gm;

          // gamma[jT][iT][iF][jF]
          if (is[i]!=js[j]) {
            ig = (cl-1)*nres_cl*nres_cl*nres_cl*nres_cl + js[j]*nres_cl*nres_cl*nres_cl + is[i]*nres_cl*nres_cl + ks[k]*nres_cl + ls[l];
            gamma[ig] = gm;
          }

          // gamma[iT][jT][jF][iF]
          if (ks[k]!=ls[l]) {
            ig = (cl-1)*nres_cl*nres_cl*nres_cl*nres_cl + is[i]*nres_cl*nres_cl*nres_cl + js[j]*nres_cl*nres_cl + ls[l]*nres_cl + ks[k];
            gamma[ig] = gm;
          }

          // gamma[jT][iT][jF][iF]
          if (is[i]!=js[j] && ks[k]!=ls[l]) {
            ig = (cl-1)*nres_cl*nres_cl*nres_cl*nres_cl + js[j]*nres_cl*nres_cl*nres_cl + is[i]*nres_cl*nres_cl + ls[l]*nres_cl + ks[k];
            gamma[ig] = gm;
          }
        }
      }
    }
  }
}

bool Gamma_Array::isEmptyString(char *str)
{
  int len = strlen(str);

  if (len==0) return true;

  for (int i=0;i<len;++i) {
    if (str[i]!=' ' && str[i]!='\t' && str[i]!='\n') return false;
  }

  return true;
}

int Gamma_Array::minSep()
{
  return i_sep[0];
}

int Gamma_Array::maxSep()
{
  return i_sep[nseq_cl];
}

// --------------------------------------------------------------------//

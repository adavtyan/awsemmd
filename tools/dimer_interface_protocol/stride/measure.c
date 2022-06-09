#include "stride.h"

void Measure(CHAIN **Chain, int NChain, int El, COMMAND *Cmd, FILE *Out)
{

  QUALITY QualDssp_H, QualDssp_E, Qual_H, Qual_E;   
  int i, Cn, Flag1, Flag2;
  float Q2_, Q2_Dssp, O_, O_Dssp;
  float Content, AlphaCont, BetaCont, PerCor, PerCorDssp;
  int HelAlp, HelPI, Hel310, Sheet, Turn;
  int YYN_H, NYY_H, YNN_H, NNY_H, YYN_E, NYY_E, YNN_E, NNY_E;
  int Total_H, Better_H, Worse_H, Total_E, Better_E, Worse_E;
  char Tmp[3][MAX_RES], *Asn, *PdbAsn, *DsspAsn;
  CHAIN *c;


  for( Cn=0; Cn<NChain; Cn++ ) {

    c = Chain[Cn];

    if( !c->Valid )
      continue;

    Asn     = (char *)ckalloc(c->NRes*sizeof(char));
    PdbAsn  = (char *)ckalloc(c->NRes*sizeof(char));
    DsspAsn = (char *)ckalloc(c->NRes*sizeof(char));
    
    ExtractAsn(Chain,Cn,Asn);
    ExtractPdbAsn(Chain,Cn,PdbAsn);
    ExtractDsspAsn(Chain,Cn,DsspAsn);

    for( i=0; i<c->NRes; i++ ) {
      if( Asn[i] != 'H' && Asn[i] != 'E' && Asn[i] != 'B' )
	Asn[i] = 'C';
      if( PdbAsn[i] != 'H' && PdbAsn[i] != 'E' )
	PdbAsn[i] = 'C';
      if( DsspAsn[i] != 'H' && DsspAsn[i] != 'E' && DsspAsn[i] != 'B' )
	DsspAsn[i] = 'C';
      if( Asn[i] == 'B' )
	Asn[i] = 'E';
      if( DsspAsn[i] == 'B' )
	DsspAsn[i] = 'E';
    }

    CorrectAsnDouble(Asn,DsspAsn,PdbAsn,c->NRes,'E','C');


    if( El ) {
      FullElement(Asn,PdbAsn,DsspAsn,c->NRes,'H',4,'z',&YYN_H,&NYY_H,&YNN_H,&NNY_H);
      FullElement(Asn,PdbAsn,DsspAsn,c->NRes,'E',1,'y',&YYN_E,&NYY_E,&YNN_E,&NNY_E);
      fprintf(Out,"DIFF_EL ");
    }
    else
      fprintf(Out,"DIFF ");
    
    Content = SecStrContent(c,&HelAlp,&HelPI,&Hel310,&Sheet,&Turn);
    AlphaCont = (float)HelAlp/(float)c->NRes;
    BetaCont  = (float)Sheet/(float)c->NRes;
    
    Presnell(PdbAsn,c->NRes,DsspAsn,c->NRes,'E',0.5,&Q2_Dssp,&O_Dssp);
    Presnell(PdbAsn,c->NRes,Asn,c->NRes,'E',0.5,&Q2_,&O_);
    
    Difference(DsspAsn,PdbAsn,c->NRes,'H',&QualDssp_H);
    Difference(DsspAsn,PdbAsn,c->NRes,'E',&QualDssp_E);
    Difference(Asn,PdbAsn,c->NRes,'H',&Qual_H);
    Difference(Asn,PdbAsn,c->NRes,'E',&Qual_E);
    
    PerCor     = PercentCorrect(Asn,PdbAsn,c->NRes);
    PerCorDssp = PercentCorrect(DsspAsn,PdbAsn,c->NRes);
    
    Total_H = CompareElements(Asn,PdbAsn,DsspAsn,c->NRes,'H',&Better_H,&Worse_H);
    Total_E = CompareElements(Asn,PdbAsn,DsspAsn,c->NRes,'E',&Better_E,&Worse_E);
    
    if( (Flag1=AssessCorr(&Qual_H)) ) 
      fprintf(Out,"%6.4f ",Qual_H.Corr);
    else fprintf(Out,"  None ");
    fprintf(Out,"%6.4f ",Qual_H.Perc);
    
    if( (Flag2=AssessCorr(&QualDssp_H)) ) 
      fprintf(Out,"%6.4f ",QualDssp_H.Corr);
    else fprintf(Out,"  None ");
    fprintf(Out,"%6.4f ",QualDssp_H.Perc);
    
    if( Flag1 && Flag2 ) 
      fprintf(Out,"%7.4f ",Qual_H.Corr-QualDssp_H.Corr);
    else fprintf(Out,"   None ");
    fprintf(Out,"%7.4f | ",Qual_H.Perc-QualDssp_H.Perc);
    
    
    if( (Flag1=AssessCorr(&Qual_E)) ) 
      fprintf(Out,"%6.4f ",Qual_E.Corr);
    else fprintf(Out,"  None ");
    fprintf(Out,"%6.4f ",Qual_E.Perc);
    
    if( (Flag2=AssessCorr(&QualDssp_E)) ) 
      fprintf(Out,"%6.4f ",QualDssp_E.Corr);
    else fprintf(Out,"  None ");
    fprintf(Out,"%6.4f ",QualDssp_E.Perc);
    
    if( Flag1 && Flag2 ) 
      fprintf(Out,"%7.4f ",Qual_E.Corr-QualDssp_E.Corr);
    else fprintf(Out,"   None ");
    fprintf(Out,"%7.4f | ",Qual_E.Perc-QualDssp_E.Perc);
    
    fprintf(Out,"%6.4f %6.4f | %6.4f %6.4f | ",Q2_,Q2_Dssp,O_,O_Dssp);
    
    fprintf(Out,"%4d %4d %4d %4d | %4d %4d %4d %4d |%4d %4d %4d %4d |%4d %4d %4d %4d | %s%c %4d %4.2f | %5.3f %5.3f %5.3f",
	    Qual_H.TP,Qual_H.TN,Qual_H.FP,Qual_H.FN,
	    QualDssp_H.TP,QualDssp_H.TN,QualDssp_H.FP,QualDssp_H.FN,
	    Qual_E.TP,Qual_E.TN,Qual_E.FP,Qual_E.FN,
	    QualDssp_E.TP,QualDssp_E.TN,QualDssp_E.FP,QualDssp_E.FN,
	    c->File,c->Id,c->NRes,c->Resolution,Content,AlphaCont,BetaCont);
    
    if( El ) {
      fprintf(Out," FullH: YYN %2d NYY %2d YNN %2d NNY %2d ",YYN_H,NYY_H,YNN_H,NNY_H);
      fprintf(Out," FullE: YYN %2d NYY %2d YNN %2d NNY %2d ",YYN_E,NYY_E,YNN_E,NNY_E);
    }
    
    memset(Tmp[0],'C',c->NRes);
    memset(Tmp[1],'C',c->NRes);
    memset(Tmp[2],'C',c->NRes);
    
    for( i=0; i<c->NRes; i++ ) {
      if( Asn[i] == 'H' ) 
	Tmp[0][i] = Asn[i];
      if( PdbAsn[i] == 'H' )
	Tmp[1][i] = PdbAsn[i];
      if( DsspAsn[i] == 'H' )
	Tmp[2][i] = DsspAsn[i];
    }
    ExcludeObvious(Tmp[0],Tmp[2],Tmp[1],c->NRes);
    Difference(Tmp[0],Tmp[1],c->NRes,'H',&Qual_H);
    fprintf(Out," | ResH: YYN %3d NNY %3d YNN %3d NYY %3d",
	    Qual_H.TP,Qual_H.TN,Qual_H.FP,Qual_H.FN);
    
    memset(Tmp[0],'C',c->NRes);
    memset(Tmp[1],'C',c->NRes);
    memset(Tmp[2],'C',c->NRes);
    
    for( i=0; i<c->NRes; i++ ) {
      if( Asn[i] == 'E' ) 
	Tmp[0][i] = Asn[i];
      if( PdbAsn[i] == 'E' )
	Tmp[1][i] = PdbAsn[i];
      if( DsspAsn[i] == 'E' )
	Tmp[2][i] = DsspAsn[i];
    }
    ExcludeObvious(Tmp[0],Tmp[2],Tmp[1],c->NRes);
    Difference(Tmp[0],Tmp[1],c->NRes,'E',&Qual_E);
    fprintf(Out," | ResE: YYN %3d NNY %3d YNN %3d NYY %3d",
	    Qual_E.TP,Qual_E.TN,Qual_E.FP,Qual_E.FN);
    
    fprintf(Out," | ToH: %2d BtH: %2d WsH: %2d ",Total_H,Better_H,Worse_H);
    fprintf(Out," | ToE: %2d BtE: %2d WsE: %2d ",Total_E,Better_E,Worse_E);
    fprintf(Out,"PerCor %7.4f PerCorDssp %7.4f\n",PerCor,PerCorDssp);
    
    fprintf(Out,"\n");
    
/*   int Bound[MAX_ASSIGN][2], NElem; */
/*     if( !El ) {
 *       NElem = Boundaries(Asn,c->NRes,'H',Bound);
 *       for( i=0; i<NElem; i++ )
 * 	fprintf(Out,"HELIXSTRIDE %4d residues long in %s%c\n",
 * 		Bound[i][1]-Bound[i][0]+1,c->File,c->Id);
 * 
 *       NElem = Boundaries(Asn,c->NRes,'E',Bound);
 *       for( i=0; i<NElem; i++ )
 * 	fprintf(Out,"STRANDSTRIDE %4d residues long in %s%c\n",
 * 		Bound[i][1]-Bound[i][0]+1,c->File,c->Id);
 * 
 *       NElem = Boundaries(PdbAsn,c->NRes,'H',Bound);
 *       for( i=0; i<NElem; i++ )
 * 	fprintf(Out,"HELIXPDB %4d residues long in %s%c\n",
 * 		Bound[i][1]-Bound[i][0]+1,c->File,c->Id);
 * 
 *       NElem = Boundaries(PdbAsn,c->NRes,'E',Bound);
 *       for( i=0; i<NElem; i++ )
 * 	fprintf(Out,"STRANDPDB %4d residues long in %s%c\n",
 * 		Bound[i][1]-Bound[i][0]+1,c->File,c->Id);
 * 
 *       NElem = Boundaries(DsspAsn,c->NRes,'H',Bound);
 *       for( i=0; i<NElem; i++ )
 * 	fprintf(Out,"HELIXSDSSP %4d residues long in %s%c\n",
 * 		Bound[i][1]-Bound[i][0]+1,c->File,c->Id);
 * 
 *       NElem = Boundaries(DsspAsn,c->NRes,'E',Bound);
 *       for( i=0; i<NElem; i++ )
 * 	fprintf(Out,"STRANDDSSP %4d residues long in %s%c\n",
 * 		Bound[i][1]-Bound[i][0]+1,c->File,c->Id);
 *     }
 */
    free(Asn);
    free(PdbAsn);
    free(DsspAsn);
  }
}


int Presnell(char *Asn1, int L1, char *Asn2, int L2, char SecStr, float Threshold, 
	     float *Q2, float *O)
{
  int Boundaries(char *Asn, int L, char SecondStr, int (*Bound)[2]);
  int Bound1[MAX_ASSIGN][2], Bound2[MAX_ASSIGN][2], Length1[MAX_ASSIGN], 
      Length2[MAX_ASSIGN], NSeg1, NSeg2;
  int Overlap, MaxOverlap, TP=0, FP=0, FN=0;
  register int i, j;

  NSeg1 = Boundaries(Asn1,L1,SecStr,Bound1);
  NSeg2 = Boundaries(Asn2,L2,SecStr,Bound2);
  

  for(i=0; i<NSeg1; i++)
    Length1[i] = Bound1[i][1]-Bound1[i][0]+1;

  for(j=0; j<NSeg2; j++)
    Length2[j] = Bound2[j][1]-Bound2[j][0]+1;
									    
  for(i=0; i<NSeg1; i++) {
    MaxOverlap = 0;
    for(j=0; j<NSeg2; j++) {
      Overlap = Minimum(Bound1[i][1],Bound2[j][1])-Maximum(Bound1[i][0],Bound2[j][0])+1;
      if( Overlap > MaxOverlap ) MaxOverlap = Overlap;
    }
    if( (float)MaxOverlap/(float)Length1[i] >= Threshold ) 
      TP++;
    else
      FN++;
  }

  for(i=0; i<NSeg2; i++) {
    MaxOverlap = 0;
    for(j=0; j<NSeg1; j++) {
      Overlap = Minimum(Bound2[i][1],Bound1[j][1])-Maximum(Bound2[i][0],Bound1[j][0])+1;
      if( Overlap > MaxOverlap ) MaxOverlap = Overlap;
    }
    if( (float)MaxOverlap/(float)Length2[i] < Threshold ) 
      FP++;
  }

  if( TP+FN != 0 ) {
    *Q2 = (float)TP/((float)TP+(float)FN);
    *O = (float)FP/((float)TP+(float)FN);
  }
  else {
    *Q2 = -1.0;
    *O = -1.0;
  }

  return(1);
}			      


int Sov(char *Asn1, int L1, char *Asn2, int L2, char SecStr, float Threshold, float *Q2)
{
  int Bound1[MAX_ASSIGN][2], Bound2[MAX_ASSIGN][2], Length1[MAX_ASSIGN], 
      Length2[MAX_ASSIGN], NSeg1, NSeg2;
  int Overlap, MaxOverlap, TP=0, FN=0;
  register int i, j;

  NSeg1 = Boundaries(Asn1,L1,SecStr,Bound1);
  NSeg2 = Boundaries(Asn2,L2,SecStr,Bound2);
  

  for(i=0; i<NSeg1; i++)
    Length1[i] = Bound1[i][1]-Bound1[i][0]+1;

  for(j=0; j<NSeg2; j++)
    Length2[j] = Bound2[j][1]-Bound2[j][0]+1;
									    
  for(i=0; i<NSeg1; i++) {
    MaxOverlap = 0;
    for(j=0; j<NSeg2; j++) {
      
      if( Bound1[i][0] > Bound2[j][1] || Bound1[i][1] < Bound2[j][0] ) /*No overlap */
	Overlap = 0;
      else if ( Bound2[j][0] >= Bound1[i][0] && Bound2[j][1] <= Bound1[i][1] ) /* j inside i */
	Overlap = Length2[j];
      else if ( Bound1[i][0] <= Bound2[j][0] ) /* j shifted to the right with respect to i */ 
	Overlap = Minimum(Bound2[j][1]-Bound2[j][0],Bound1[i][1]-Bound2[j][0])+1;
      else if( Bound1[i][0] >= Bound2[j][0] ) /* i shifted to the right with respect to j */
	Overlap = Minimum(Bound2[j][1]-Bound2[j][0],Bound2[j][1]-Bound1[i][0])+1;
      else
	Overlap = Length1[i]; /* i inside j */
      if( Overlap > MaxOverlap ) MaxOverlap = Overlap;
    }
    if( (float)MaxOverlap/(float)Length1[i] >= Threshold ) 
      TP++;
    else
      FN++;
  }
  if( TP+FN != 0 )
    *Q2 = (float)TP/((float)TP+(float)FN);
  else
    *Q2 = -1.0;
/*  O  = (float)TP/((float)TP+(float)FN);*/
  return(1);
}			      



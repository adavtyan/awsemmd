#include "stride.h"


void Sheet(CHAIN **Chain, int Cn1, int Cn2, HBOND **HBond, COMMAND *Cmd, float **PhiPsiMap)
{
  PATTERN **PatN, **PatP;
  RESIDUE *Res1, *Res3, *Res2, *Res4, *ResA, *ResB, *Res1m1, *Res3p1;
  int R1, R3, R2, R4, RA, RB, PatCntN = 0, PatCntP = 0, Beg;
  char *AntiPar1, *Par1, *AntiPar2, *Par2;
  register int i;

  PatN = (PATTERN **)ckalloc(MAXHYDRBOND*sizeof(PATTERN *));
  PatP = (PATTERN **)ckalloc(MAXHYDRBOND*sizeof(PATTERN *));

  AntiPar1  = (char *)ckalloc(Chain[Cn1]->NRes*sizeof(char)); /* Antiparallel strands */
  Par1      = (char *)ckalloc(Chain[Cn1]->NRes*sizeof(char)); /* Parallel strands */
  AntiPar2  = (char *)ckalloc(Chain[Cn2]->NRes*sizeof(char)); /* Antiparallel strands */
  Par2      = (char *)ckalloc(Chain[Cn2]->NRes*sizeof(char)); /* Parallel strands */

  for( i=0; i<Chain[Cn1]->NRes; i++ ) { 
    AntiPar1[i] = 'C'; 
    Par1[i] = 'C'; 
  }

  for( i=0; i<Chain[Cn2]->NRes; i++ ) { 
    AntiPar2[i] = 'C'; 
    Par2[i] = 'C'; 
  }

  for( R1=0; R1<Chain[Cn1]->NRes; R1++ ) {

    Res1   = Chain[Cn1]->Rsd[R1];

    if( (!Res1->Inv->NBondDnr && !Res1->Inv->NBondAcc) ||
        ((Cn1 != Cn2) && !Res1->Inv->InterchainHBonds) )
      continue;

    RA     = R1+1;
    R2     = R1+2;
    Res1m1 = Chain[Cn1]->Rsd[R1-1];
    ResA   = Chain[Cn1]->Rsd[RA];
    Res2   = Chain[Cn1]->Rsd[R2]; 

    if( R2 >= Chain[Cn1]->NRes || 
        Res1->Prop->PhiZn == ERR || Res1->Prop->PsiZn == ERR ||
        Res2->Prop->PhiZn == ERR || Res2->Prop->PsiZn == ERR ||
        ResA->Prop->PhiZn == ERR || ResA->Prop->PsiZn == ERR ) 
      continue;

    if( Cn1 != Cn2 ) 
      Beg = 0;
    else
      Beg = R1+1;

    for( R3=Beg; R3<Chain[Cn2]->NRes; R3++ ) {

      /* Process anti-parallel strands */

      Res3   = Chain[Cn2]->Rsd[R3];

      if( (!Res3->Inv->NBondAcc && !Res3->Inv->NBondDnr ) || 
	  ((Cn1 != Cn2) && !Res3->Inv->InterchainHBonds) ) 
	continue;

      RB     = R3-1;
      R4     = R3-2;
      Res3p1 = Chain[Cn2]->Rsd[R3+1];
      ResB   = Chain[Cn2]->Rsd[RB];
      Res4   = Chain[Cn2]->Rsd[R4];
      
      if( Cn1 != Cn2 || R3 - R1 >= 3  )
	Link(HBond,Chain,Cn1,Cn2,Res1,Res3,Res3,Res1,Res1,Res3,
	     PhiPsiMap,PatN,&PatCntN,"1331",Cmd->Treshold_E1,Cmd,0);
      
      if( R2 < Chain[Cn1]->NRes && ((Cn1 != Cn2 && R4 >= 0) || R4-R2 >=2 ) )
	Link(HBond,Chain,Cn2,Cn1,Res3,Res1,Res2,Res4,ResB,ResA,
	     PhiPsiMap,PatN,&PatCntN,"3124",Cmd->Treshold_E1,Cmd,0);
      
      if( ((Cn1 != Cn2 && RB >= 0 ) || RB-R1 > 4) && 
	 ( RA >= Chain[Cn1]->NRes || (Cn1 == Cn2 && R3-RA <= 4 ) ||
	  !Link(HBond,Chain,Cn1,Cn2,Res1,Res3,Res3,ResA,NULL,Res3,
		PhiPsiMap,PatN,&PatCntN,"133A",Cmd->Treshold_E1,Cmd,1))
	 &&
	 ( R1-1 < 0 ||
	  !Link(HBond,Chain,Cn1,Cn2,Res1m1,ResB,ResB,Res1,NULL,ResB,
		PhiPsiMap,PatN,&PatCntN,"1-BB1",Cmd->Treshold_E1,Cmd,1)))
	Link(HBond,Chain,Cn1,Cn2,Res1,Res3,ResB,Res1,Res1,NULL,
	     PhiPsiMap,PatN,&PatCntN,"13B1",Cmd->Treshold_E1,Cmd,0);
      
      if( (RA < Chain[Cn1]->NRes && (Cn1 != Cn2 || R3-RA > 4)) &&
	 ( (Cn1 == Cn2 && RB-R1 <= 4 ) || (Cn1 != Cn2 && RB < 0 ) || 
	  !Link(HBond,Chain,Cn1,Cn2,Res1,Res3,ResB,Res1,Res1,NULL,
	        PhiPsiMap,PatN,&PatCntN,"13B1",Cmd->Treshold_E1,Cmd,1))
	 &&
	 ( R3+1 >= Chain[Cn2]->NRes ||
	  !Link(HBond,Chain,Cn1,Cn2,ResA,Res3p1,Res3,ResA,ResA,NULL,
		PhiPsiMap,PatN,&PatCntN,"A3+3A",Cmd->Treshold_E1,Cmd,1)))
	Link(HBond,Chain,Cn1,Cn2,Res1,Res3,Res3,ResA,NULL,Res3,
	     PhiPsiMap,PatN,&PatCntN,"133A",Cmd->Treshold_E1,Cmd,0);
      
      /* Process parallel strands */

      R4 = R3+2; 
      RB = R3+1;
      ResB   = Chain[Cn2]->Rsd[RB];
      Res4   = Chain[Cn2]->Rsd[R4];

      if( (Cn1 == Cn2 && abs(R3-R1) <= 3) || R4 >= Chain[Cn2]->NRes ) continue;
      
      if( R2 < Chain[Cn1]->NRes && (Cn1 != Cn2 || abs(R2-R3) > 3) )
	Link(HBond,Chain,Cn2,Cn1,Res3,Res1,Res2,Res3,Res3,ResA,
	     PhiPsiMap,PatP,&PatCntP,"3123",Cmd->Treshold_E2,Cmd,0);
      
      if( R4 < Chain[Cn2]->NRes && (Cn1 != Cn2 || abs(R4-R1) > 3) )
	Link(HBond,Chain,Cn1,Cn2,Res1,Res3,Res4,Res1,Res1,ResB,
	     PhiPsiMap,PatP,&PatCntP,"1341",Cmd->Treshold_E2,Cmd,0);
    }
  }
    
  FilterAntiPar(PatN,PatCntN);
  FilterPar(PatP,PatCntP);

  MergePatternsAntiPar(PatN,PatCntN);
  MergePatternsPar(PatP,PatCntP);

  if( Cmd->Info )  {
    PrintPatterns(PatN,PatCntN,Chain,Cn1,Cn2);
    PrintPatterns(PatP,PatCntP,Chain,Cn1,Cn2);
  }

  FillAsnAntiPar(AntiPar1,AntiPar2,Chain,Cn1,Cn2,PatN,PatCntN,Cmd);
  FillAsnPar(Par1,Par2,Chain,Cn1,Cn2,PatP,PatCntP,Cmd);

  Bridge(AntiPar1,AntiPar2,Chain,Cn1,Cn2,PatN,PatCntN);
  Bridge(Par1,Par2,Chain,Cn1,Cn2,PatP,PatCntP);

  for( i=0; i<Chain[Cn1]->NRes; i++ )
    if( AntiPar1[i] == 'N' || Par1[i] == 'P' ) 
      Chain[Cn1]->Rsd[i]->Prop->Asn = 'E';
    else
    if( AntiPar1[i] == 'B' || Par1[i] == 'B' )
      Chain[Cn1]->Rsd[i]->Prop->Asn = 'B';
    else
    if( AntiPar1[i] == 'b' || Par1[i] == 'b' )
      Chain[Cn1]->Rsd[i]->Prop->Asn = 'b';

  for( i=0; i<Chain[Cn2]->NRes; i++ )
    if( Chain[Cn2]->Rsd[i]->Prop->Asn == 'E' )
      continue;
    else
    if( AntiPar2[i] == 'N' || Par2[i] == 'P' ) 
      Chain[Cn2]->Rsd[i]->Prop->Asn = 'E';
    else
    if( AntiPar2[i] == 'B' || Par2[i] == 'B' ) 
      Chain[Cn2]->Rsd[i]->Prop->Asn = 'B';
    else
    if( AntiPar2[i] == 'b' || Par2[i] == 'b' ) 
      Chain[Cn2]->Rsd[i]->Prop->Asn = 'b';

/*
  for( i=0; i<PatCntN; i++ )
    free(PatN[i]);
  for( i=0; i<PatCntP; i++ )
    free(PatP[i]);
*/
  free(PatN);
  free(PatP);
  free(AntiPar1);
  free(Par1);
  free(AntiPar2);
  free(Par2);

}

int Link(HBOND **HBond, CHAIN **Chain, int Cn1, int Cn2, RESIDUE *Res1_1, 
	 RESIDUE *Res1_2, RESIDUE *Res2_2, RESIDUE *Res2_1, RESIDUE *CRes1, 
	 RESIDUE *CRes2, float **PhiPsiMap, PATTERN **Pattern, int *NumPat, 
	 char *Text, float Treshold, COMMAND *Cmd, int Test)
{

  int BondNumber1, BondNumber2, Flag = 0;
  static char *Result[2] = {" NO \n"," YES \n"};
  float Prob1, Prob2, Conf, Coeff;

  if( (BondNumber1 = FindPolInt(HBond,Res1_1,Res1_2)) == ERR ) 
    return(FAILURE);

  if( (BondNumber2 = FindPolInt(HBond,Res2_2,Res2_1)) == ERR ) 
    return(FAILURE);

  if( CRes1 == NULL ) {
    if( CRes2->Prop->PhiZn == ERR || CRes2->Prop->PsiZn == ERR )
      return(FAILURE);
    Conf = PhiPsiMap[CRes2->Prop->PhiZn][CRes2->Prop->PsiZn];
  }
  else
  if( CRes2 == NULL ) {
    if( CRes1->Prop->PhiZn == ERR || CRes1->Prop->PsiZn == ERR )
      return(FAILURE);
    Conf = PhiPsiMap[CRes1->Prop->PhiZn][CRes1->Prop->PsiZn];
  }
  else {
    if( CRes2->Prop->PhiZn == ERR || CRes2->Prop->PsiZn == ERR ||
        CRes1->Prop->PhiZn == ERR || CRes1->Prop->PsiZn == ERR )
      return(FAILURE);
    Conf = 
      0.5*(PhiPsiMap[CRes1->Prop->PhiZn][CRes1->Prop->PsiZn]+
	   PhiPsiMap[CRes2->Prop->PhiZn][CRes2->Prop->PsiZn]);
  }
  Coeff = 1+Cmd->C1_E+Cmd->C2_E*Conf;
  Prob1 = HBond[BondNumber1]->Energy*Coeff;
  Prob2 = HBond[BondNumber2]->Energy*Coeff;

  if( Prob1 < Treshold && Prob2 < Treshold ) {

    if( !Test ) {
      Pattern[*NumPat] = (PATTERN *)ckalloc(sizeof(PATTERN));
      Pattern[*NumPat]->ExistPattern = YES;
      Pattern[*NumPat]->Hb1 = HBond[BondNumber1];
      Pattern[*NumPat]->Hb2 = HBond[BondNumber2];
      Pattern[*NumPat]->Nei1 = NULL;
      Pattern[*NumPat]->Nei2 = NULL;
      strcpy(Pattern[*NumPat]->Type,Text);
      (*NumPat)++;
    }
    Flag = 1;
  }

  if( Cmd->Info && Flag ) {
    fprintf(stdout,"%s %c: %3s %c: %3s | %c: %3s %c: %3s | ",
	    Text,
	    Chain[Cn1]->Id,Res1_1->PDB_ResNumb,
	    Chain[Cn2]->Id,Res1_2->PDB_ResNumb,
	    Chain[Cn2]->Id,Res2_2->PDB_ResNumb,
	    Chain[Cn1]->Id,Res2_1->PDB_ResNumb);
    fprintf(stdout,"%8.6f %6.4f | ", Prob1,HBond[BondNumber1]->Energy);
    fprintf(stdout,"%8.6f %6.4f | ", Prob2,HBond[BondNumber2]->Energy);

    if( CRes1 != NULL && 
        CRes1->Prop->PhiZn != ERR && CRes1->Prop->PsiZn != ERR )
      fprintf(stdout,"%6.4f %2d %2d | ",   
	      PhiPsiMap[CRes1->Prop->PhiZn][CRes1->Prop->PsiZn],
	      CRes1->Prop->PhiZn,CRes1->Prop->PsiZn);
    else
      fprintf(stdout,"000000 00 00 | ");

    if( CRes2 != NULL &&
        CRes2->Prop->PhiZn != ERR && CRes2->Prop->PsiZn != ERR )
      fprintf(stdout,"%6.4f %2d %2d | ",   
	      PhiPsiMap[CRes2->Prop->PhiZn][CRes2->Prop->PsiZn],
	      CRes2->Prop->PhiZn,CRes2->Prop->PsiZn);
    else
      fprintf(stdout,"000000 00 00 | ");
    
    fprintf(stdout,"%s",Result[Flag]);
  }
  
  return(Flag);
  
}

void PrintPatterns(PATTERN **Pat, int NPat, CHAIN **Chain, int Cn1, int Cn2)
{

  register int i;
  int D1, A1, D2, A2;


  for( i=0; i<NPat; i++ ) {
    if( !Pat[i]->ExistPattern ) continue;

    D1 = Pat[i]->Hb1->Dnr->D_Res;
    A1 = Pat[i]->Hb1->Acc->A_Res;
    D2 = Pat[i]->Hb2->Dnr->D_Res;
    A2 = Pat[i]->Hb2->Acc->A_Res;
    
    fprintf(stdout,"%3d %c %c ",
	    i,Pat[i]->Hb1->Dnr->Chain->Id,Pat[i]->Hb2->Dnr->Chain->Id);
    if( Pat[i]->Hb1->Dnr->Chain->Id == Chain[Cn1]->Id )
      fprintf(stdout,"%3s(%3d) %3s(%3d) %3s(%3d) %3s(%3d)",
	      Chain[Cn1]->Rsd[D1]->PDB_ResNumb,D1,
	      Chain[Cn2]->Rsd[A1]->PDB_ResNumb,A1,
	      Chain[Cn2]->Rsd[D2]->PDB_ResNumb,D2,
	      Chain[Cn1]->Rsd[A2]->PDB_ResNumb,A2);
    else
      fprintf(stdout,"%3s(%3d) %3s(%3d) %3s(%3d) %3s(%3d)",
	      Chain[Cn2]->Rsd[D1]->PDB_ResNumb,D1,
	      Chain[Cn1]->Rsd[A1]->PDB_ResNumb,A1,
	      Chain[Cn1]->Rsd[D2]->PDB_ResNumb,D2,
	      Chain[Cn2]->Rsd[A2]->PDB_ResNumb,A2);
    
    if( Pat[i]->Nei1 != NULL ) {
      D1 = Pat[i]->Nei1->Hb1->Dnr->D_Res;
      A1 = Pat[i]->Nei1->Hb1->Acc->A_Res;
      D2 = Pat[i]->Nei1->Hb2->Dnr->D_Res;
      A2 = Pat[i]->Nei1->Hb2->Acc->A_Res;
      
      fprintf(stdout," N1 %c %c ",
	      Pat[i]->Nei1->Hb1->Dnr->Chain->Id,Pat[i]->Nei1->Hb2->Dnr->Chain->Id);
      if( Pat[i]->Nei1->Hb1->Dnr->Chain->Id == Chain[Cn1]->Id )
	fprintf(stdout,"%3s(%3d) %3s(%3d) %3s(%3d) %3s(%3d) ",
		Chain[Cn1]->Rsd[D1]->PDB_ResNumb,D1,
		Chain[Cn2]->Rsd[A1]->PDB_ResNumb,A1,
		Chain[Cn2]->Rsd[D2]->PDB_ResNumb,D2,
		Chain[Cn1]->Rsd[A2]->PDB_ResNumb,A2);
      else
	fprintf(stdout,"%3s(%3d) %3s(%3d) %3s(%3d) %3s(%3d) ",
		Chain[Cn2]->Rsd[D1]->PDB_ResNumb,D1,
		Chain[Cn1]->Rsd[A1]->PDB_ResNumb,A1,
		Chain[Cn1]->Rsd[D2]->PDB_ResNumb,D2,
		Chain[Cn2]->Rsd[A2]->PDB_ResNumb,A2);
    }
    
    if( Pat[i]->Nei2 != NULL ) {
      D1 = Pat[i]->Nei2->Hb1->Dnr->D_Res;
      A1 = Pat[i]->Nei2->Hb1->Acc->A_Res;
      D2 = Pat[i]->Nei2->Hb2->Dnr->D_Res;
      A2 = Pat[i]->Nei2->Hb2->Acc->A_Res;
      fprintf(stdout," N2 %c %c ",
	      Pat[i]->Nei2->Hb1->Dnr->Chain->Id,Pat[i]->Nei2->Hb2->Dnr->Chain->Id);
      if( Pat[i]->Nei2->Hb1->Dnr->Chain->Id == Chain[Cn1]->Id )
	fprintf(stdout,"%3s(%3d) %3s(%3d) %3s(%3d) %3s(%3d) ",
		Chain[Cn1]->Rsd[D1]->PDB_ResNumb,D1,
		Chain[Cn2]->Rsd[A1]->PDB_ResNumb,A1,
		Chain[Cn2]->Rsd[D2]->PDB_ResNumb,D2,
		Chain[Cn1]->Rsd[A2]->PDB_ResNumb,A2);
      else
	fprintf(stdout,"%3s(%3d) %3s(%3d) %3s(%3d) %3s(%3d) ",
		Chain[Cn2]->Rsd[D1]->PDB_ResNumb,D1,
		Chain[Cn1]->Rsd[A1]->PDB_ResNumb,A1,
		Chain[Cn1]->Rsd[D2]->PDB_ResNumb,D2,
		Chain[Cn2]->Rsd[A2]->PDB_ResNumb,A2);
    }
    fprintf(stdout,"\n");
  }
  
}  
    

void Bridge(char *Asn1, char *Asn2, CHAIN **Chain, int Cn1, int Cn2, PATTERN **Pat, int NPat)
{

  register int i;
  int B_Res;
  
  for( i=0; i<NPat; i++ ) {
    if( Pat[i]->Nei1 != NULL || Pat[i]->Nei2 != NULL ) continue;

    if( !strcmp(Pat[i]->Type,"1331") && 
       ( Cn1 != Cn2 || abs(Pat[i]->Hb1->Dnr->D_Res-Pat[i]->Hb1->Acc->A_Res) >= 3 ) ) {
      
      if( Pat[i]->Hb1->Dnr->Chain->Id == Chain[Cn1]->Id ) {
	if( Asn1[Pat[i]->Hb1->Dnr->D_Res] == 'C' )
	  Asn1[Pat[i]->Hb1->Dnr->D_Res] = 'B';
	if( Asn2[Pat[i]->Hb1->Acc->A_Res] == 'C' )
	  Asn2[Pat[i]->Hb1->Acc->A_Res] = 'B';
      }
      else {
	if( Asn2[Pat[i]->Hb1->Dnr->D_Res] == 'C' )
	  Asn2[Pat[i]->Hb1->Dnr->D_Res] = 'B';
	if( Asn1[Pat[i]->Hb1->Acc->A_Res] == 'C' )
	  Asn1[Pat[i]->Hb1->Acc->A_Res] = 'B';
      }

    }
    else
      if( !strcmp(Pat[i]->Type,"3124") && 
	 ( Cn1 != Cn2 || 
	  (abs(Pat[i]->Hb1->Dnr->D_Res-Pat[i]->Hb1->Acc->A_Res) >= 2 &&
	   abs(Pat[i]->Hb2->Dnr->D_Res-Pat[i]->Hb2->Acc->A_Res) >= 2 ) ) ) {
      
	if( Pat[i]->Hb1->Dnr->Chain->Id == Chain[Cn1]->Id ) {
	
	  if( Pat[i]->Hb1->Dnr->D_Res > Pat[i]->Hb2->Acc->A_Res )
	    B_Res = Pat[i]->Hb1->Dnr->D_Res-1;
	  else
	    B_Res = Pat[i]->Hb1->Dnr->D_Res+1;

	  if( Asn1[B_Res] == 'C' )
	    Asn1[B_Res] = 'B';

	  if( Pat[i]->Hb2->Dnr->D_Res > Pat[i]->Hb1->Acc->A_Res )
	    B_Res = Pat[i]->Hb2->Dnr->D_Res-1;
	  else
	    B_Res = Pat[i]->Hb2->Dnr->D_Res+1;

	  if( Asn2[B_Res] == 'C' )
	    Asn2[B_Res] = 'B';
	}
	else {
	  if( Pat[i]->Hb1->Dnr->D_Res > Pat[i]->Hb2->Acc->A_Res )
	    B_Res = Pat[i]->Hb1->Dnr->D_Res-1;
	  else
	    B_Res = Pat[i]->Hb1->Dnr->D_Res+1;

	  if( Asn2[B_Res] == 'C' )
	    Asn2[B_Res] = 'B';

	  if( Pat[i]->Hb2->Dnr->D_Res > Pat[i]->Hb1->Acc->A_Res )
	    B_Res = Pat[i]->Hb2->Dnr->D_Res-1;
	  else
	    B_Res = Pat[i]->Hb2->Dnr->D_Res+1;

	  if( Asn1[B_Res] == 'C' )
	    Asn1[B_Res] = 'B';
	}
      }
      else
	if( ( ( !strcmp(Pat[i]->Type,"3123") || !strcmp(Pat[i]->Type,"1341") ) &&
	     ( Cn1 != Cn2 ||
	      (abs(Pat[i]->Hb1->Dnr->D_Res-Pat[i]->Hb1->Acc->A_Res) > 3 &&
	       abs(Pat[i]->Hb2->Dnr->D_Res-Pat[i]->Hb2->Acc->A_Res) > 3 ) ) ) ) { 

	  if( Pat[i]->Hb1->Dnr->Chain->Id == Chain[Cn1]->Id ) {
	  
	    if( Pat[i]->Hb1->Dnr->D_Res == Pat[i]->Hb2->Acc->A_Res ) {
	    
	      if( Asn1[Pat[i]->Hb1->Dnr->D_Res] == 'C' )
		Asn1[Pat[i]->Hb1->Dnr->D_Res] = 'B';
	    
	      if( Pat[i]->Hb2->Dnr->D_Res > Pat[i]->Hb1->Acc->A_Res )
		B_Res = Pat[i]->Hb2->Dnr->D_Res-1;
	      else
		B_Res = Pat[i]->Hb2->Dnr->D_Res+1;
	    
	      if( Asn2[B_Res] == 'C' )
		Asn2[B_Res] = 'B';
	    }
	    else {
	      if( Pat[i]->Hb2->Dnr->D_Res == Pat[i]->Hb1->Acc->A_Res )
	      
		if( Asn2[Pat[i]->Hb2->Dnr->D_Res] == 'C' )
		  Asn2[Pat[i]->Hb2->Dnr->D_Res] = 'B';
	    
	      if( Pat[i]->Hb1->Dnr->D_Res > Pat[i]->Hb2->Acc->A_Res )
		B_Res = Pat[i]->Hb1->Dnr->D_Res-1;
	      else
		B_Res = Pat[i]->Hb1->Dnr->D_Res+1;
	    
	      if( Asn1[B_Res] == 'C' )
		Asn1[B_Res] = 'B';
	    }
	  }
	}
	else
	  if( ( !strcmp(Pat[i]->Type,"13B1") || !strcmp(Pat[i]->Type,"133A") ) &&
	     ( Cn1 != Cn2 ||
	      (abs(Pat[i]->Hb1->Dnr->D_Res-Pat[i]->Hb1->Acc->A_Res) > 4 &&
	       abs(Pat[i]->Hb2->Dnr->D_Res-Pat[i]->Hb2->Acc->A_Res) > 4 ) ) ) { 

	    if( Pat[i]->Hb1->Dnr->Chain->Id == Chain[Cn1]->Id ) {
	  
	      if( Pat[i]->Hb1->Dnr->D_Res == Pat[i]->Hb2->Acc->A_Res ) {
	    
		if( Asn1[Pat[i]->Hb1->Dnr->D_Res] == 'C' )
		  Asn1[Pat[i]->Hb1->Dnr->D_Res] = 'B';
	    
		if( Pat[i]->Hb2->Dnr->D_Res > Pat[i]->Hb1->Acc->A_Res )
		  B_Res = Pat[i]->Hb2->Dnr->D_Res-1;
		else
		  B_Res = Pat[i]->Hb2->Dnr->D_Res+1;
	    
		if( Asn2[B_Res] == 'C' )
		  Asn2[B_Res] = 'B';
	      }
	      else {
		if( Pat[i]->Hb2->Dnr->D_Res == Pat[i]->Hb1->Acc->A_Res )
	      
		  if( Asn2[Pat[i]->Hb2->Dnr->D_Res] == 'C' )
		    Asn2[Pat[i]->Hb2->Dnr->D_Res] = 'b';
	    
		if( Pat[i]->Hb1->Dnr->D_Res > Pat[i]->Hb2->Acc->A_Res )
		  B_Res = Pat[i]->Hb1->Dnr->D_Res-1;
		else
		  B_Res = Pat[i]->Hb1->Dnr->D_Res+1;
	    
		if( Asn1[B_Res] == 'C' )
		  Asn1[B_Res] = 'b';
	      }
	    }
	  }
    
  }
}
















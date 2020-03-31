#include "stride.h"

void FillAsnAntiPar(char *Asn1, char *Asn2, CHAIN **Chain, int Cn1, int Cn2, 
		    PATTERN **Pat, int NPat, COMMAND *Cmd)
{
  register int i, j;
  int Beg1, Beg2, End1, End2;
  int B1D, B1A, B2D, B2A, E1D, E1A, E2D, E2A;
  char B1DCn, B1ACn, B2DCn, B2ACn, E1DCn, E1ACn, E2DCn, E2ACn, Beg1Cn, Beg2Cn; 
  PATTERN *CurrPat, *PrevPat;;

  for( i=0; i<NPat; i++ ) {

    if( Pat[i]->Nei1 != NULL && Pat[i]->Nei2 == NULL )
      CurrPat = Pat[i]->Nei1;
    else
    if( Pat[i]->Nei2 != NULL && Pat[i]->Nei1 == NULL ) 
      CurrPat = Pat[i]->Nei2;
    else 
      continue;
    
    if( Cmd->Info ) {
      fprintf(stdout,"From: %c %c ",
	      Pat[i]->Hb1->Dnr->Chain->Id,Pat[i]->Hb2->Dnr->Chain->Id);
      if( Pat[i]->Hb1->Dnr->Chain->Id == Chain[Cn1]->Id )
	fprintf(stdout,"%s %s %s %s \n",
	    Chain[Cn1]->Rsd[Pat[i]->Hb1->Dnr->D_Res]->PDB_ResNumb,
	    Chain[Cn2]->Rsd[Pat[i]->Hb1->Acc->A_Res]->PDB_ResNumb,
	    Chain[Cn2]->Rsd[Pat[i]->Hb2->Dnr->D_Res]->PDB_ResNumb,
	    Chain[Cn1]->Rsd[Pat[i]->Hb2->Acc->A_Res]->PDB_ResNumb);
      else
	fprintf(stdout,"%s %s %s %s \n",
	    Chain[Cn2]->Rsd[Pat[i]->Hb1->Dnr->D_Res]->PDB_ResNumb,
	    Chain[Cn1]->Rsd[Pat[i]->Hb1->Acc->A_Res]->PDB_ResNumb,
	    Chain[Cn1]->Rsd[Pat[i]->Hb2->Dnr->D_Res]->PDB_ResNumb,
	    Chain[Cn2]->Rsd[Pat[i]->Hb2->Acc->A_Res]->PDB_ResNumb);
    }

    PrevPat = Pat[i];
    while( CurrPat->Nei1 != NULL && CurrPat->Nei2 != NULL ) {
      
      if( (CurrPat->Nei1->Nei1 == CurrPat || CurrPat->Nei1->Nei2 == CurrPat) && 
	 CurrPat->Nei1 != PrevPat ) {
	PrevPat = CurrPat;
	CurrPat = CurrPat->Nei1;
      }
      else 
      if( (CurrPat->Nei2->Nei1 == CurrPat || CurrPat->Nei2->Nei2 == CurrPat) && 
	 CurrPat->Nei2 != PrevPat ) {
	PrevPat = CurrPat;
	CurrPat = CurrPat->Nei2;
      }
      else {
	fprintf(stdout,"Cycle Anti%s%c i = %d \n",Chain[Cn1]->File,Chain[Cn1]->Id,i);
	break;
      }
    }  
    
    if( Cmd->Info ) {
      fprintf(stdout,"To: %c %c ",
	      CurrPat->Hb1->Dnr->Chain->Id,CurrPat->Hb2->Dnr->Chain->Id);
      if( CurrPat->Hb1->Dnr->Chain->Id == Chain[Cn1]->Id )
	fprintf(stdout,"%s %s %s %s \n",
	    Chain[Cn1]->Rsd[CurrPat->Hb1->Dnr->D_Res]->PDB_ResNumb,
	    Chain[Cn2]->Rsd[CurrPat->Hb1->Acc->A_Res]->PDB_ResNumb,
	    Chain[Cn2]->Rsd[CurrPat->Hb2->Dnr->D_Res]->PDB_ResNumb,
	    Chain[Cn1]->Rsd[CurrPat->Hb2->Acc->A_Res]->PDB_ResNumb);
      else
	fprintf(stdout,"%s %s %s %s \n",
	    Chain[Cn2]->Rsd[CurrPat->Hb1->Dnr->D_Res]->PDB_ResNumb,
	    Chain[Cn1]->Rsd[CurrPat->Hb1->Acc->A_Res]->PDB_ResNumb,
	    Chain[Cn1]->Rsd[CurrPat->Hb2->Dnr->D_Res]->PDB_ResNumb,
	    Chain[Cn2]->Rsd[CurrPat->Hb2->Acc->A_Res]->PDB_ResNumb);
    }

    Alias(&B1D,&B1A,&B2D,&B2A,&B1DCn,&B1ACn,&B2DCn,&B2ACn,Pat[i]);
    Alias(&E1D,&E1A,&E2D,&E2A,&E1DCn,&E1ACn,&E2DCn,&E2ACn,CurrPat);

    if( (Cn1 != Cn2 || E1D - B2A <  E2D - B2A ) &&
        ( MakeEnds(&Beg1,B1D,B2A,&Beg1Cn,B1DCn,&End1,E2A,E1D,E2ACn,&Beg2,E2D,E1A,&Beg2Cn,E2DCn,
		   &End2,B1A,B2D,B1ACn,Pat,NPat) ||
          MakeEnds(&Beg1,B1D,B2A,&Beg1Cn,B1DCn,&End1,E1D,E2A,E1DCn,&Beg2,E1A,E2D,&Beg2Cn,E1ACn,
		   &End2,B1A,B2D,B1ACn,Pat,NPat) ) )
      ;
    else
    if( ( Cn1 != Cn2 || E2D - B2A <  E1D - B2A ) && 
        ( MakeEnds(&Beg1,B1D,B2A,&Beg1Cn,B1DCn,&End1,E1A,E2D,E1ACn,&Beg2,E1D,E2A,&Beg2Cn,E1DCn,
		   &End2,B1A,B2D,B1ACn,Pat,NPat) ||
          MakeEnds(&Beg1,B1D,B2A,&Beg1Cn,B1DCn,&End1,E2D,E1A,E2DCn,&Beg2,E2A,E1D,&Beg2Cn,E2ACn,
		   &End2,B1A,B2D,B1ACn,Pat,NPat) ) )
      ;
    else
    if( ( Cn1 != Cn2 || B2A - E1D < B2A - E2D ) && 
        ( MakeEnds(&Beg1,B1A,B2D,&Beg1Cn,B1ACn,&End1,E2D,E1A,E2DCn,&Beg2,E2A,E1D,&Beg2Cn,E2ACn,
		   &End2,B1D,B2A,B1DCn,Pat,NPat) ||
          MakeEnds(&Beg1,B1A,B2D,&Beg1Cn,B1ACn,&End1,E1A,E2D,E1ACn,&Beg2,E1D,E2A,&Beg2Cn,E1DCn,
		   &End2,B1D,B2A,B1DCn,Pat,NPat) ) )
      ;
    else
    if( ( Cn1 != Cn2 || B2A - E2D < B2A - E1D ) && 
        ( MakeEnds(&Beg1,B1A,B2D,&Beg1Cn,B1ACn,&End1,E1D,E2A,E1DCn,&Beg2,E1A,E2D,&Beg2Cn,E1ACn,
		   &End2,B1D,B2A,B1DCn,Pat,NPat) ||
          MakeEnds(&Beg1,B1A,B2D,&Beg1Cn,B1ACn,&End1,E2A,E1D,E2ACn,&Beg2,E2D,E1A,&Beg2Cn,E2DCn,
		   &End2,B1D,B2A,B1DCn,Pat,NPat) ) )
      ;
    else
    if( ( Cn1 != Cn2 || B1D - E2A <  B2D - E2A ) && 
        ( MakeEnds(&Beg1,E1D,E2A,&Beg1Cn,E1DCn,&End1,B2A,B1D,B2ACn,&Beg2,B2D,B1A,&Beg2Cn,B2DCn,
		   &End2,E1A,E2D,E1ACn,Pat,NPat) ||
          MakeEnds(&Beg1,E1D,E2A,&Beg1Cn,E1DCn,&End1,B1D,B2A,B1DCn,&Beg2,B1A,B2D,&Beg2Cn,B1ACn,
		   &End2,E1A,E2D,E1ACn,Pat,NPat) ) )
      ;
    else
    if( ( Cn1 != Cn2 || B2D - E2A <  B1D - E2A ) && 
        ( MakeEnds(&Beg1,E1D,E2A,&Beg1Cn,E1DCn,&End1,B1A,B2D,B1ACn,&Beg2,B1D,B2A,&Beg2Cn,B1DCn,
		   &End2,E1A,E2D,E1ACn,Pat,NPat) ||
          MakeEnds(&Beg1,E1D,E2A,&Beg1Cn,E1DCn,&End1,B2D,B1A,B2DCn,&Beg2,B2A,B1D,&Beg2Cn,B2ACn,
		   &End2,E1A,E2D,E1ACn,Pat,NPat) ) )
      ;
    else
    if( ( Cn1 != Cn2 || E2A - B1D < E2A - B2D ) && 
        ( MakeEnds(&Beg1,E1A,E2D,&Beg1Cn,E1ACn,&End1,B2D,B1A,B2DCn,&Beg2,B2A,B1D,&Beg2Cn,B2ACn,
		   &End2,E1D,E2A,E1DCn,Pat,NPat) ||
          MakeEnds(&Beg1,E1A,E2D,&Beg1Cn,E1ACn,&End1,B1A,B2D,B1ACn,&Beg2,B1D,B2A,&Beg2Cn,B1DCn,
		   &End2,E1D,E2A,E1DCn,Pat,NPat) ) )
      ;
    else
    if( ( Cn1 != Cn2 || E2A - B2D < E2A - B1D ) && 
        ( MakeEnds(&Beg1,E1A,E2D,&Beg1Cn,E1ACn,&End1,B1D,B2A,B1DCn,&Beg2,B1A,B2D,&Beg2Cn,B1ACn,
		   &End2,E1D,E2A,E1DCn,Pat,NPat) ||
          MakeEnds(&Beg1,E1A,E2D,&Beg1Cn,E1ACn,&End1,B2A,B1D,B2ACn,&Beg2,B2D,B1A,&Beg2Cn,B2DCn,
		   &End2,E1D,E2A,E1DCn,Pat,NPat) ) )
      ;
    else {
/*      fprintf(stdout,"Ne tot variant.. Anti.. %s%c\n",Chain[Cn1]->File,Chain[Cn1]->Id);*/
      continue;
    }


    if( Beg1Cn == Chain[Cn1]->Id ) {
      for( j=Beg1; j<=End1; j++ ) 
	Asn1[j] = 'N';
      for( j=Beg2; j<=End2; j++ ) 
	Asn2[j] = 'N';
    }
    else {
      for( j=Beg1; j<=End1; j++ ) 
	Asn2[j] = 'N';
      for( j=Beg2; j<=End2; j++ ) 
	Asn1[j] = 'N';
    }

    Pat[i]->Nei1 = NULL;
    Pat[i]->Nei2 = NULL;
    CurrPat->Nei1 = NULL;
    CurrPat->Nei2 = NULL;

  }
}
  

void FillAsnPar(char *Asn1, char *Asn2, CHAIN **Chain, int Cn1, int Cn2, 
		PATTERN **Pat, int NPat, COMMAND *Cmd)
{
  register int i, j;
  int Beg1, Beg2, End1, End2;
  int B1D, B1A, B2D, B2A, E1D, E1A, E2D, E2A; 
  char B1DCn, B1ACn, B2DCn, B2ACn, E1DCn, E1ACn, E2DCn, E2ACn, Beg1Cn, Beg2Cn; 
  PATTERN *CurrPat, *PrevPat;;

  for( i=0; i<NPat; i++ ) {
    
    if( Pat[i]->Nei1 != NULL && Pat[i]->Nei2 == NULL )
      CurrPat = Pat[i]->Nei1;
    else
    if( Pat[i]->Nei2 != NULL && Pat[i]->Nei1 == NULL ) 
      CurrPat = Pat[i]->Nei2;
    else 
      continue;
    
    if( Cmd->Info ) {
      fprintf(stdout,"From: %c %c ",
	      Pat[i]->Hb1->Dnr->Chain->Id,Pat[i]->Hb2->Dnr->Chain->Id);
      if( Pat[i]->Hb1->Dnr->Chain->Id == Chain[Cn1]->Id )
	fprintf(stdout,"%s %s %s %s \n",
	    Chain[Cn1]->Rsd[Pat[i]->Hb1->Dnr->D_Res]->PDB_ResNumb,
	    Chain[Cn2]->Rsd[Pat[i]->Hb1->Acc->A_Res]->PDB_ResNumb,
	    Chain[Cn2]->Rsd[Pat[i]->Hb2->Dnr->D_Res]->PDB_ResNumb,
	    Chain[Cn1]->Rsd[Pat[i]->Hb2->Acc->A_Res]->PDB_ResNumb);
      else
	fprintf(stdout,"%s %s %s %s \n",
	    Chain[Cn2]->Rsd[Pat[i]->Hb1->Dnr->D_Res]->PDB_ResNumb,
	    Chain[Cn1]->Rsd[Pat[i]->Hb1->Acc->A_Res]->PDB_ResNumb,
	    Chain[Cn1]->Rsd[Pat[i]->Hb2->Dnr->D_Res]->PDB_ResNumb,
	    Chain[Cn2]->Rsd[Pat[i]->Hb2->Acc->A_Res]->PDB_ResNumb);
    }

    PrevPat = Pat[i];
    while( CurrPat->Nei1 != NULL && CurrPat->Nei2 != NULL ) {
      
      if( (CurrPat->Nei1->Nei1 == CurrPat || CurrPat->Nei1->Nei2 == CurrPat) && 
	 CurrPat->Nei1 != PrevPat ) {
	PrevPat = CurrPat;
	CurrPat = CurrPat->Nei1;
      }
      else {
	PrevPat = CurrPat;
	CurrPat = CurrPat->Nei2;
      }
    }  

    if( Cmd->Info ) {
      fprintf(stdout,"To: %c %c ",
	      CurrPat->Hb1->Dnr->Chain->Id,CurrPat->Hb2->Dnr->Chain->Id);
      if( CurrPat->Hb1->Dnr->Chain->Id == Chain[Cn1]->Id )
	fprintf(stdout,"%s %s %s %s \n",
	    Chain[Cn1]->Rsd[CurrPat->Hb1->Dnr->D_Res]->PDB_ResNumb,
	    Chain[Cn2]->Rsd[CurrPat->Hb1->Acc->A_Res]->PDB_ResNumb,
	    Chain[Cn2]->Rsd[CurrPat->Hb2->Dnr->D_Res]->PDB_ResNumb,
	    Chain[Cn1]->Rsd[CurrPat->Hb2->Acc->A_Res]->PDB_ResNumb);
      else
	fprintf(stdout,"%s %s %s %s \n",
	    Chain[Cn2]->Rsd[CurrPat->Hb1->Dnr->D_Res]->PDB_ResNumb,
	    Chain[Cn1]->Rsd[CurrPat->Hb1->Acc->A_Res]->PDB_ResNumb,
	    Chain[Cn1]->Rsd[CurrPat->Hb2->Dnr->D_Res]->PDB_ResNumb,
	    Chain[Cn2]->Rsd[CurrPat->Hb2->Acc->A_Res]->PDB_ResNumb);
    }

    Alias(&B1D,&B1A,&B2D,&B2A,&B1DCn,&B1ACn,&B2DCn,&B2ACn,Pat[i]);
    Alias(&E1D,&E1A,&E2D,&E2A,&E1DCn,&E1ACn,&E2DCn,&E2ACn,CurrPat);

    if( ( Cn1 != Cn2 || abs(E1D-B2A) < abs(E2D-B2A) ) && 
        ( MakeEnds(&Beg1,B1D,B2A,&Beg1Cn,B1DCn,&End1,E2A,E1D,E2ACn,&Beg2,B1A,B2D,&Beg2Cn,B1ACn,
		   &End2,E2D,E1A,E2DCn,Pat,NPat) ||
          MakeEnds(&Beg1,B1D,B2A,&Beg1Cn,B1DCn,&End1,E1D,E2A,E1DCn,&Beg2,B1A,B2D,&Beg2Cn,B1ACn,
		   &End2,E1A,E2D,E1ACn,Pat,NPat) ) )
      ;
    else
    if( ( Cn1 != Cn2 || abs(E2D-B2A) < abs(E1D-B2A) ) && 
        ( MakeEnds(&Beg1,B1D,B2A,&Beg1Cn,B1DCn,&End1,E1A,E2D,E1ACn,&Beg2,B1A,B2D,&Beg2Cn,B1ACn,
		   &End2,E1D,E2A,E1DCn,Pat,NPat) ||
          MakeEnds(&Beg1,B1D,B2A,&Beg1Cn,B1DCn,&End1,E2D,E1A,E2DCn,&Beg2,B1A,B2D,&Beg2Cn,B1ACn,
		   &End2,E2A,E1D,E2ACn,Pat,NPat) ) )
      ;
    else
    if( ( Cn1 != Cn2 || abs(B2A-E1D) < abs(B2A-E2D) ) && 
        ( MakeEnds(&Beg1,B1A,B2D,&Beg1Cn,B1ACn,&End1,E2D,E1A,E2DCn,&Beg2,B1D,B2A,&Beg2Cn,B1DCn,
		   &End2,E2A,E1D,E2ACn,Pat,NPat) ||
          MakeEnds(&Beg1,B1A,B2D,&Beg1Cn,B1ACn,&End1,E1A,E2D,E1ACn,&Beg2,B1D,B2A,&Beg2Cn,B1DCn,
		   &End2,E1D,E2A,E1DCn,Pat,NPat) ) )
      ;
    else
    if( ( Cn1 != Cn2 || abs(B2A-E2D) < abs(B2A-E1D) ) && 
        ( MakeEnds(&Beg1,B1A,B2D,&Beg1Cn,B1ACn,&End1,E1D,E2A,E1DCn,&Beg2,B1D,B2A,&Beg2Cn,B1DCn,
		   &End2,E1A,E2D,E1ACn,Pat,NPat) ||
          MakeEnds(&Beg1,B1A,B2D,&Beg1Cn,B1ACn,&End1,E2A,E1D,E2ACn,&Beg2,B1D,B2A,&Beg2Cn,B1DCn,
		   &End2,E2D,E1A,E2DCn,Pat,NPat) ) )
      ;
    else
    if( ( Cn1 != Cn2 || abs(B1D-E2A) < abs(B2D-E2A) ) && 
        ( MakeEnds(&Beg1,E1D,E2A,&Beg1Cn,E1DCn,&End1,B2A,B1D,B2ACn,&Beg2,E1A,E2D,&Beg2Cn,E1ACn,
		   &End2,B2D,B1A,B2DCn,Pat,NPat) ||
          MakeEnds(&Beg1,E1D,E2A,&Beg1Cn,E1DCn,&End1,B1D,B2A,B1DCn,&Beg2,E1A,E2D,&Beg2Cn,E1ACn,
		   &End2,B1A,B2D,B1ACn,Pat,NPat) ) )
      ;
    else
    if( ( Cn1 != Cn2 || abs(B2D-E2A) < abs(B1D-E2A) ) && 
        ( MakeEnds(&Beg1,E1D,E2A,&Beg1Cn,E1DCn,&End1,B1A,B2D,B1ACn,&Beg2,E1A,E2D,&Beg2Cn,E1ACn,
		   &End2,B1D,B2A,B1DCn,Pat,NPat) ||
          MakeEnds(&Beg1,E1D,E2A,&Beg1Cn,E1DCn,&End1,B2D,B1A,B2DCn,&Beg2,E1A,E2D,&Beg2Cn,E1ACn,
		   &End2,B2A,B1D,B2ACn,Pat,NPat) ) )
      ;
    else
    if( ( Cn1 != Cn2 || abs(E2A-B1D) < abs(E2A-B2D) ) && 
        ( MakeEnds(&Beg1,E1A,E2D,&Beg1Cn,E1ACn,&End1,B2D,B1A,B2DCn,&Beg2,E1D,E2A,&Beg2Cn,E1DCn,
		   &End2,B2A,B1D,B2ACn,Pat,NPat) ||
          MakeEnds(&Beg1,E1A,E2D,&Beg1Cn,E1ACn,&End1,B1A,B2D,B1ACn,&Beg2,E1D,E2A,&Beg2Cn,E1DCn,
		   &End2,B1D,B2A,B1DCn,Pat,NPat) ) )
      ;
    else
    if( ( Cn1 != Cn2 || abs(E2A-B2D) < abs(E2A-B1D) ) && 
        ( MakeEnds(&Beg1,E1A,E2D,&Beg1Cn,E1ACn,&End1,B1D,B2A,B1DCn,&Beg2,E1D,E2A,&Beg2Cn,E1DCn,
		   &End2,B1A,B2D,B1ACn,Pat,NPat) ||
          MakeEnds(&Beg1,E1A,E2D,&Beg1Cn,E1ACn,&End1,B2A,B1D,B2ACn,&Beg2,E1D,E2A,&Beg2Cn,E1DCn,
		   &End2,B2D,B1A,B2DCn,Pat,NPat) ) )
      ;
    else {
/*      fprintf(stdout,"Ne tot variant.. Par %s%c\n",Chain[Cn1]->File,Chain[Cn1]->Id);*/
      continue;
    }

    if( Beg1Cn == Chain[Cn1]->Id ) {
      for( j=Beg1; j<=End1; j++ ) Asn1[j] = 'P';
      for( j=Beg2; j<=End2; j++ ) Asn2[j] = 'P';
    }
    else {
      for( j=Beg1; j<=End1; j++ ) Asn2[j] = 'P';
      for( j=Beg2; j<=End2; j++ ) Asn1[j] = 'P';
    }

    Pat[i]->Nei1 = NULL;
    Pat[i]->Nei2 = NULL;
    CurrPat->Nei1 = NULL;
    CurrPat->Nei2 = NULL;

  }
}
  

int MakeEnds(int *Beg1, int ResBeg1, int NeiBeg1, char *Beg1Cn, char ResBeg1Cn, int *End1, 
	     int ResEnd1, int NeiEnd1, char ResEnd1Cn, int *Beg2, int ResBeg2, int NeiBeg2, 
	     char *Beg2Cn, char ResBeg2Cn, int *End2, int ResEnd2, int NeiEnd2, 
	     char ResEnd2Cn, PATTERN **Pat, int NPat)
{

  register int i;
  int Flag1 = 0, Flag2 = 0;


  if( ResBeg1 <= NeiBeg1 && NeiBeg1 <= NeiEnd1 && NeiEnd1 <= ResEnd1 &&
      ResBeg2 <= NeiBeg2 && NeiBeg2 <= NeiEnd2 && NeiEnd2 <= ResEnd2 &&
      ResBeg1Cn == ResEnd1Cn && ResBeg2Cn == ResEnd2Cn ) {

    *Beg1 = ResBeg1;
    *End1 = ResEnd1;
    *Beg2 = ResBeg2;
    *End2 = ResEnd2;
    *Beg1Cn = ResBeg1Cn;
    *Beg2Cn = ResBeg2Cn;
    
    for( i=0; i<NPat && (Flag1 == 0 || Flag2 == 0); i++ ) {
      if( ( (Pat[i]->Hb1->Dnr->D_Res == (*Beg1) 
	     && Pat[i]->Hb1->Acc->A_Res == (*End2)
	     && Pat[i]->Hb1->Dnr->Chain->Id == (*Beg1Cn)
	     && Pat[i]->Hb1->Acc->Chain->Id == (*Beg2Cn) ) 
	   ||
	   (Pat[i]->Hb1->Acc->A_Res == (*Beg1) 
	    && Pat[i]->Hb1->Dnr->D_Res == (*End2) 
	    && Pat[i]->Hb1->Acc->Chain->Id == (*Beg1Cn)
	    && Pat[i]->Hb1->Dnr->Chain->Id == (*Beg2Cn) ) ) 
	 && Pat[i]->Hb1->Dnr->D_Res == Pat[i]->Hb2->Acc->A_Res 
	 && Pat[i]->Hb2->Dnr->D_Res == Pat[i]->Hb1->Acc->A_Res )
	Flag1 = 1; 
      if( ( (Pat[i]->Hb1->Dnr->D_Res == (*Beg2) 
	     && Pat[i]->Hb1->Acc->A_Res == (*End1) 
	     && Pat[i]->Hb1->Dnr->Chain->Id == (*Beg2Cn)
	     && Pat[i]->Hb1->Acc->Chain->Id == (*Beg1Cn) ) 
	   ||
	   (Pat[i]->Hb1->Acc->A_Res == (*Beg2) 
	    && Pat[i]->Hb1->Dnr->D_Res == (*End1) 
	    && Pat[i]->Hb1->Acc->Chain->Id == (*Beg2Cn)
	    && Pat[i]->Hb1->Dnr->Chain->Id == (*Beg1Cn) ) ) 
	 && Pat[i]->Hb1->Dnr->D_Res == Pat[i]->Hb2->Acc->A_Res 
	 && Pat[i]->Hb2->Dnr->D_Res == Pat[i]->Hb1->Acc->A_Res )
	Flag2 = 1; 
    }
    
    if( !Flag1 ) {
      if( *Beg1 != NeiBeg1 ) (*Beg1)++;
      if( *End2 != NeiEnd2 ) (*End2)--;
    }
    
    if( !Flag2 ) {
      if( *End1 != NeiEnd1 ) (*End1)--;
      if( *Beg2 != NeiBeg2 ) (*Beg2)++;
    }
    return(SUCCESS);
  }
  
  return(FAILURE);
}
  
  
void FilterAntiPar(PATTERN **Pat, int NPat)
{
 
  register int i, j;
  int I1A, I1D, I2A, I2D, J1A, J1D, J2A, J2D;
  char I1ACn, I1DCn, I2ACn, I2DCn, J1ACn, J1DCn, J2ACn, J2DCn;

  for( i=0; i<NPat; i++ ) {
    
    if( !Pat[i]->ExistPattern ) continue;

    Alias(&I1D,&I1A,&I2D,&I2A,&I1DCn,&I1ACn,&I2DCn,&I2ACn,Pat[i]);
    
    for( j=0; j<NPat; j++ ) {
      
      if( j == i || !Pat[j]->ExistPattern ) continue;
      
      Alias(&J1D,&J1A,&J2D,&J2A,&J1DCn,&J1ACn,&J2DCn,&J2ACn,Pat[j]);
      
      if( J1D == J2A && J2D == J1A && I1D != I2A && I2D != I1A &&
	 ( (J1D == I1D && J1A == I1A) ||  (J1D == I1A && J1A == I1D) || 
	   (J1D == I2A && J1A == I2D) ||  (J1D == I2D && J1A == I2A) ) ) continue;

      if( ( ( I1D < I2A || I2D < I1A ) && 
	   ( (J1A <= I2A && J1A >= I1D && J2D <= I2A && J2D >= I1D && J2DCn == I1DCn &&
	      J2A <= I1A && J2A >= I2D && J1D <= I1A && J1D >= I2D && J1DCn == I2DCn) ||
	     (J2A <= I2A && J2A >= I1D && J1D <= I2A && J1D >= I1D && J1DCn == I1DCn &&
	      J1A <= I1A && J1A >= I2D && J2D <= I1A && J2D >= I2D && J2DCn == I2DCn) ) ) || 
	  ( ( I1D > I2A || I2D > I1A ) && 
	   ( (J1A >= I2A && J1A <= I1D && J2D >= I2A && J2D <= I1D && J2DCn == I1DCn &&
	      J2A >= I1A && J2A <= I2D && J1D >= I1A && J1D <= I2D && J1DCn == I2DCn) ||
	     (J2A >= I2A && J2A <= I1D && J1D >= I2A && J1D <= I1D && J1DCn == I1DCn &&
	      J1A >= I1A && J1A <= I2D && J2D >= I1A && J2D <= I2D && J2DCn == I2DCn) ) ) ) {
	Pat[j]->ExistPattern = NO;
      }
    }
  }
}

void FilterPar(PATTERN **Pat, int NPat)
{
 
  register int i, j;
  int I1A, I1D, I2A, I2D, J1A, J1D, J2A, J2D;
  char I1ACn, I1DCn, I2ACn, I2DCn, J1ACn, J1DCn, J2ACn, J2DCn;

  for( i=0; i<NPat; i++ ) {
    
    if( !Pat[i]->ExistPattern ) continue;

    Alias(&I1D,&I1A,&I2D,&I2A,&I1DCn,&I1ACn,&I2DCn,&I2ACn,Pat[i]);
    
    for( j=0; j<NPat; j++ ) {
      
      if( j == i || !Pat[j]->ExistPattern ) continue;
      
      Alias(&J1D,&J1A,&J2D,&J2A,&J1DCn,&J1ACn,&J2DCn,&J2ACn,Pat[j]);
      
      if( ( ( I1A >= I2D && I1D >= I2A ) && 
	   ( (J1A >= I2A && J1A <= I1D && J2D >= I2A && J2D <= I1D && J2DCn == I1DCn &&
	      J2A <= I1A && J2A >= I2D && J1D <= I1A && J1D >= I2D && J1DCn == I2DCn) ||
	     (J2A >= I2A && J2A <= I1D && J1D >= I2A && J1D <= I1D && J1DCn == I1DCn &&
	      J1A <= I1A && J1A >= I2D && J2D <= I1A && J2D >= I2D && J2DCn == I2DCn) ) ) || 

	  ( I2A >= I1D && I2D >= I1A  && 
	   ( (J1A <= I2A && J1A >= I1D && J2D <= I2A && J2D >= I1D && J2DCn == I1DCn &&
	      J2A >= I1A && J2A <= I2D && J1D >= I1A && J1D <= I2D && J1DCn == I2DCn) ||

	     (J2A <= I2A && J2A >= I1D && J1D <= I2A && J1D >= I1D && J1DCn == I1DCn &&
	      J1A >= I1A && J1A <= I2D && J2D >= I1A && J2D <= I2D && J2DCn == I2DCn) ) ) ) {
	Pat[j]->ExistPattern = NO;
      }
    }
  }
}
      
void Alias(int *D1,int *A1,int *D2,int *A2,char *D1Cn,char *A1Cn,char *D2Cn,char *A2Cn,
	  PATTERN *Pat)
{
    *D1 = Pat->Hb1->Dnr->D_Res;
    *A1 = Pat->Hb1->Acc->A_Res; 
    *D2 = Pat->Hb2->Dnr->D_Res;
    *A2 = Pat->Hb2->Acc->A_Res;
    *D1Cn = Pat->Hb1->Dnr->Chain->Id;
    *A1Cn = Pat->Hb1->Acc->Chain->Id;
    *D2Cn = Pat->Hb2->Dnr->Chain->Id;
    *A2Cn = Pat->Hb2->Acc->Chain->Id;
}		      


#include "stride.h"

void MergePatternsAntiPar(PATTERN **Pat, int NPat)
{
  register int i, j;
  int DB, DW, MinDB1, MinDB2, MinDW1, MinDW2, Min, Lnk1A, Lnk1D;
  int I1A, I1D, I2A, I2D, J1A, J1D, J2A, J2D;
  char I1ACn, I1DCn, I2ACn, I2DCn, J1ACn, J1DCn, J2ACn, J2DCn;

  for( i=0; i<NPat; i++ ) {

    if( !Pat[i]->ExistPattern ) continue;

    MinDB1 = MinDB2 = MinDW1 = MinDW2 = 1000;
    Min = ERR;
    Lnk1D = Lnk1A = ERR;

    Alias(&I1D,&I1A,&I2D,&I2A,&I1DCn,&I1ACn,&I2DCn,&I2ACn,Pat[i]);

    for( j=0; j<NPat; j++ ) {

      if( i == j || !Pat[j]->ExistPattern ) continue;

      Alias(&J1D,&J1A,&J2D,&J2A,&J1DCn,&J1ACn,&J2DCn,&J2ACn,Pat[j]);
      
      if( Near(I1D,J1D,J1A,I1A,J2A,J2D,I2A,I2D,I1DCn,J1DCn,J1ACn,I1ACn,&DB,&DW) && 
	  ((DB < MinDB1 && DW <= MinDW1) || (DB <= MinDB1 && DW < MinDW1) ) && 
	 RightSide(J1A,J1D,I1A,I1D,I2A,I2D) )
	JoinNeighbours(&Lnk1A,J2D,&Lnk1D,J2A,&Pat[i]->Nei1,Pat[j],&MinDB1,DB,&MinDW1,DW,&Min,j);

      if( Near(I1D,J1A,J1D,I1A,J2D,J2A,I2A,I2D,I1DCn,J1ACn,J1DCn,I1ACn,&DB,&DW) && 
	  ((DB < MinDB1 && DW <= MinDW1) || (DB <= MinDB1 && DW < MinDW1) ) && 
	 RightSide(J1D,J1A,I1A,I1D,I2A,I2D) )
	JoinNeighbours(&Lnk1A,J2A,&Lnk1D,J2D,&Pat[i]->Nei1,Pat[j],&MinDB1,DB,&MinDW1,DW,&Min,j);

      if( Near(I1D,J2D,J2A,I1A,J1A,J1D,I2A,I2D,I1DCn,J2DCn,J2ACn,I1ACn,&DB,&DW) && 
	  ((DB < MinDB1 && DW <= MinDW1) || (DB <= MinDB1 && DW < MinDW1) ) && 
	 RightSide(J2A,J2D,I1A,I1D,I2A,I2D) )
	JoinNeighbours(&Lnk1A,J1D,&Lnk1D,J1A,&Pat[i]->Nei1,Pat[j],&MinDB1,DB,&MinDW1,DW,&Min,j);

      if( Near(I1D,J2A,J2D,I1A,J1D,J1A,I2A,I2D,I1DCn,J2ACn,J2DCn,I1ACn,&DB,&DW) && 
	  ((DB < MinDB1 && DW <= MinDW1) || (DB <= MinDB1 && DW < MinDW1) ) && 
	 RightSide(J2D,J2A,I1A,I1D,I2A,I2D) )
	JoinNeighbours(&Lnk1A,J1A,&Lnk1D,J1D,&Pat[i]->Nei1,Pat[j],&MinDB1,DB,&MinDW1,DW,&Min,j);

      if( Near(I1A,J1D,J1A,I1D,J2A,J2D,I2D,I2A,I1ACn,J1DCn,J1ACn,I1DCn,&DB,&DW) && 
	  ((DB < MinDB1 && DW <= MinDW1) || (DB <= MinDB1 && DW < MinDW1) ) && 
	 RightSide(J1A,J1D,I1A,I1D,I2A,I2D) )
	JoinNeighbours(&Lnk1A,J2A,&Lnk1D,J2D,&Pat[i]->Nei1,Pat[j],&MinDB1,DB,&MinDW1,DW,&Min,j);

      if( Near(I1A,J1A,J1D,I1D,J2D,J2A,I2D,I2A,I1ACn,J1ACn,J1DCn,I1DCn,&DB,&DW) && 
	  ((DB < MinDB1 && DW <= MinDW1) || (DB <= MinDB1 && DW < MinDW1) ) && 
	 RightSide(J1A,J1D,I1A,I1D,I2A,I2D) )
	JoinNeighbours(&Lnk1A,J2D,&Lnk1D,J2A,&Pat[i]->Nei1,Pat[j],&MinDB1,DB,&MinDW1,DW,&Min,j);

      if( Near(I1A,J2D,J2A,I1D,J1A,J1D,I2D,I2A,I1ACn,J2DCn,J2ACn,I1DCn,&DB,&DW) && 
	  ((DB < MinDB1 && DW <= MinDW1) || (DB <= MinDB1 && DW < MinDW1) ) && 
	 RightSide(J2D,J2A,I1A,I1D,I2A,I2D) )
	JoinNeighbours(&Lnk1A,J1A,&Lnk1D,J1D,&Pat[i]->Nei1,Pat[j],&MinDB1,DB,&MinDW1,DW,&Min,j);

      if( Near(I1A,J2A,J2D,I1D,J1D,J1A,I2D,I2A,I1ACn,J2ACn,J2DCn,I1DCn,&DB,&DW) && 
	  ((DB < MinDB1 && DW <= MinDW1) || (DB <= MinDB1 && DW < MinDW1) ) && 
	 RightSide(J2A,J2D,I1A,I1D,I2A,I2D) )
	JoinNeighbours(&Lnk1A,J1D,&Lnk1D,J2D,&Pat[i]->Nei1,Pat[j],&MinDB1,DB,&MinDW1,DW,&Min,j);

    }

    for( j=0; j<NPat; j++ ) {

      if( j == Min || j == i || !Pat[j]->ExistPattern ) continue;

      Alias(&J1D,&J1A,&J2D,&J2A,&J1DCn,&J1ACn,&J2DCn,&J2ACn,Pat[j]);

      if( Near(I2D,J1D,J1A,I2A,J2A,J2D,I1A,I1D,I2DCn,J1DCn,J1ACn,I2ACn,&DB,&DW) && 
	  (( DB < MinDB2 && DW <= MinDW2) || (DB <= MinDB2 && DW < MinDW2) ) &&
	 RightSide2(Lnk1A,Lnk1D,J2A,J2D,I1A,I1D,I2A,I2D) )
	JoinNeighb(&Pat[i]->Nei2,Pat[j],&MinDB2,DB,&MinDW2,DW);
      
      if( Near(I2D,J1A,J1D,I2A,J2D,J2A,I1A,I1D,I2DCn,J1ACn,J1DCn,I2ACn,&DB,&DW) && 
	  (( DB < MinDB2 && DW <= MinDW2) || (DB <= MinDB2 && DW < MinDW2) ) &&
	 RightSide2(Lnk1A,Lnk1D,J2D,J2A,I1A,I1D,I2A,I2D) )
	JoinNeighb(&Pat[i]->Nei2,Pat[j],&MinDB2,DB,&MinDW2,DW);
      
      if( Near(I2D,J2D,J2A,I2A,J1A,J1D,I1A,I1D,I2DCn,J2DCn,J2ACn,I2ACn,&DB,&DW) && 
	  (( DB < MinDB2 && DW <= MinDW2) || (DB <= MinDB2 && DW < MinDW2) ) &&
	 RightSide2(Lnk1A,Lnk1D,J1A,J1D,I1A,I1D,I2A,I2D) ) 
	JoinNeighb(&Pat[i]->Nei2,Pat[j],&MinDB2,DB,&MinDW2,DW);
      
      if( Near(I2D,J2A,J2D,I2A,J1D,J1A,I1A,I1D,I2DCn,J2ACn,J2DCn,I2ACn,&DB,&DW) && 
	  (( DB < MinDB2 && DW <= MinDW2) || (DB <= MinDB2 && DW < MinDW2) ) &&
	 RightSide2(Lnk1A,Lnk1D,J1D,J1A,I1A,I1D,I2A,I2D) )
	JoinNeighb(&Pat[i]->Nei2,Pat[j],&MinDB2,DB,&MinDW2,DW);
      
      if( Near(I2A,J1D,J1A,I2D,J2A,J2D,I1D,I1A,I2ACn,J1DCn,J1ACn,I2DCn,&DB,&DW) && 
	  (( DB < MinDB2 && DW <= MinDW2) || (DB <= MinDB2 && DW < MinDW2) ) &&
	 RightSide2(Lnk1A,Lnk1D,J2D,J2A,I1A,I1D,I2A,I2D) )
	JoinNeighb(&Pat[i]->Nei2,Pat[j],&MinDB2,DB,&MinDW2,DW);
      
      if( Near(I2A,J1A,J1D,I2D,J2D,J2A,I1D,I1A,I2ACn,J1ACn,J1DCn,I2DCn,&DB,&DW) && 
	  (( DB < MinDB2 && DW <= MinDW2) || (DB <= MinDB2 && DW < MinDW2) ) &&
	 RightSide2(Lnk1A,Lnk1D,J2A,J2D,I1A,I1D,I2A,I2D) )
	JoinNeighb(&Pat[i]->Nei2,Pat[j],&MinDB2,DB,&MinDW2,DW);
      
      if( Near(I2A,J2D,J2A,I2D,J1A,J1D,I1D,I1A,I2ACn,J2DCn,J2ACn,I2DCn,&DB,&DW) && 
	  (( DB < MinDB2 && DW <= MinDW2) || (DB <= MinDB2 && DW < MinDW2) ) &&
	 RightSide2(Lnk1A,Lnk1D,J1D,J1A,I1A,I1D,I2A,I2D) )
	JoinNeighb(&Pat[i]->Nei2,Pat[j],&MinDB2,DB,&MinDW2,DW);
      
      if( Near(I2A,J2A,J2D,I2D,J1D,J1A,I1D,I1A,I2ACn,J2ACn,J2DCn,I2DCn,&DB,&DW) && 
	  (( DB < MinDB2 && DW <= MinDW2) || (DB <= MinDB2 && DW < MinDW2) ) &&
	 RightSide2(Lnk1A,Lnk1D,J1A,J1D,I1A,I1D,I2A,I2D) )
	JoinNeighb(&Pat[i]->Nei2,Pat[j],&MinDB2,DB,&MinDW2,DW);
      }
  }
}

void MergePatternsPar(PATTERN **Pat, int NPat)
{
  register int i, j;
  int DB, DW, MinDB1, MinDB2, MinDW1, MinDW2, Min, Lnk1A, Lnk1D;
  int I1A, I1D, I2A, I2D, J1A, J1D, J2A, J2D;
  char I1ACn, I1DCn, I2ACn, I2DCn, J1ACn, J1DCn, J2ACn, J2DCn;

  for( i=0; i<NPat; i++ ) {

    if( !Pat[i]->ExistPattern ) continue;

    MinDB1 = MinDB2 = MinDW1 = MinDW2 = 1000;
    Min = ERR;
    Lnk1D = Lnk1A = ERR;

    Alias(&I1D,&I1A,&I2D,&I2A,&I1DCn,&I1ACn,&I2DCn,&I2ACn,Pat[i]);

    for( j=0; j<NPat; j++ ) {

      if( i == j || !Pat[j]->ExistPattern ) continue;

      Alias(&J1D,&J1A,&J2D,&J2A,&J1DCn,&J1ACn,&J2DCn,&J2ACn,Pat[j]);
      
      if( NearPar(I1D,J1D,J1A,I1A,J2A,J2D,I2A,I2D,I1DCn,J1DCn,J1ACn,I1ACn,&DB,&DW) && 
	  ((DB < MinDB1 && DW <= MinDW1) || (DB <= MinDB1 && DW < MinDW1) ) && 
	 RightSidePar(J1A,J1D,I1A,I1D,I2A,I2D) )
	JoinNeighbours(&Lnk1A,J2D,&Lnk1D,J2A,&Pat[i]->Nei1,Pat[j],&MinDB1,DB,&MinDW1,DW,&Min,j);

      if( NearPar(I1D,J1A,J1D,I1A,J2D,J2A,I2A,I2D,I1DCn,J1ACn,J1DCn,I1ACn,&DB,&DW) && 
	  ((DB < MinDB1 && DW <= MinDW1) || (DB <= MinDB1 && DW < MinDW1) ) && 
	 RightSidePar(J1D,J1A,I1A,I1D,I2A,I2D) )
	JoinNeighbours(&Lnk1A,J2A,&Lnk1D,J2D,&Pat[i]->Nei1,Pat[j],&MinDB1,DB,&MinDW1,DW,&Min,j);

      if( NearPar(I1D,J2D,J2A,I1A,J1A,J1D,I2A,I2D,I1DCn,J2DCn,J2ACn,I1ACn,&DB,&DW) && 
	  ((DB < MinDB1 && DW <= MinDW1) || (DB <= MinDB1 && DW < MinDW1) ) && 
	 RightSidePar(J2A,J2D,I1A,I1D,I2A,I2D) )
	JoinNeighbours(&Lnk1A,J1D,&Lnk1D,J1A,&Pat[i]->Nei1,Pat[j],&MinDB1,DB,&MinDW1,DW,&Min,j);

      if( NearPar(I1D,J2A,J2D,I1A,J1D,J1A,I2A,I2D,I1DCn,J2ACn,J2DCn,I1ACn,&DB,&DW) && 
	  ((DB < MinDB1 && DW <= MinDW1) || (DB <= MinDB1 && DW < MinDW1) ) && 
	 RightSidePar(J2D,J2A,I1A,I1D,I2A,I2D) )
	JoinNeighbours(&Lnk1A,J1A,&Lnk1D,J1D,&Pat[i]->Nei1,Pat[j],&MinDB1,DB,&MinDW1,DW,&Min,j);

      if( NearPar(I1A,J1D,J1A,I1D,J2A,J2D,I2D,I2A,I1ACn,J1DCn,J1ACn,I1DCn,&DB,&DW) && 
	  ((DB < MinDB1 && DW <= MinDW1) || (DB <= MinDB1 && DW < MinDW1) ) && 
	 RightSidePar(J1A,J1D,I1A,I1D,I2A,I2D) )
	JoinNeighbours(&Lnk1A,J2A,&Lnk1D,J2D,&Pat[i]->Nei1,Pat[j],&MinDB1,DB,&MinDW1,DW,&Min,j);

      if( NearPar(I1A,J1A,J1D,I1D,J2D,J2A,I2D,I2A,I1ACn,J1ACn,J1DCn,I1DCn,&DB,&DW) && 
	  ((DB < MinDB1 && DW <= MinDW1) || (DB <= MinDB1 && DW < MinDW1) ) && 
	 RightSidePar(J1A,J1D,I1A,I1D,I2A,I2D) )
	JoinNeighbours(&Lnk1A,J2D,&Lnk1D,J2A,&Pat[i]->Nei1,Pat[j],&MinDB1,DB,&MinDW1,DW,&Min,j);

      if( NearPar(I1A,J2D,J2A,I1D,J1A,J1D,I2D,I2A,I1ACn,J2DCn,J2ACn,I1DCn,&DB,&DW) && 
	  ((DB < MinDB1 && DW <= MinDW1) || (DB <= MinDB1 && DW < MinDW1) ) && 
	 RightSidePar(J2D,J2A,I1A,I1D,I2A,I2D) )
	JoinNeighbours(&Lnk1A,J1A,&Lnk1D,J1D,&Pat[i]->Nei1,Pat[j],&MinDB1,DB,&MinDW1,DW,&Min,j);

      if( NearPar(I1A,J2A,J2D,I1D,J1D,J1A,I2D,I2A,I1ACn,J2ACn,J2DCn,I1DCn,&DB,&DW) && 
	  ((DB < MinDB1 && DW <= MinDW1) || (DB <= MinDB1 && DW < MinDW1) ) && 
	 RightSidePar(J2A,J2D,I1A,I1D,I2A,I2D) )
	JoinNeighbours(&Lnk1A,J1D,&Lnk1D,J2D,&Pat[i]->Nei1,Pat[j],&MinDB1,DB,&MinDW1,DW,&Min,j);

    }

    for( j=0; j<NPat; j++ ) {

      if( j == Min || j == i || !Pat[j]->ExistPattern ) continue;

      Alias(&J1D,&J1A,&J2D,&J2A,&J1DCn,&J1ACn,&J2DCn,&J2ACn,Pat[j]);

      if( NearPar(I2D,J1D,J1A,I2A,J2A,J2D,I1A,I1D,I2DCn,J1DCn,J1ACn,I2ACn,&DB,&DW) && 
	  (( DB < MinDB2 && DW <= MinDW2) || (DB <= MinDB2 && DW < MinDW2) ) &&
	 RightSide2(Lnk1A,Lnk1D,J2A,J2D,I1A,I1D,I2A,I2D) )
	JoinNeighb(&Pat[i]->Nei2,Pat[j],&MinDB2,DB,&MinDW2,DW);
      
      if( NearPar(I2D,J1A,J1D,I2A,J2D,J2A,I1A,I1D,I2DCn,J1ACn,J1DCn,I2ACn,&DB,&DW) && 
	  (( DB < MinDB2 && DW <= MinDW2) || (DB <= MinDB2 && DW < MinDW2) ) &&
	 RightSide2(Lnk1A,Lnk1D,J2D,J2A,I1A,I1D,I2A,I2D) )
	JoinNeighb(&Pat[i]->Nei2,Pat[j],&MinDB2,DB,&MinDW2,DW);
      
      if( NearPar(I2D,J2D,J2A,I2A,J1A,J1D,I1A,I1D,I2DCn,J2DCn,J2ACn,I2ACn,&DB,&DW) && 
	  (( DB < MinDB2 && DW <= MinDW2) || (DB <= MinDB2 && DW < MinDW2) ) &&
	 RightSide2(Lnk1A,Lnk1D,J1A,J1D,I1A,I1D,I2A,I2D) ) 
	JoinNeighb(&Pat[i]->Nei2,Pat[j],&MinDB2,DB,&MinDW2,DW);
      
      if( NearPar(I2D,J2A,J2D,I2A,J1D,J1A,I1A,I1D,I2DCn,J2ACn,J2DCn,I2ACn,&DB,&DW) && 
	  (( DB < MinDB2 && DW <= MinDW2) || (DB <= MinDB2 && DW < MinDW2) ) &&
	 RightSide2(Lnk1A,Lnk1D,J1D,J1A,I1A,I1D,I2A,I2D) )
	JoinNeighb(&Pat[i]->Nei2,Pat[j],&MinDB2,DB,&MinDW2,DW);
      
      if( NearPar(I2A,J1D,J1A,I2D,J2A,J2D,I1D,I1A,I2ACn,J1DCn,J1ACn,I2DCn,&DB,&DW) && 
	  (( DB < MinDB2 && DW <= MinDW2) || (DB <= MinDB2 && DW < MinDW2) ) &&
	 RightSide2(Lnk1A,Lnk1D,J2D,J2A,I1A,I1D,I2A,I2D) )
	JoinNeighb(&Pat[i]->Nei2,Pat[j],&MinDB2,DB,&MinDW2,DW);
      
      if( NearPar(I2A,J1A,J1D,I2D,J2D,J2A,I1D,I1A,I2ACn,J1ACn,J1DCn,I2DCn,&DB,&DW) && 
	  (( DB < MinDB2 && DW <= MinDW2) || (DB <= MinDB2 && DW < MinDW2) ) &&
	 RightSide2(Lnk1A,Lnk1D,J2A,J2D,I1A,I1D,I2A,I2D) )
	JoinNeighb(&Pat[i]->Nei2,Pat[j],&MinDB2,DB,&MinDW2,DW);
      
      if( NearPar(I2A,J2D,J2A,I2D,J1A,J1D,I1D,I1A,I2ACn,J2DCn,J2ACn,I2DCn,&DB,&DW) && 
	  (( DB < MinDB2 && DW <= MinDW2) || (DB <= MinDB2 && DW < MinDW2) ) &&
	 RightSide2(Lnk1A,Lnk1D,J1D,J1A,I1A,I1D,I2A,I2D) )
	JoinNeighb(&Pat[i]->Nei2,Pat[j],&MinDB2,DB,&MinDW2,DW);
      
      if( NearPar(I2A,J2A,J2D,I2D,J1D,J1A,I1D,I1A,I2ACn,J2ACn,J2DCn,I2DCn,&DB,&DW) && 
	  (( DB < MinDB2 && DW <= MinDW2) || (DB <= MinDB2 && DW < MinDW2) ) &&
	 RightSide2(Lnk1A,Lnk1D,J1A,J1D,I1A,I1D,I2A,I2D) )
	JoinNeighb(&Pat[i]->Nei2,Pat[j],&MinDB2,DB,&MinDW2,DW);
      }
  }
}

int RightSide2(int L_A1, int L_D1, int LnkD, int LnkA, int I1A, int I1D, int I2A, int I2D)
{

  if( ( I2A < I1D && LnkA <= I1D && LnkA <= I2A ) || 
      ( I2A > I1D && LnkA >= I1D && LnkA >= I2A ) ||
      ( I2D < I1A && LnkD <= I1A && LnkA <= I2D ) || 
      ( I2D > I1A && LnkD >= I1A && LnkD >= I2D )  )
    return(SUCCESS);
  else
    if( I2A == I1D && I2D == I1A ) {
      if( L_A1 != ERR && 
	 ( ( LnkD <= I2D && L_A1 <= I2D && LnkA >= I2A && L_D1 >= I2A ) || 
	   ( LnkD >= I2D && L_A1 >= I2D && LnkA <= I2A && L_D1 <= I2A ) ) )
	return(FAILURE);
      else
	return(SUCCESS);
    }

  return(FAILURE);
}

int RightSide(int LnkA, int LnkD, int I1A, int I1D, int I2A, int I2D )
{

  if( ( I1A == I2D && I1D == I2A ) ||
      ( I1A < I2D && LnkA <= I2D && LnkA <= I1A ) || 
      ( I1A > I2D && LnkA >= I2D && LnkA >= I1A ) ||
      ( I1D < I2A && LnkD <= I2A && LnkD <= I1D ) ||
      ( I1D > I2A && LnkD >= I2A && LnkD >= I1D ) )
    return(SUCCESS);
  
  return(FAILURE);
}

int RightSidePar(int LnkA, int LnkD, int I1A, int I1D, int I2A, int I2D )
{

  if( ( I1A == I2D && I1D == I2A ) ||
      ( I1A < I2D && LnkA < I2D && LnkA <= I1A && I1D <= I2A && LnkD <= I2A && LnkD <= I1D ) || 
      ( I1A > I2D && LnkA > I2D && LnkA >= I1A && I1D >= I2A && LnkD >= I2A && LnkD >= I1D ) ||
      ( I1D < I2A && LnkD < I2A && LnkD <= I1D && I1A <= I2D && LnkA <= I2D && LnkA <= I1A ) ||
      ( I1D > I2A && LnkD > I2A && LnkD >= I1D && I1A >= I2D && LnkA >= I2D && LnkA >= I1A) )
    return(SUCCESS);
  
  return(FAILURE);
}

void JoinNeighbours(int *Lnk1A, int Res1, int *Lnk1D, int Res2, PATTERN **Nei, 
		    PATTERN *Pat, int *MinDB1, int DB, int *MinDW1, int DW, int *Min, int j)
{
  *Lnk1A = Res1;
  *Lnk1D = Res2;
  (*Nei) = Pat;
  *MinDB1 = DB;
  *MinDW1 = DW;
  *Min = j;
}

void JoinNeighb(PATTERN **Nei, PATTERN *Pat, int *MinDB2, int DB, int *MinDW2, int DW)
{
  (*Nei) = Pat;
  *MinDB2 = DB;
  *MinDW2 = DW;
}

int NearPar(int Res1, int Res2, int Res3, int Res4, int Res5, int Res6, int Res7, int Res8,
	 char Cn1, char Cn2, char Cn3, char Cn4, int *DistBest, int *DistWorst)
{

/*
   Res5 Res2 Res1
   Res6 Res3 Res4
*/

  int a, b, c1, d1, c, d, Nei1, Nei2;

  if( Cn1 != Cn2 || Cn3 != Cn4 ) return(FAILURE);

  if( Res1 >= Res2 && Res2 >= Res5 && Res7 >= Res1 && 
      Res4 >= Res3 && Res4 >= Res6 && Res8 >= Res4 ) {

    if( Res5 == Res2 ) 
      Nei1 = Res2;
    else
      Nei1 = Res2-1;

    if( Res1 == Res7 ) 
      Nei2 = Res1;
    else
      Nei2 = Res1+1;

    a = Nei2-Nei1;
    c1 = Nei2-Res5;

    if( Res3 == Res6 ) 
      Nei1 = Res3;
    else
      Nei1 = Res3-1;

    if( Res4 == Res8 ) 
      Nei2 = Res4;
    else
      Nei2 = Res4+1;

    b = Nei2-Nei1;
    d1 = Nei2-Res6;

  }
  else
  if( Res1 <= Res2 && Res2 <= Res5 && Res7 <= Res1 && 
      Res4 <= Res3 && Res4 <= Res6 && Res8 <= Res4 ) {

    if( Res5 == Res2 ) 
      Nei1 = Res2;
    else
      Nei1 = Res2+1;

    if( Res1 == Res7 ) 
      Nei2 = Res1;
    else
      Nei2 = Res1-1;

    a = Nei1-Nei2;
    c1 = Res1-Res7;

    if( Res3 == Res6 ) 
      Nei1 = Res3;
    else
      Nei1 = Res3+1;

    if( Res4 == Res8 ) 
      Nei2 = Res4;
    else
      Nei2 = Res4-1;

    b = Nei1-Nei2;
    d1 = Nei1-Res8;


  }
  else
   return(FAILURE);

  c = Maximum(c1,a);
  d = Maximum(d1,b);

  if( a >= 0 && b >= 0 && c >= 0 && d >= 0 && 
      ( (a <= 2 && b <= 5) || (a <= 5 && b <= 2) ) ) {
    *DistBest  = Minimum(a,b);
    *DistWorst = Maximum(c,d);
    if( *DistBest <= *DistWorst ) 
      return(SUCCESS);
    else
      return(FAILURE);
  }
  
  return(FAILURE);
}

int Near(int Res1, int Res2, int Res3, int Res4, int Res5, int Res6, int Res7, int Res8,
	 char Cn1, char Cn2, char Cn3, char Cn4, int *DistBest, int *DistWorst)
{

/*
   Res5 Res2 Res1
   Res6 Res3 Res4
*/

  int a, b, c1, d1, c, d, Nei1, Nei2;


  if( Cn1 != Cn2 || Cn3 != Cn4 ) return(FAILURE);


  if( Res1 >= Res2 && Res2 >= Res5 && Res7 >= Res1 && 
      Res4 <= Res3 && Res4 <= Res6 && Res8 <= Res4 ) {

    if( Res5 == Res2 ) 
      Nei1 = Res2;
    else
      Nei1 = Res2-1;
    
    if( Res1 == Res7 ) 
      Nei2 = Res1;
    else
      Nei2 = Res1+1;
    
    a = Nei2-Nei1;
    c1 = Nei2-Res5;
    
    if( Res3 == Res6 ) 
      Nei1 = Res3;
    else
      Nei1 = Res3+1;
    
    if( Res4 == Res8 ) 
      Nei2 = Res4;
    else
      Nei2 = Res4-1;
    
    b = Nei1-Nei2;
    d1 = Res6-Nei2;
  }
  else
    return(FAILURE);
  
  c = Maximum(c1,a);
  d = Maximum(d1,b);
  
  if( a >= 0 && b >= 0 && c >= 0 && d >= 0 && 
      ( (a <= 2 && b <= 5) || (a <= 5 && b <= 2) ) ) {
    *DistBest  = Minimum(a,b);
    *DistWorst = Maximum(c,d);
    if( *DistBest <= *DistWorst ) 
      return(SUCCESS);
    else
      return(FAILURE);
  }
  
  return(FAILURE);
}



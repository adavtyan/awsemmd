#include "stride.h"
    
void BetaTurn(CHAIN **Chain, int Cn)
{

  register int i;
  RESIDUE **r;
  TURN *t;
  int CA1, CA4, Tn;
  float Phi2, Phi3, Psi2, Psi3, Range1 = 30.0, Range2 = 45.0;
  char TurnType;

  for( i=0; i<Chain[Cn]->NRes-4; i++ ) {

    r = &Chain[Cn]->Rsd[i];
    
    if( r[1]->Prop->Asn == 'H' || r[2]->Prop->Asn == 'H' || 
        r[1]->Prop->Asn == 'G' || r[2]->Prop->Asn == 'G' || 
        r[1]->Prop->Asn == 'I' || r[2]->Prop->Asn == 'G' || 
       !FindAtom(Chain[Cn],i,"CA",&CA1) || !FindAtom(Chain[Cn],i+3,"CA",&CA4) ||
       Dist(r[0]->Coord[CA1],r[3]->Coord[CA4]) > 7.0 )
      continue;
    
    Phi2 = r[1]->Prop->Phi;
    Psi2 = r[1]->Prop->Psi;
    Phi3 = r[2]->Prop->Phi;
    Psi3 = r[2]->Prop->Psi;
    
    if( TurnCondition(Phi2,-60.0,Psi2,-30,Phi3,-90.0,Psi3,0,Range1,Range2) )
      TurnType = '1';
    else
    if( TurnCondition(Phi2,60.0,Psi2,30,Phi3,90.0,Psi3,0,Range1,Range2) )
      TurnType = '2';
    else
    if( TurnCondition(Phi2,-60.0,Psi2,120,Phi3,80.0,Psi3,0,Range1,Range2) )
      TurnType = '3';
    else
    if( TurnCondition(Phi2,60.0,Psi2,-120,Phi3,-80.0,Psi3,0,Range1,Range2) )
      TurnType = '4';
    else
    if( TurnCondition(Phi2,-60.0,Psi2,120,Phi3,-90.0,Psi3,0,Range1,Range2) )
      TurnType = '5';
    else
    if( TurnCondition(Phi2,-120.0,Psi2,120,Phi3,-60.0,Psi3,0,Range1,Range2) )
      TurnType = '6';
    else
    if( TurnCondition(Phi2,-60.0,Psi2,-30,Phi3,-120.0,Psi3,120,Range1,Range2) )
      TurnType = '7';
    else
      TurnType = '8';
  
    if( r[0]->Prop->Asn == 'C' ) 
      r[0]->Prop->Asn = 'T';
    
    if( r[1]->Prop->Asn == 'C' )
      r[1]->Prop->Asn = 'T';
    
    if( r[2]->Prop->Asn == 'C' )
      r[2]->Prop->Asn = 'T';
    
    if( r[3]->Prop->Asn == 'C' )
      r[3]->Prop->Asn = 'T';
    
    Tn = Chain[Cn]->NAssignedTurn;
    Chain[Cn]->AssignedTurn[Tn] = (TURN *)ckalloc(sizeof(TURN));
    t = Chain[Cn]->AssignedTurn[Tn];
    strcpy(t->Res1,r[0]->ResType);
    strcpy(t->Res2,r[3]->ResType);
    strcpy(t->PDB_ResNumb1,r[0]->PDB_ResNumb);
    strcpy(t->PDB_ResNumb2,r[3]->PDB_ResNumb);
    t->TurnType = TurnType;
    Chain[Cn]->NAssignedTurn++;

  }
}


void GammaTurn(CHAIN **Chain, int Cn, HBOND **HBond)
{

  register int i;
  RESIDUE **r;
  TURN *t;
  int Tn;
  float Phi2, Psi2;
  char TurnType, Asn;

  for( i=0; i<Chain[Cn]->NRes-2; i++ ) {

    r = &Chain[Cn]->Rsd[i-1];

    Asn = r[2]->Prop->Asn;

    if( Asn == 'H' || Asn == 'T' || Asn == 'G' || Asn == 'I' ||
        FindBnd(HBond,r[3],r[1]) == ERR ||
        (i > 0 && FindBnd(HBond,r[3],r[0]) != ERR) || 
        (i < Chain[Cn]->NRes-3 && FindBnd(HBond,r[4],r[1]) != ERR) )
      continue;
    
    Phi2 = r[2]->Prop->Phi;
    Psi2 = r[2]->Prop->Psi;
    
    if( Phi2 > 0.0 && Psi2 < 0.0 )
      TurnType = '@';
    else
    if( Phi2 < 0.0 && Psi2 > 0.0 )
      TurnType = '&';
    else 
      continue;

    if( r[1]->Prop->Asn == 'C' )
      r[1]->Prop->Asn = 'T';
    
    if( r[2]->Prop->Asn == 'C' )
      r[2]->Prop->Asn = 'T';
    
    if( r[3]->Prop->Asn == 'C' )
      r[3]->Prop->Asn = 'T';
    
    Tn = Chain[Cn]->NAssignedTurn;
    Chain[Cn]->AssignedTurn[Tn] = (TURN *)ckalloc(sizeof(TURN));
    t = Chain[Cn]->AssignedTurn[Tn];
    strcpy(t->Res1,r[1]->ResType);
    strcpy(t->Res2,r[3]->ResType);
    strcpy(t->PDB_ResNumb1,r[1]->PDB_ResNumb);
    strcpy(t->PDB_ResNumb2,r[3]->PDB_ResNumb);
    t->TurnType = TurnType;
    Chain[Cn]->NAssignedTurn++;
  }
}


int TurnCondition(float Phi2,float Phi2S,float Psi2,float Psi2S,
		  float Phi3,float Phi3S,float Psi3,float Psi3S,
		  float Range1,float Range2)
{
  if((IN(Phi2,Phi2S,Range2)==YES && IN(Psi2,Psi2S,Range1)==YES && 
      IN(Phi3,Phi3S,Range1)==YES && IN(Psi3,Psi3S,Range1)==YES)
     ||
     (IN(Phi2,Phi2S,Range1)==YES && IN(Psi2,Psi2S,Range2)==YES && 
      IN(Phi3,Phi3S,Range1)==YES && IN(Psi3,Psi3S,Range1)==YES)
     ||
     (IN(Phi2,Phi2S,Range1)==YES && IN(Psi2,Psi2S,Range1)==YES && 
      IN(Phi3,Phi3S,Range2)==YES && IN(Psi3,Psi3S,Range1)==YES)
     ||
     (IN(Phi2,Phi2S,Range1)==YES && IN(Psi2,Psi2S,Range1)==YES && 
      IN(Phi3,Phi3S,Range1)==YES && IN(Psi3,Psi3S,Range2)==YES)
     )
    return(SUCCESS);
  
  return(FAILURE);
}
    



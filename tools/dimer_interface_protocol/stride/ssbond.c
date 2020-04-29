#include "stride.h"
    
int SSBond(CHAIN **Chain, int NChain)
{

  register int Res1, Res2, Cn1, Cn2;
  int S1, S2, Bn, Cnt=0;
  
  for( Cn1=0; Cn1<NChain; Cn1++ )
    for( Res1=0; Res1<Chain[Cn1]->NRes; Res1++ ) {
      if( strcmp(Chain[Cn1]->Rsd[Res1]->ResType,"CYS") )
	continue;
      for( Cn2=Cn1; Cn2<NChain; Cn2++ )
	for( Res2 = ( (Cn2 == Cn1)? Res1+1 : 0) ; Res2<Chain[Cn2]->NRes; Res2++ ) {
	  if( strcmp(Chain[Cn2]->Rsd[Res2]->ResType,"CYS") )
	    continue;
	  
	  if( !ExistSSBond(Chain,NChain,Cn1,Cn2,
			   Chain[Cn1]->Rsd[Res1]->PDB_ResNumb,
			   Chain[Cn2]->Rsd[Res2]->PDB_ResNumb) && 
	     FindAtom(Chain[Cn1],Res1,"SG",&S1) && FindAtom(Chain[Cn2],Res2,"SG",&S2) &&
	     Dist(Chain[Cn1]->Rsd[Res1]->Coord[S1],
		  Chain[Cn2]->Rsd[Res2]->Coord[S2]) <= SSDIST ) {
	    Bn = Chain[0]->NBond;
	    Chain[0]->SSbond[Bn] =  (SSBOND *)ckalloc(sizeof(SSBOND));
	    strcpy(Chain[0]->SSbond[Bn]->PDB_ResNumb1,Chain[Cn1]->Rsd[Res1]->PDB_ResNumb);
	    strcpy(Chain[0]->SSbond[Bn]->PDB_ResNumb2,Chain[Cn2]->Rsd[Res2]->PDB_ResNumb);
	    Chain[0]->SSbond[Bn]->ChainId1 = Chain[Cn1]->Id;
	    Chain[0]->SSbond[Bn]->ChainId2 = Chain[Cn2]->Id;
	    Chain[0]->SSbond[Bn]->AsnSource = Stride;
	    Chain[0]->NBond++;
	    Cnt++;
	  }
	}
    }
  
  return(Cnt);
}

BOOLEAN ExistSSBond(CHAIN **Chain,int NChain, int Cn1,int Cn2,char *Res1,char *Res2)
{

  register int i;
  SSBOND *ptr;

  for( i=0; i<Chain[0]->NBond; i++ ) {
    ptr = Chain[0]->SSbond[i];
    if( ( !strcmp(Res1,ptr->PDB_ResNumb1) && 
	  !strcmp(Res2,ptr->PDB_ResNumb2) &&
	  FindChain(Chain,NChain,ptr->ChainId1) == Cn1 &&
	  FindChain(Chain,NChain,ptr->ChainId2) == Cn2 ) ||
        ( !strcmp(Res2,ptr->PDB_ResNumb1) && 
          !strcmp(Res1,ptr->PDB_ResNumb2) &&
	  FindChain(Chain,NChain,ptr->ChainId1) == Cn2 &&
	  FindChain(Chain,NChain,ptr->ChainId2) == Cn1 ) )
      return(SUCCESS);
  }

  return(FAILURE);
}

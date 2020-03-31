#include "stride.h"

void BackboneAngles(CHAIN **Chain, int NChain)
{

  register int Res, Cn;

  for( Cn=0; Cn<NChain; Cn++ ) {

    for( Res=0; Res<Chain[Cn]->NRes; Res++ ) {
      PHI(Chain[Cn],Res);
      PSI(Chain[Cn],Res);
    }
  }
}

void DiscrPhiPsi(CHAIN **Chain, int NChain, COMMAND *Cmd)
{

  register int i, Res, Cn;
  RESIDUE *r;

  for( Cn=0; Cn<NChain; Cn++ ) {
    
    for( Res=0; Res<Chain[Cn]->NRes; Res++ ) {
      
      r = Chain[Cn]->Rsd[Res];
      
      r->Prop->PhiZn = ERR;
      r->Prop->PsiZn = ERR;
      
      if( Res != 0 ) {
	for( i=0; i<Cmd->NPixel; i++ )
	  if( r->Prop->Phi  >  MINPHIPSI+(float)(i)*Cmd->PhiPsiStep && 
	      r->Prop->Phi <=  MINPHIPSI+(float)(i+1)*Cmd->PhiPsiStep ) {
	    r->Prop->PhiZn = i;
	    break;
	  }
      }
      
      if( Res != Chain[Cn]->NRes-1 ) {
	for( i=0; i<Cmd->NPixel; i++ )
	  if( r->Prop->Psi  >  MINPHIPSI+(float)(i)*Cmd->PhiPsiStep && 
	      r->Prop->Psi <=  MINPHIPSI+(float)(i+1)*Cmd->PhiPsiStep ) {
	    r->Prop->PsiZn = i;
	    break;
	  }
      }
      
    }
    
    for(Res=0; Res<Chain[Cn]->NRes; Res++ )  {
      r = Chain[Cn]->Rsd[Res];
      if( Res != 0 && r->Prop->PsiZn == ERR )
	r->Prop->PsiZn = Chain[Cn]->Rsd[Res-1]->Prop->PsiZn;
      if( Res != Chain[Cn]->NRes-1 && r->Prop->PhiZn == ERR )
	r->Prop->PhiZn = Chain[Cn]->Rsd[Res+1]->Prop->PhiZn;
    }

  }
}



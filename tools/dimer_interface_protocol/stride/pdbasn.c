#include "stride.h"

/*************************************************************************
**                                                                      **
** Get PDB secondary structure assignment for every residue             **
**                                                                      **
*************************************************************************/
void GetPdbAsn(CHAIN **Chain, int NChain)
{

  register int i, j, k;
  int Cn, Beg, End;
  char SecondStr;
  CHAIN *c;

  for( Cn=0; Cn<NChain; Cn++ ) {

    c = Chain[Cn];

    for( i=0; i<c->NHelix; i++ ) {

      switch( c->Helix[i]->Class ) {
      case 1:  SecondStr = 'H';
	break;
      case 3:  SecondStr = 'I';
	break;
      case 5:  SecondStr = 'G';
	break;
      }

      if( PdbN2SeqN(c,c->Helix[i]->PDB_ResNumb1,&Beg) &&
	  PdbN2SeqN(c,c->Helix[i]->PDB_ResNumb2,&End) )
	for( j=Beg; j<=End; j++ )
	  if( c->Rsd[j]->Prop->PdbAsn != 'H' )
	    c->Rsd[j]->Prop->PdbAsn = SecondStr;
    }
    
    for( i=0; i<c->NSheet; i++ )
      for( j=0; j<c->Sheet[i]->NStrand; j++ ) {
	  if( PdbN2SeqN(c,c->Sheet[i]->PDB_ResNumb1[j],&Beg) &&
	      PdbN2SeqN(c,c->Sheet[i]->PDB_ResNumb2[j],&End) )
	    for( k=Beg; k<=End; k++ )
	      if( c->Rsd[k]->Prop->PdbAsn != 'H' )
		c->Rsd[k]->Prop->PdbAsn = 'E';
	}

    for( i=0; i<c->NTurn; i++ ) { 
      if( PdbN2SeqN(c,c->Turn[i]->PDB_ResNumb1,&Beg) &&
	  PdbN2SeqN(c,c->Turn[i]->PDB_ResNumb2,&End) )
	for( j=Beg; j<=End; j++ )
	  if( c->Rsd[j]->Prop->PdbAsn != 'H' &&  c->Rsd[j]->Prop->PdbAsn != 'E' )
	    c->Rsd[j]->Prop->PdbAsn = 'T';
    }
  }
}




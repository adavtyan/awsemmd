#include "stride.h"
    
void ContactMap(CHAIN **Chain, int NChain, COMMAND *Cmd )
{

  register int Res1, Res2;
  int Cn1, Cn2, CA1, CA2;
  float D;
  
  for( Cn1=0; Cn1<NChain; Cn1++ ) {
    
    if( !Chain[Cn1]->Valid )
      continue;


    for( Cn2=0; Cn2<NChain; Cn2++ ) {

      if( !Chain[Cn2]->Valid )
	continue;

      for( Res1=0; Res1<Chain[Cn1]->NRes; Res1++ ) {

	if( !FindAtom(Chain[Cn1],Res1,"CA",&CA1) )
	  continue;

	for( Res2=0; Res2<Chain[Cn2]->NRes; Res2++ ) {
	  
	  if( !FindAtom(Chain[Cn2],Res2,"CA",&CA2) )
	    continue;

	  D = Dist(Chain[Cn1]->Rsd[Res1]->Coord[CA1],
		   Chain[Cn2]->Rsd[Res2]->Coord[CA2]);
	  printf("%c %d <-> %c %d : %7.3f\n",
		 SpaceToDash(Chain[Cn1]->Id),Res1,
		 SpaceToDash(Chain[Cn2]->Id),Res2,D);
	}
      }
    }
  }
  exit(0);
}

#include "stride.h"

int PlaceHydrogens(CHAIN *Chain)
{

  int Res, i, N, C, CA, H, PlacedCnt=0;
  float Length_N_C, Length_N_CA, Length_N_H;
  RESIDUE *r, *rr;
  
  for( Res=1; Res<Chain->NRes; Res++ ) {
    
    r  = Chain->Rsd[Res];
    rr = Chain->Rsd[Res-1];

      if( !strcmp(r->ResType,"PRO") ) continue;

      /* Replace deiterium atoms by hydrogens */
      if( FindAtom(Chain,Res,"D",&H) ) 
	strcmp(r->AtomType[H],"H");

      if( !FindAtom(Chain,Res,"H",&H)  && FindAtom(Chain,Res,"N",&N)   &&
	   FindAtom(Chain,Res-1,"C",&C) && FindAtom(Chain,Res,"CA",&CA) ) {

	H = r->NAtom;
	
	Length_N_C   = Dist(r->Coord[N],rr->Coord[C]);
	Length_N_CA  = Dist(r->Coord[N],r->Coord[CA]);
	
	for( i=0; i<3; i++ )
	  r->Coord[H][i] = r->Coord[N][i] - 
	      ( (rr->Coord[C][i] -  r->Coord[N][i])/Length_N_C +
	        (r->Coord[CA][i]  -  r->Coord[N][i])/Length_N_CA );

	Length_N_H = Dist(r->Coord[N],r->Coord[H]);

	for( i=0; i<3; i++ )
	  r->Coord[H][i] = r->Coord[N][i] + 
	    DIST_N_H*(r->Coord[H][i]-r->Coord[N][i])/Length_N_H;

	strcpy(r->AtomType[H],"H");
	r->NAtom++;
	PlacedCnt++;
      }
    }
  return(PlacedCnt);
}




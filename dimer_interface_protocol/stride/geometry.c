/*************************************************************************
**                                                                      **
**  Calculate torsion angle                                             **
**                                                                      **
** INPUT:  *Coord1, *Coord2, *Coord3, *Coord4 Coordinates of four atoms **
**                                                                      **
** RETURNS: Calculate torsion angle                                     **
**                                                                      **
** Adapted  from  the  program  of D.S.Moss.                            **
** Reference:   Moss,  D. S. (1992)  Molecular  geometry.  In:          **
** Computer modelling of biomolecular processess,  Goodfellow, J.M,     **
** Moss,D.S., eds, pp. 5-18.                                            **
**                                                                      **
*************************************************************************/

#include "stride.h"

float Torsion(float *Coord1, float *Coord2, float *Coord3, float *Coord4)
{
  double Comp[3][3], ScalarProd, TripleScalarProd, AbsTorsAng;
  double Perp_123[3], Perp_234[3], Len_Perp_123, Len_Perp_234;
  int i, j, k;

  /* Find the components of the three bond vectors */
  for( i=0; i<3; i++ ) {
    Comp[0][i] = (double)(Coord2[i]-Coord1[i]);
    Comp[1][i] = (double)(Coord3[i]-Coord2[i]);
    Comp[2][i] = (double)(Coord4[i]-Coord3[i]);
  }

  /* Calculate vectors perpendicular to the planes 123 and 234 */
  Len_Perp_123 = 0.0; Len_Perp_234 = 0.0;
  for( i=0; i<3; i++ ) {
    j = (i+1)%3;
    k = (j+1)%3;
    Perp_123[i] = Comp[0][j]*Comp[1][k] - Comp[0][k]*Comp[1][j];
    Perp_234[i] = Comp[1][j]*Comp[2][k] - Comp[1][k]*Comp[2][j];
    Len_Perp_123 += Perp_123[i]*Perp_123[i];
    Len_Perp_234 += Perp_234[i]*Perp_234[i];
  }
  
  Len_Perp_123 = sqrt(Len_Perp_123);
  Len_Perp_234 = sqrt(Len_Perp_234);

  /* Normalize the vectors perpendicular to 123 and 234 */
  for( i=0; i<3; i++ ) {
    Perp_123[i] /= Len_Perp_123;
    Perp_234[i] /= Len_Perp_234;
  }

  /* Find the scalar product of the unit normals */
  ScalarProd = 0.0;
  for( i=0; i<3; i++ ) 
    ScalarProd += Perp_123[i]*Perp_234[i];

  /* Find the absolute value of the torsion angle */
  if( ScalarProd > 0.0 && fabs(ScalarProd - 1.0) < Eps ) 
    ScalarProd -= Eps;
  else
  if( ScalarProd < 0.0 && fabs(ScalarProd + 1.0) < Eps ) 
    ScalarProd += Eps;
  AbsTorsAng = RADDEG*acos(ScalarProd);

  /* Find the triple scalar product of the three bond vectors */
  TripleScalarProd = 0.0;
  for( i=0; i<3; i++ ) 
    TripleScalarProd += Comp[0][i]*Perp_234[i];

  /* Torsion angle has the sign of the triple scalar product */
  return( (TripleScalarProd > 0.0) ? (float)AbsTorsAng : (float)(-AbsTorsAng) );

}
/*************************************************************************
**                                                                      **
** INPUT:  *Coord1   Coordinates of the first point                     **
**         *Coord2   Coordinates of the second point                    **
**                                                                      **
** RETURNS:          Distance between two points                        **
**                                                                      **
*************************************************************************/
float Dist(float *Coord1, float *Coord2)
{
  register int i;
  float Dist=0;

  for( i=0; i<3; i++ ) 
    Dist += (Coord1[i]-Coord2[i])*(Coord1[i]-Coord2[i]);

  return( sqrt(Dist) );
}

/*************************************************************************
**                                                                      **
** INPUT:  *Coord1   Coordinates of the first point                     **
**         *Coord2   Coordinates of the second point                    **
**         *Coord3   Coordinates of the third point                     **
**                                                                      **
** RETURNS:          Angle 1-2-3                                        **
**                                                                      **
*************************************************************************/
float Ang(float *Coord1, float *Coord2, float *Coord3)
{
  float Vector1[3], Vector2[3];
  double A, B, C, D;

  Vector1[0] = Coord1[0] - Coord2[0];
  Vector1[1] = Coord1[1] - Coord2[1];
  Vector1[2] = Coord1[2] - Coord2[2];

  Vector2[0] = Coord3[0] - Coord2[0];
  Vector2[1] = Coord3[1] - Coord2[1];
  Vector2[2] = Coord3[2] - Coord2[2];

  A = Vector1[0]*Vector2[0]+Vector1[1]*Vector2[1]+Vector1[2]*Vector2[2];
  B = sqrt( Vector1[0]*Vector1[0]+Vector1[1]*Vector1[1]+Vector1[2]*Vector1[2]);
  C = sqrt( Vector2[0]*Vector2[0]+Vector2[1]*Vector2[1]+Vector2[2]*Vector2[2]);

  D = A/(B*C);
  if( D > 0.0 && fabs(D - 1.0) < Eps ) 
    D -= Eps;
  else
  if( D < 0.0 && fabs(D + 1.0) < Eps ) 
    D += Eps;
  
  return((float)(RADDEG*acos(D)));
}

/*************************************************************************
**                                                                      **
** INPUT:  *Chain    Protein chain                                      **
**         *Res      Residue number                                     **
**                                                                      **
** OUTPUT: Chain->Rsd[Res]->Prop->Phi Phi torsional angle               **
**                                                                      **
*************************************************************************/
void PHI(CHAIN *Chain, int Res)
{

  int C_Prev, N_Curr, CA_Curr, C_Curr;
  RESIDUE *r, *rr;

  r = Chain->Rsd[Res];
  r->Prop->Phi = 360.0;

  if( Res == 0 )
    return;

  rr = Chain->Rsd[Res-1];

  if( FindAtom(Chain,Res-1,"C",&C_Prev) && FindAtom(Chain,Res,"N",&N_Curr)   &&
      FindAtom(Chain,Res,"CA",&CA_Curr) && FindAtom(Chain,Res,"C",&C_Curr)   &&
      Dist(rr->Coord[C_Prev],r->Coord[N_Curr]) < BREAKDIST ) {
    r->Prop->Phi = Torsion(rr->Coord[C_Prev],r->Coord[N_Curr],
			   r->Coord[CA_Curr],r->Coord[C_Curr]);
  }
}

/*************************************************************************
**                                                                      **
** INPUT:  *Chain    Protein chain                                      **
**         *Res      Residue number                                     **
**                                                                      **
** OUTPUT: Chain->Rsd[Res]->Prop->Psi  Psi torsional angle              **
**                                                                      **
*************************************************************************/
void PSI(CHAIN *Chain, int Res)
{

  int N_Curr, CA_Curr, C_Curr, N_Next;
  RESIDUE *r, *rr;

  r = Chain->Rsd[Res];
  r->Prop->Psi = 360.0;

  if( Res == Chain->NRes-1 )
    return;

  rr = Chain->Rsd[Res+1];

  if( FindAtom(Chain,Res,"N",&N_Curr) && FindAtom(Chain,Res,"CA",&CA_Curr) &&
      FindAtom(Chain,Res,"C",&C_Curr) && FindAtom(Chain,Res+1,"N",&N_Next) &&
      Dist(r->Coord[C_Curr],rr->Coord[N_Next]) < BREAKDIST ){

    r->Prop->Psi = Torsion(r->Coord[N_Curr],r->Coord[CA_Curr],
			   r->Coord[C_Curr],rr->Coord[N_Next]);
  }
}

/*************************************************************************
**                                                                      **
** INPUT:  *Chain    Protein chain                                      **
**         *Res      Residue number                                     **
**                                                                      **
** OUTPUT: *Omega    Omega torsional angle                              **
**                                                                      **
*************************************************************************/
void OMEGA(CHAIN *Chain, int Res)
{

  int CA_Prev, C_Prev, N_Curr, CA_Curr;
  RESIDUE *r, *rr;

  
  r = Chain->Rsd[Res];
  r->Prop->Omega = 360.0;

  if( Res == 0 )
    return;

  rr = Chain->Rsd[Res-1];

  if( FindAtom(Chain,Res-1,"CA",&CA_Prev) && FindAtom(Chain,Res-1,"C",&C_Prev)   &&
      FindAtom(Chain,Res,"N",&N_Curr)     && FindAtom(Chain,Res,"CA",&CA_Curr) ) {

    r->Prop->Omega = Torsion(rr->Coord[CA_Prev],rr->Coord[C_Prev],
			     r->Coord[N_Curr],r->Coord[CA_Curr]);
  }
}

/*************************************************************************
**                                                                      **
** Place atom X in the plane of atoms 1,2 and 3 given the               **
** distance |X-3| and angle 2-3-X                                       **
**                                                                      **
** INPUT: *Coord1, *Coord2, *Coord3  Coordinates of three atoms  in the **
**                                   plane                              **
**        Dist3X                     Distance between atom 3 and the    **
**                                   atom to be placed                  **
**        Ang23X                     Angle between atoms 2,3 and the    **
**                                   atom to be placed                  **
**                                                                      **
** OUTPUT: *Coordx                   Coordinates of the placed atom     **
**                                                                      **
*************************************************************************/


void Place123_X(float *Coord1, float *Coord2, float *Coord3, float Dist3X, float Ang23X, 
	       float *CoordX)
{

/*

   Atom1
     \                     AtomX
      \ ^UnVect2            /
       \|                  /
       Atom2----------Atom3->UnVect1

*/

  float Length_23, Length_12;
  float Proj3X_1, Proj3X_2, Proj12_1, Proj12_2, Rad1, Rad2;
  float UnVect1[3], UnVect2[3];
  int i;

  Length_23 = Dist(Coord3,Coord2);
  Length_12 = Dist(Coord2,Coord1);
  Rad1 = RAD(180.0-Ang23X);
  Rad2 = RAD(Ang(Coord1,Coord2,Coord3)-90.0); 
  Proj3X_1 = Dist3X*cos(Rad1);
  Proj3X_2 = Dist3X*sin(Rad1);
  Proj12_2 = cos(Rad2)*Length_12;
  Proj12_1 = sin(Rad2)*Length_12;

  for( i=0; i<3; i++ ) {
    UnVect1[i] = (Coord3[i]-Coord2[i])/Length_23;
    UnVect2[i] = ((Coord1[i]-Coord2[i]) - ( -UnVect1[i]*Proj12_1))/Proj12_2;
  }

  for( i=0; i<3; i++ ) 
    CoordX[i] = Proj3X_1*UnVect1[i]+Proj3X_2*UnVect2[i]+Coord3[i];
}


/*************************************************************************
**                                                                      **
** INPUT: *Vector1, Vector2                                             **
**                                                                      **
** OUTPUT: *Product          Vector pruduct of Vector1 and Vector2      **
**                                                                      **
*************************************************************************/
float VectorProduct(float *Vector1, float *Vector2, float *Product)
{

  int i, j, k;
  float ProductLength;

  ProductLength = 0.0; 

  for( i=0; i<3; i++ ) {
    j = (i+1)%3;
    k = (j+1)%3;
    Product[i] = Vector1[j]*Vector2[k] - Vector1[k]*Vector2[j];
    ProductLength += Product[i]*Product[i];
  }

  return(sqrt(ProductLength));
}

/*************************************************************************
**                                                                      **
** Find projection of an atom to a plane                                **
**                                                                      **
** INPUT: *Coord1, *Coord2, *Coord3  Coordinates of three atoms in a    **
**                                   plance                             **
**        *Coord4                    Coordinates of the fourth atom     **
**                                                                      **
** OUTPUT: *Coord_Proj4_123          Coordinates of the fourth atom's   **
**                                   projection to the place            **
**                                                                      **
*************************************************************************/
void Project4_123(float *Coord1, float *Coord2, float *Coord3, float *Coord4, 
		 float *Coord_Proj4_123)
{

/*
                          Atom4
   Atom3                  .
   \                     .
    \                   .    
     \                 .  .Proj4_123
      \               .. 
      Atom2-------Atom1

*/


  float Vector21[3], Vector23[3], Vector14[3], VectorNormal_123[3];
  float Length_21 = 0.0, Length_23 = 0.0, Length_14 = 0.0, NormalLength;
  float COS_Norm_14, Proj_14_Norm;
  int i;

  for( i=0; i<3; i++ ) {
    Vector21[i] = Coord1[i] - Coord2[i];
    Vector23[i] = Coord3[i] - Coord2[i];
    Vector14[i] = Coord4[i] - Coord1[i];
    Length_21 += Vector21[i]*Vector21[i];
    Length_23 += Vector23[i]*Vector23[i];
    Length_14 += Vector14[i]*Vector14[i];
  }
  
  Length_21 = sqrt(Length_21);
  Length_23 = sqrt(Length_23);
  Length_14 = sqrt(Length_14);

  NormalLength = VectorProduct(Vector21,Vector23,VectorNormal_123);

  for( i=0; i<3; i++ ) 
    VectorNormal_123[i] /= NormalLength;

  COS_Norm_14 = 0.0;

  for( i=0; i<3; i++ ) 
    COS_Norm_14 += VectorNormal_123[i]*Vector14[i];

  COS_Norm_14 /= (Length_14*NormalLength);

  if( COS_Norm_14 < 0.0 ) {
    COS_Norm_14 = fabs(COS_Norm_14);
    for( i=0; i<3; i++ ) 
      VectorNormal_123[i] = -VectorNormal_123[i];
  }

  Proj_14_Norm = Length_14*COS_Norm_14;

  for( i=0; i<3; i++ ) {
    VectorNormal_123[i] *= Proj_14_Norm;
    Coord_Proj4_123[i] = (Vector14[i] - VectorNormal_123[i]) + Coord1[i];
  }
}

double GetAtomRadius(char *AtomType)
{

  if( !strcmp(AtomType,"O") )
    return(1.40);
  else
  if( !strcmp(AtomType,"N") )
    return(1.65);
  else
  if( !strcmp(AtomType,"CA") )
    return(1.87);
  else
  if( !strcmp(AtomType,"C") )
    return(1.76);
  else
    return(1.80);
}


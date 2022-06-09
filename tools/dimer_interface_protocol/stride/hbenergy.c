#include "stride.h"

/*************************************************
 Calculate the hydrogen bond energy as defined by
 Boobbyer et al., 1989
**************************************************/

void GRID_Energy(float *CA2, float *C, float *O, float *H, float *N, COMMAND *Cmd, HBOND *HBond)
{

  float ProjH[3];

 /***** Distance dependence ( 8-6 potential ) ****/

  if( Cmd->Truncate && HBond->AccDonDist < RmGRID ) 
    HBond->AccDonDist = RmGRID;
  HBond->Er = CGRID/pow(HBond->AccDonDist,8.0) - DGRID/pow(HBond->AccDonDist,6.0);

 /************** Angular dependance ****************/

 /* Find projection of the hydrogen on the O-C-CA plane */
  Project4_123(O,C,CA2,H,ProjH); 


 /* Three angles determining the direction of the hydrogen bond */
  HBond->ti = fabs(180.0 - Ang(ProjH,O,C)); 
  HBond->to = Ang(H,O,ProjH);             
  HBond->p  = Ang(N,H,O);

 /* Calculate both angle-dependent HB energy components Et and Ep */ 
  if( HBond->ti >= 0.0 && HBond->ti < 90.0 )   
    HBond->Et = cos(RAD(HBond->to))*(0.9+0.1*sin(RAD(2*HBond->ti)));
  else
  if( HBond->ti >= 90.0 && HBond->ti < 110.0 ) 
    HBond->Et = K1GRID*cos(RAD(HBond->to))*
      (pow((K2GRID-pow(cos(RAD(HBond->ti)),2.0)),3.0));
  else
    HBond->Et = 0.0;

  if( HBond->p > 90.0 && HBond->p < 270.0 )
    HBond->Ep = pow(cos(RAD(HBond->p)),2.0);
  else
    HBond->Ep = 0.0;

    /******** Full hydrogen bond energy *********************/
  HBond->Energy = 1000.0*HBond->Er*HBond->Et*HBond->Ep;
}

#define Q -27888.0

/********************************************************
 Calculate the energy of polar interaction as defined by
 Kabsch and Sander (1983) 
*********************************************************/

void DSSP_Energy(float *Dummy, float *C, float *O, float *H, float *N, COMMAND *Cmd, 
		 HBOND *HBond)
                     
/* Dummy not used, for compatibility with GRID_Energy */
{ HBond->Energy = Q/Dist(O,H) + Q/Dist(C,N) - Q/HBond->AccDonDist - Q/Dist(C,H); }













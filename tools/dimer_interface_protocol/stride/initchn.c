#include "stride.h"


void InitChain(CHAIN **Chain)
{
  
  *Chain = (CHAIN *)ckalloc(sizeof(CHAIN));

  (*Chain)->NRes                = 0;
  (*Chain)->NHetRes             = 0;
  (*Chain)->NonStandRes         = 0;
  (*Chain)->NHet                = 0;
  (*Chain)->NonStandAtom        = 0;
  (*Chain)->NHelix              = 0;
  (*Chain)->NSheet              = -1;
  (*Chain)->NTurn               = 0;
  (*Chain)->NAssignedTurn       = 0;
  (*Chain)->NBond               = 0;
  (*Chain)->NHydrBond           = 0;
  (*Chain)->NHydrBondTotal      = 0;
  (*Chain)->NHydrBondInterchain = 0;
  (*Chain)->Method              = XRay;
  (*Chain)->Ter                 = 0;
  (*Chain)->Resolution          = 0.0;

  (*Chain)->File                = (char           *)ckalloc(BUFSZ*sizeof(char));
  (*Chain)->Rsd                 = (RESIDUE       **)ckalloc(MAX_RES*sizeof(RESIDUE *));
  (*Chain)->HetRsd              = (HETERORESIDUE **)ckalloc(MAX_HETRES*sizeof(HETERORESIDUE *));
  (*Chain)->Het                 = (HET           **)ckalloc(MAX_HET*sizeof(HET *));
  (*Chain)->Helix               = (HELIX         **)ckalloc(MAX_HELIX*sizeof(HELIX *));
  (*Chain)->Sheet               = (SHEET         **)ckalloc(MAX_SHEET*sizeof(SHEET *));
  (*Chain)->Turn                = (TURN          **)ckalloc(MAX_TURN*sizeof(TURN *));
  (*Chain)->AssignedTurn        = (TURN          **)ckalloc(MAX_TURN*sizeof(TURN *));
  (*Chain)->SSbond              = (SSBOND        **)ckalloc(MAX_BOND*sizeof(SSBOND *));
  (*Chain)->Info                = (char          **)ckalloc(MAX_INFO*sizeof(char *));

  (*Chain)->Valid               = YES;
}


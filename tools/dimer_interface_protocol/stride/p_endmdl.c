#include "stride.h"

int Process_ENDMDL(BUFFER Buffer, CHAIN **Chain, int *ChainNumber)
{

  int CC;

  for( CC=0; CC < *ChainNumber; CC++ )
    Chain[CC]->Ter = 1;
  
  return(SUCCESS);
}

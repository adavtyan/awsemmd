#include "stride.h"

int Process_HELIX(BUFFER Buffer, CHAIN **Chain, int *ChainNumber, COMMAND *Cmd)
{
  int CC, HC;
  char *Field[MAX_FIELD];
  BUFFER Tmp;

  if( Cmd->NActive && !ChInStr(Cmd->Active,SpaceToDash(Buffer[19])) )
     return(SUCCESS);
     
  for( CC=0; CC < *ChainNumber && Chain[CC]->Id != Buffer[19]; CC++ );

  if( CC == *ChainNumber ) {
    InitChain(&Chain[CC]); 
    Chain[CC]->Id = Buffer[19];
    (*ChainNumber)++;
  }

  HC = Chain[CC]->NHelix;
  Chain[CC]->Helix[HC] =  (HELIX *)ckalloc(sizeof(HELIX));

  SplitString(Buffer+15,Field,1);

  strcpy(Chain[CC]->Helix[HC]->Res1,Field[0]);

  SplitString(Buffer+27,Field,1);

  strcpy(Chain[CC]->Helix[HC]->Res2,Field[0]);

  strcpy(Tmp,Buffer);
  Tmp[25] = ' ';
  Tmp[37] = ' ';
  SplitString(Tmp+21,Field,1);
  strcpy(Chain[CC]->Helix[HC]->PDB_ResNumb1,Field[0]);
  SplitString(Tmp+33,Field,1);
  strcpy(Chain[CC]->Helix[HC]->PDB_ResNumb2,Field[0]);

  Chain[CC]->Helix[HC]->InsCode1 = Buffer[25];
  Chain[CC]->Helix[HC]->InsCode2 = Buffer[37];

  Tmp[40] = ' ';
  SplitString(Tmp+38,Field,1);
  Chain[CC]->Helix[HC]->Class = atoi(Field[0]);

  Chain[CC]->NHelix++;

  return(SUCCESS);
}


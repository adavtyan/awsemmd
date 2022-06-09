#include "stride.h"

int Process_TURN(BUFFER Buffer, CHAIN **Chain, int *ChainNumber, COMMAND *Cmd)
{
  int CC, TC;
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

  TC = Chain[CC]->NTurn;
  Chain[CC]->Turn[TC] = (TURN *)ckalloc(sizeof(TURN));

  SplitString(Buffer+15,Field,1);

  strcpy(Chain[CC]->Turn[TC]->Res1,Field[0]);

  SplitString(Buffer+26,Field,1);
  strcpy(Chain[CC]->Turn[TC]->Res2,Field[0]);

  strcpy(Tmp,Buffer);
  Tmp[24] = ' ';
  Tmp[35] = ' ';
  SplitString(Tmp+20,Field,1);
  strcpy(Chain[CC]->Turn[TC]->PDB_ResNumb1,Field[0]);
  SplitString(Tmp+31,Field,1);
  strcpy(Chain[CC]->Turn[TC]->PDB_ResNumb2,Field[0]);

  Chain[CC]->Turn[TC]->InsCode1 = Buffer[24];
  Chain[CC]->Turn[TC]->InsCode2 = Buffer[35];

  Chain[CC]->NTurn++;

  return(SUCCESS);
}


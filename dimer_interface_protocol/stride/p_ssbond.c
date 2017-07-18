#include "stride.h"

int Process_SSBOND(BUFFER Buffer, CHAIN **Chain, int *ChainNumber, COMMAND *Cmd)
{
  int CC, BC;
  char *Field[MAX_FIELD];
  BUFFER Tmp;

  if( Cmd->NActive && !ChInStr(Cmd->Active,SpaceToDash(Buffer[15])) )
     return(SUCCESS);

  CC = 0;

  if( *ChainNumber == 0 ) {
    InitChain(&Chain[CC]); 
    Chain[CC]->Id = Buffer[15];
    (*ChainNumber)++;
  }

  BC = Chain[CC]->NBond;
  Chain[CC]->SSbond[BC] =  (SSBOND *)ckalloc(sizeof(SSBOND));

  strcpy(Tmp,Buffer);
  Tmp[21] = ' ';
  Tmp[35] = ' ';
  SplitString(Tmp+17,Field,1);
  strcpy(Chain[CC]->SSbond[BC]->PDB_ResNumb1,Field[0]);
  SplitString(Tmp+31,Field,1);
  strcpy(Chain[CC]->SSbond[BC]->PDB_ResNumb2,Field[0]);

  Chain[CC]->SSbond[BC]->ChainId1 = Buffer[15];
  Chain[CC]->SSbond[BC]->ChainId2 = Buffer[29];

  Chain[CC]->SSbond[BC]->InsCode1 = Buffer[21];
  Chain[CC]->SSbond[BC]->InsCode2 = Buffer[35];

  Chain[CC]->SSbond[BC]->AsnSource = Pdb;

  Chain[CC]->NBond++;

  return(SUCCESS);
}


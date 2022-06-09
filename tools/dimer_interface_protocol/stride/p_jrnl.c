#include "stride.h"

int Process_JRNL(BUFFER Buffer, BOOLEAN *Published)
{
  char *Field[MAX_FIELD];

  SplitString(Buffer,Field,10);
  
  if( !strncmp(Field[1],"REF",3) && !strncmp(Field[2],"TO",2) && 
      !strncmp(Field[3],"BE",2)  && !strncmp(Field[4],"PUBLISHED",9) )
    *Published = NO;

  return(SUCCESS);
}


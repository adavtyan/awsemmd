#include "stride.h"

int Process_EXPDTA(BUFFER Buffer, enum METHOD *Method)
{
  char *Field[MAX_FIELD];
  int i, NFields;

  NFields = SplitString(Buffer,Field,10);
  
  for( i=0; i<NFields; i++ ) {
    if( !strcmp(Field[i],"MODEL") ) 
      *Method = Model;
    else 
    if( !strcmp(Field[i],"NMR") || !strcmp(Field[i],"/NMR$")  ) 
      *Method = NMR;
  }
  return(SUCCESS);
}


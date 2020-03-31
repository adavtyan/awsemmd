#include "stride.h"

int Process_REMARK(BUFFER Buffer, enum METHOD *Method, float *Resolution, BOOLEAN *DsspAssigned)
{
  char *Field[MAX_FIELD];
  int NFields;

  NFields = SplitString(Buffer,Field,10);
  
  if( NFields >= 5 &&  !strncmp(Field[2],"RESOLUTION",10) && 
      !strncmp(Field[4],"ANGSTROMS",9) && isdigit(*Field[3]) ) 
    *Resolution = atof(Field[3]);

  if( NFields >= 9 && !strcmp(Field[2],"THESE") && !strcmp(Field[3],"COORDINATES") &&
      !strcmp(Field[4],"WERE")  && !strcmp(Field[5],"GENERATED") &&
      !strcmp(Field[6],"FROM")  && !strcmp(Field[7],"SOLUTION") &&
      ( !strcmp(Field[8],"NMR") || !strcmp(Field[8],"/NMR$") ) ) *Method = NMR;

  if( strstr(Buffer,"SANDER ") || strstr(Buffer,"SANDER,") || strstr(Buffer,"SANDER:") || 
      strstr(Buffer,"SANDER;") || strstr(Buffer,"SANDER.") || strstr(Buffer,"SANDER(") || 
      strstr(Buffer,"SANDER)") || strstr(Buffer,"DSSP") )
    *DsspAssigned = YES;
  
  return(SUCCESS);
}


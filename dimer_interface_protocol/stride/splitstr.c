/* Split a char string into text fields */

#include "stride.h"

int SplitString(char *Buffer, char **Fields, int MaxField)
{
  int FieldCnt, SymbCnt, FieldFlag, BuffLen;
  static char LocalBuffer[BUFSZ];


  FieldCnt =0; FieldFlag = 0;
  BuffLen = (int)strlen(Buffer) - 1;

  strcpy(LocalBuffer,Buffer);

  for(SymbCnt=0; SymbCnt<BuffLen; SymbCnt++) {
    if( (isspace(LocalBuffer[SymbCnt])) && FieldFlag == 0 && SymbCnt != BuffLen-1 ) continue;
    if( (!isspace(LocalBuffer[SymbCnt])) && FieldFlag == 1 && SymbCnt == BuffLen-1 ) {
      LocalBuffer[SymbCnt+1] = '\0';
      return(FieldCnt);
    }
    else
    if( (isspace(LocalBuffer[SymbCnt])) && FieldFlag == 1 ) {
      LocalBuffer[SymbCnt] = '\0';
      FieldFlag = 0;
      if( FieldCnt == MaxField ) return(FieldCnt);
    }
    else
    if( (!isspace(LocalBuffer[SymbCnt])) && FieldFlag == 0 ) {
      FieldFlag = 1;
      Fields[FieldCnt] = LocalBuffer+SymbCnt;
      FieldCnt++;
    }
  }
  
     return(FieldCnt);
 }

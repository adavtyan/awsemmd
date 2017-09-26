#include <stdio.h>
#include <string.h>

void GetFileNameFromPath(char *Path, char *FileName)
{

  int i;
  static char DirDelim[5] = { ':','/','\\',']','\0'};

  for( i = (int)strlen(Path)-1; i>=0; i-- )
    if( strchr(DirDelim,Path[i]) ) break;

  strcpy(FileName,Path+i+1);
}


void StripPathFromLastExtention(char *Path, char *StrippedPath)
{
  int i;

  strcpy(StrippedPath,Path);

  for( i = (int)strlen(StrippedPath); i>=0; i-- )
    if( StrippedPath[i] == '.' ) {
      StrippedPath[i] = '\0';
      break;
    }
}

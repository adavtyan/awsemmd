#include <stdio.h>
#include <string.h>
#include <ctype.h>

#define BUFSZ 1024

char *tolostr(char *InputString)
{
  register int i;
  int Length;
  static char OutputString[BUFSZ];

  strcpy(OutputString,InputString);

  Length = (int)strlen(OutputString);

  for( i=0; i<Length; i++ )
    OutputString[i] = tolower(OutputString[i]);
  return(OutputString);
}

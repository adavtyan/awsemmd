#include "stride.h"

FILE *efopen(char *file, char *mode, char *progname)    /* fopen file, die if can't */
{
  FILE *fp;
  
  if( (fp=fopen(file,mode)) ) 
    return fp;
  else 
    die("%s: can't open file %s mode %s\n",progname,file,mode);
  return(FAILURE);
}


int Uniq(char **List, int ListLength)
{
    int i, j;

    for( i=1; i<ListLength-1; i++ ) {
      if( *List[i] != '-' ) continue;
      for( j=i+1; j<ListLength; j++ ) {
	if( *List[j] != '-' ) continue;
	if( !strcmp(List[i],List[j] ) ) return(0);
      }
    }
    
    return(1);
}

BOOLEAN Specified(char **List, int ListLength, char Option)
{
    int i;

    for( i=1; i<ListLength; i++ )
      if( *List[i] == '-' && *(List[i]+1) == Option ) 
	return(YES);

    return(NO);
}

int Parse(char **List, int ListLength, char *Option)
{
    int i;

    for( i=1; i<ListLength; i++ ) {
      if( *List[i] != '-' ) continue;
      if( !strcmp(List[i],Option) ) return(i);
    }
    
    return(0);
}

int CollectOptions(char **List, int ListLength, int Stream, int *Options)
{
    int OptCnt, i;

    OptCnt = 0;

    for( i=1; i<ListLength; i++ )
	if( *List[i] == '-' && !isdigit( *(List[i]+1) ) && atoi( List[i]+2 ) == Stream )
	    Options[OptCnt++] = i;

    return(OptCnt);

}




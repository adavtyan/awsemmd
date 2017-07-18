int Replace(char *String, char From, char To)
{

  int Replaced=0;

  if( From == '\0' )
    return(Replaced);

  for( ; *String != '\0'; String++ )
    if( *String == From ) {
      *String = To;
    Replaced++;
    }

  return(Replaced);
}

int Delete(char *String, char From)
{
  
  int Deleted = 0;
  char *c;
  
  if( From == '\0' )
    return(Deleted);
  
  for( ; *String != '\0'; String++ )
    if( *String == From ) {
      c = String;
      for( ;; c++ ) {
	*c = *(c+1);
	if( *c == '\0' ) 
	  break;
      }
      Deleted++;
      String--;
    }
  return(Deleted);
}
  

#include "stride.h"

int Process_COMPND(BUFFER Buffer, enum METHOD *Method)
{
  if( strstr(Buffer,"NMR") ) 
    *Method = NMR;

  if( strstr(Buffer,"MODEL") ) 
    if( *Method == XRay ) 
      *Method = Model;
  
  return(SUCCESS);
}



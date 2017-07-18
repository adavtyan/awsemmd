#include "stride.h"

int Process_MODEL(enum METHOD *Method)
{

  if( *Method == XRay ) 
    *Method = Model;

  return(SUCCESS);
}



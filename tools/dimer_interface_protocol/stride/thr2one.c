#include <string.h>
#include <stdio.h>

char ThreeToOne(char *Three)
{
  if( !strcmp(Three,"ALA") ) return('A');
  else if( !strcmp(Three,"ARG") ) return('R');
  else if( !strcmp(Three,"ASN") ) return('N');
  else if( !strcmp(Three,"ASP") ) return('D');
  else if( !strcmp(Three,"ASX") ) return('B');
  else if( !strcmp(Three,"CYS") ) return('C');
  else if( !strcmp(Three,"GLN") ) return('Q');
  else if( !strcmp(Three,"GLU") ) return('E');
  else if( !strcmp(Three,"GLX") ) return('Z');
  else if( !strcmp(Three,"GLY") ) return('G');
  else if( !strcmp(Three,"HIS") ) return('H');
  else if( !strcmp(Three,"ILE") ) return('I');
  else if( !strcmp(Three,"LEU") ) return('L');
  else if( !strcmp(Three,"LYS") ) return('K');
  else if( !strcmp(Three,"MET") ) return('M');
  else if( !strcmp(Three,"PRO") ) return('P');
  else if( !strcmp(Three,"PHE") ) return('F');
  else if( !strcmp(Three,"SER") ) return('S');
  else if( !strcmp(Three,"THR") ) return('T');
  else if( !strcmp(Three,"TRP") ) return('W');
  else if( !strcmp(Three,"TYR") ) return('Y');
  else if( !strcmp(Three,"VAL") ) return('V');
  else    return('X');
}
      

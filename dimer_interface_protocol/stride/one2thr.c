#include <string.h>
#include <stdio.h>

char *OneToThree(char One)
{

  if( One == 'A' ) return("ALA");
  else if( One == 'R' ) return("ARG");
  else if( One == 'N' ) return("ASN");
  else if( One == 'D' ) return("ASP");
  else if( One == 'B' ) return("ASX");
  else if( One == 'C' ) return("CYS");
  else if( One == 'Q' ) return("GLN");
  else if( One == 'E' ) return("GLU");
  else if( One == 'Z' ) return("GLX");
  else if( One == 'G' ) return("GLY");
  else if( One == 'H' ) return("HIS");
  else if( One == 'I' ) return("ILE");
  else if( One == 'L' ) return("LEU");
  else if( One == 'K' ) return("LYS");
  else if( One == 'M' ) return("MET");
  else if( One == 'P' ) return("PRO");
  else if( One == 'F' ) return("PHE");
  else if( One == 'S' ) return("SER");
  else if( One == 'T' ) return("THR");
  else if( One == 'W' ) return("TRP");
  else if( One == 'Y' ) return("TYR");
  else if( One == 'V' ) return("VAL");
  else if( One == 'X' ) return("UNK");
  else return("***");
}

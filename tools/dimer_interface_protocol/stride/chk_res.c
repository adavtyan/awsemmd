#include "stride.h"

int CheckRes(char *Res)
{
  register int ResTypeCnt;
  int ResTypeNumber = 45;
  static char *Rsd[MAX_ResType] = {
    "ACE", "ALA", "ARG", "ASN", "ASP", "ASX", "CYS", "GLN", "GLU", "GLX", "GLY", "HIS", "ILE", 
    "LEU", "LYS", "MET", "PRO", "PHE", "SER", "THR", "TRP", "TYR", "VAL", "FOR", "UNK", "HOH",
    /* Residues found empirically in the protein files */
/*  1gp1   1gp1   1hne   1tmn   2mcg   5hvp   6cha   1bbo   1ctg   act1   act1   aom1   rom7 */
    "SEC", "ALM", "MSU", "CLT", "PCA", "STA", "EXC", "ABU", "HYP", "HMS", "ASS", "OCT", "CYH",
/*  sod0  7adh   2psg   tim1   tim2   2pia */
    "MN", "INI", "PO3", "SUL", "WAT", "FMN"
    };

  for( ResTypeCnt=0; ResTypeCnt<ResTypeNumber; ResTypeCnt++ )
      if( !strcmp(Res,Rsd[ResTypeCnt]) )
	 return(SUCCESS);

  return(FAILURE);
}
      

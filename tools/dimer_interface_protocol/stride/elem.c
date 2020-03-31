#include "stride.h"

/*************************************************************************
**                                                                      **
** Calculate the number of individual secondary structure elements of   **
** type SecStrType and length not less than ElemLength that are:        **
**  - Present in Asn1 and Asn2 and absent  in Asn3 (YYN)                **
**  - Present in Asn2 and Asn3 and absent  in Asn1 (NYY)                **
**  - Absent  in Asn2 and Asn3 and present in Asn1 (YYN)                **
**  - Absent  in Asn1 and Asn2 and present in Asn3 (YYN)                **
**                                                                      **
*************************************************************************/

int FullElement(char *Asn1, char *Asn2, char *Asn3, int Length, char SecStrType, int ElemLength,
		 char EditChar, int *YYN, int *NYY, int *YNN, int *NNY)
{

  register int i, j, Count1, Count2, Count3;
  int Beg, ElLength;

  *YYN = 0;
  *NYY = 0;
  *YNN = 0;
  *NNY = 0;

  if( ElemLength >= Length )
    return(0);

  ElLength = ElemLength-1;
  Count1 = 0;
  Count2 = 0;
  Count3 = 0;

  Beg = -1;
  
  for( i=1; i<Length; i++ ) {
    if( ( i == 0 && 
	  ( Asn1[i] == SecStrType && Asn2[i] == SecStrType && Asn3[i] == SecStrType) ||
          ( Asn1[i] != SecStrType && Asn2[i] != SecStrType && Asn3[i] != SecStrType) )
       ||
        ( i  > 0 && 
	  ( Asn1[i] != Asn1[i-1] || Asn2[i] != Asn2[i-1] || Asn3[i] != Asn3[i-1] ) )
       ||
        i == Length-1 ) {
      
      if( Count1 >= ElLength && Count2 >= ElLength && Count3 <  ElLength ) 
	(*YYN)++;
      else
      if( Count1 <  ElLength && Count2 >= ElLength && Count3 >= ElLength ) 
	(*NYY)++;
      else
      if( Count1 >= ElLength && Count2 <  ElLength && Count3 <  ElLength ) 
	(*YNN)++;
      else
      if( Count1 <  ElLength && Count2 <  ElLength && Count3 >= ElLength ) 
	(*NNY)++;
      
/*       if( Count1 >= ElLength || Count2 >= ElLength || Count3 >= ElLength ) {
 * 	for( j=Beg-1; j<i; j++ ) {
 * 	  Asn1[j] = 'X';
 * 	  Asn2[j] = 'X';
 * 	  Asn3[j] = 'X';
 * 	}
 *       }
 * 
 */
      if( Count1 >= ElLength && ( Count2 < ElLength || Count3 < ElLength ) )
	for( j=Beg-1; j<i; j++ )
	  Asn1[j] = EditChar;

      if( Count2 >= ElLength && ( Count1 < ElLength || Count3 < ElLength ) )
	for( j=Beg-1; j<i; j++ )
	  Asn2[j] = EditChar;

      if( Count3 >= ElLength && ( Count1 < ElLength || Count2 < ElLength ) )
	for( j=Beg-1; j<i; j++ )
	  Asn3[j] = EditChar;

      Count1 = 0;
      Count2 = 0;
      Count3 = 0;
      Beg = -1;
      
    }
    else {
      if( Asn1[i] == SecStrType ) Count1++;
      if( Asn2[i] == SecStrType ) Count2++;
      if( Asn3[i] == SecStrType ) Count3++;
      if( Beg == -1 && (Count1 == 1 || Count2 == 1 || Count3 == 1) ) Beg = i;
    }
  }	

  CorrectAsn(Asn1,Length,SecStrType,EditChar,ElLength);
  CorrectAsn(Asn2,Length,SecStrType,EditChar,ElLength);
  CorrectAsn(Asn3,Length,SecStrType,EditChar,ElLength);

  return( (*YYN) * (*NYY) * (*YNN) * (*NNY) );
}


/*************************************************************************
**                                                                      **
** Calculate the number of individual secondary structure elements of   **
** type SecStrType in the known assignment Asn2 that are:               **
**  - Reproduced in Asn1 better than in Asn3 (Better)                   **
**  - Reproduced in Asn1 worse  than in Asn3 (Worse)                    **
**                                                                      **
*************************************************************************/

int CompareElements(char *Asn1, char *Asn2, char *Asn3, int Length, 
		   char SecStrType, int *Better, int *Worse)
{

  register int i, j, Count1, Count2;
  int TotalNumber = 0, Beg;

  *Better = 0;
  *Worse = 0;

  Beg = -1;
  
  for( i=0; i<Length; i++ ) {
    if( (Asn1[i] == SecStrType || Asn2[i] == SecStrType || Asn3[i] == SecStrType) &&
	(i == 0 || 
	 ( Asn1[i-1] != SecStrType && Asn2[i-1] != SecStrType && Asn3[i-1] != SecStrType) ) ) {
      TotalNumber++;
      Beg = i;
    }
    else
    if( Beg != -1 && ( i == Length-1 || 
	 ( Asn1[i] != SecStrType && Asn2[i] != SecStrType && Asn3[i] != SecStrType ) ) ) {
      Count1 = Count2 = 0;
      for( j=Beg; j<=i; j++ ) {
	if( (Asn1[j] == SecStrType || Asn2[j] == SecStrType) && Asn1[j] != Asn2[j] ) 
	  Count1++;
	if( (Asn3[j] == SecStrType || Asn2[j] == SecStrType) && Asn3[j] != Asn2[j] ) 
	  Count2++;
      }
      if( Count1 > Count2 ) 
	(*Worse)++;
      else
      if( Count2 > Count1 ) 
	(*Better)++;
      Beg = -1;
    }
  }
  return(TotalNumber);
}



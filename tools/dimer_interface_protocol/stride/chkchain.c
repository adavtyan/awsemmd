#include "stride.h"

int CheckChain(CHAIN *Chain, COMMAND *Cmd)
{

  int Res, HelAlp, HelPI, Hel310, Sheet, Turn, Bound[300][2];
  int i, j, AsnNumb=0, At, SuspCnt, Beg, End;
  float Content;

  if( Cmd->NProcessed && !ChInStr(Cmd->Processed,SpaceToDash(Chain->Id)) ) {
    Chain->Valid = NO;
    return(FAILURE);
  }

  if( Chain->NRes < 5 )
    return(NotValid(Chain,"less than 5 residues"));
    
  if( !Cmd->Stringent ) 
    return(SUCCESS);
  
  for( Res=0; Res<Chain->NRes; Res++ ){
    if( !CheckRes(Chain->Rsd[Res]->ResType) )
      Chain->NonStandRes++;
    for( At=0; At<Chain->Rsd[Res]->NAtom; At++) {
      if( !CheckAtom(Chain->Rsd[Res]->AtomType[At]) ) 
	Chain->NonStandAtom++;
      if(Chain->Rsd[Res]->Coord[At][0] < MIN_X || Chain->Rsd[Res]->Coord[At][0] > MAX_X ||
	 Chain->Rsd[Res]->Coord[At][1] < MIN_Y || Chain->Rsd[Res]->Coord[At][1] > MAX_Y ||
	 Chain->Rsd[Res]->Coord[At][2] < MIN_Z || Chain->Rsd[Res]->Coord[At][2] > MAX_Z ||
	 Chain->Rsd[Res]->Occupancy[At] < MIN_Occupancy || 
	 Chain->Rsd[Res]->Occupancy[At] > MAX_Occupancy ||
	 Chain->Rsd[Res]->TempFactor[At] < MIN_TempFactor || 
	 Chain->Rsd[Res]->TempFactor[At] > MAX_TempFactor )
	break;
    }      
    if( At < Chain->Rsd[Res]->NAtom )
      break;
  }
  
  if( Res < Chain->NRes )
    return(NotValid(Chain,"suspicious coordinates, occupancy or temperature factor"));
  
  if( 100.0*(float)Chain->NonStandRes/(float)Chain->NRes > MAXNONSTAND )
    return(NotValid(Chain,"too many non-standard residues"));
  
  if( Chain->NRes < Cmd->MinLength )
    return(NotValid(Chain,"Short chain"));
  
  if( Chain->NRes > Cmd->MaxLength )
    return(NotValid(Chain,"Long chain"));
  
  if( Chain->Method == XRay && 
     (Chain->Resolution < Cmd->MinResolution || Chain->Resolution > Cmd->MaxResolution ) )
    return(NotValid(Chain,"Resolution out of range"));
  
  if( (int)strlen(Cmd->Cond) != 0 )  {
    
    if( ChInStr(Cmd->Cond,'c') ) {
      for( Res=0; Res<Chain->NRes; Res++ )
	if( FindAtom(Chain,Res,"N",&At) || 
	    FindAtom(Chain,Res,"O",&At) || 
	    FindAtom(Chain,Res,"C",&At) ) 
	  break;
      
      if( Res == Chain->NRes )
	return(NotValid(Chain,"only CA"));
    }
    
    if( Chain->Method == NMR && !ChInStr(Cmd->Cond,'n') )
      return(NotValid(Chain,"NMR chain"));
    
    if( Chain->Method == XRay && !ChInStr(Cmd->Cond,'x') )
      return(NotValid(Chain,"XRay chain"));
    
    if( Chain->Method == Model && !ChInStr(Cmd->Cond,'m') )
      return(NotValid(Chain,"Model chain"));
    
    if( Chain->Published == NO && ChInStr(Cmd->Cond,'p') )
      return(NotValid(Chain,"Not published"));
    
    if( Chain->DsspAssigned == YES && ChInStr(Cmd->Cond,'d') )
      return(NotValid(Chain,"Assigned according to DSSP"));
    
    if( ChInStr(Cmd->Cond,'a') ) {
      
      if( Chain->Valid && Chain->NHelix == 0 && Chain->NSheet == -1 && Chain->NTurn == 0 )
	return(NotValid(Chain,"No assignment"));
      
      if( (Content = SecStrContent(Chain,&HelAlp,&HelPI,&Hel310,&Sheet,&Turn)) < 0.4 || 
	   Content > 0.9 )
	return(NotValid(Chain,"Suspicious content"));
      
      SuspCnt = 0;
      for( Res=1; Res<Chain->NRes-1; Res++ ) {
	if( ( Chain->Rsd[Res]->Prop->PdbAsn != 'H' && Chain->Rsd[Res]->Prop->PdbAsn != 'T' &&
	     Chain->Rsd[Res]->Prop->Phi > -150.0 && Chain->Rsd[Res]->Prop->Phi < 0.0 && 
	     Chain->Rsd[Res]->Prop->Psi > -100.0 && Chain->Rsd[Res]->Prop->Psi < 10.0) )
	  SuspCnt++;
	
      }
      
      if( (float)SuspCnt/(float)Chain->NRes > 0.4 )
	return(NotValid(Chain,"Suspicious assignment"));
      
      for( i=0; i<Chain->NHelix; i++ ) {
	if( !PdbN2SeqN(Chain,Chain->Helix[i]->PDB_ResNumb1,&Beg) ||
	   !PdbN2SeqN(Chain,Chain->Helix[i]->PDB_ResNumb2,&End) ||
	   /*	      !CheckRes(Chain->Helix[i]->PDB_ResNumb1) ||
		      !CheckRes(Chain->Helix[i]->PDB_ResNumb2) || */
	   Chain->Helix[i]->Class > 10 ||
	   Chain->Helix[i]->Class < 1  ||
	   End-Beg > 100 || End-Beg < 0 )
	  break;
	else 
	  if( Chain->Helix[i]->Class == 1 ) {
	    Bound[AsnNumb][0] = Beg;
	    Bound[AsnNumb][1] = End;
	    AsnNumb++;
	  }
      }
      if( i < Chain->NHelix )
	return(NotValid(Chain,"Erraneous helix assignment"));
      
      for( i=0; i<Chain->NSheet; i++ )
	for( j=0; j<Chain->Sheet[i]->NStrand; j++ ) {
	  if( !PdbN2SeqN(Chain,Chain->Sheet[i]->PDB_ResNumb1[j],&Beg) ||
	     !PdbN2SeqN(Chain,Chain->Sheet[i]->PDB_ResNumb2[j],&End) ||
	     /*	        !CheckRes(Chain->Sheet[i]->PDB_ResNumb1[j]) ||
			!CheckRes(Chain->Sheet[i]->PDB_ResNumb2[j]) || */
	     End-Beg > 100 || End-Beg < 0 )
	    break;
	  else
	    if( Chain->Sheet[i]->Sence[j] != 0 ) {
	      Bound[AsnNumb][0] = Beg;
	      Bound[AsnNumb][1] = End;
	      AsnNumb++;
	    }
	  if( j < Chain->Sheet[i]->NStrand )
	    break;
	}
      
      if( i < Chain->NSheet )
	return(NotValid(Chain,"Erraneous sheet assignment"));
      
      for( i=0; i<Chain->NTurn; i++ )
	if( !PdbN2SeqN(Chain,Chain->Turn[i]->PDB_ResNumb1,&Beg) ||
	   !PdbN2SeqN(Chain,Chain->Turn[i]->PDB_ResNumb2,&End) ||
	   End-Beg > 100 || End-Beg < 0 ) 
	  break;
      
      if( i < Chain->NTurn ) 
	NotValid(Chain,"Erraneous turn assignment");
      
      for( i=0; i<AsnNumb-1; i++ ) {
	for( j=i+1; j<AsnNumb; j++ ) {
	  if( Bound[i][0] == Bound[j][0] && Bound[i][1] == Bound[j][1] ) continue;
	  if( (Bound[j][0] > Bound[i][0] && Bound[j][0] < Bound[i][1]) ||
	     (Bound[j][1] > Bound[i][0] && Bound[j][1] < Bound[i][1]) ||
	     (Bound[i][0] > Bound[j][0] && Bound[i][0] < Bound[j][1]) ||
	     (Bound[i][1] > Bound[j][0] && Bound[i][1] < Bound[j][1]) )
	    break;
	}
	if( j < AsnNumb )
	  break;
      }
      
      if( i < AsnNumb-1 )
	return(NotValid(Chain,"Assignment overlap"));
    }
  }

  fprintf(stderr,"ACCEPTED %s %c %4d %7.3f\n",
	  Chain->File,Chain->Id,Chain->NRes,Chain->Resolution);
  return(SUCCESS);
}


int NotValid(CHAIN *Chain, char *Message)
{
  
  fprintf(stderr,"IGNORED %s %c ",Chain->File,SpaceToDash(Chain->Id));
  fprintf(stderr,"(%s)\n",Message);
  Chain->Valid = NO;
  return(FAILURE);
  
}




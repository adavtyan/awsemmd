#include "stride.h"

int ReadPDBFile(CHAIN **Chain, int *Cn, COMMAND *Cmd)
{

  int ChainCnt, InfoCnt, i;
  enum METHOD Method = XRay;
  BOOLEAN First_ATOM, Published=YES, DsspAssigned=NO;
  float Resolution = 0.0;
  FILE *pdb;
  BUFFER Buffer;
  char *Info[MAX_INFO], PdbIdent[5];
  RESIDUE *r;
  CHAIN *c;

  *Cn= 0;
  InfoCnt = 0;
  strcpy(PdbIdent,"~~~~");

  if( !(pdb = fopen(Cmd->InputFile,"r")) )
    return(FAILURE); 
     
  First_ATOM = YES;
  
  while( fgets(Buffer,BUFSZ,pdb) ) {
    
    if(!strncmp(Buffer,"HEADER",6)) {
      Info[InfoCnt] = (char *)ckalloc(BUFSZ*sizeof(char));
      strcpy(Info[InfoCnt],"HDR  ");
      strcat(Info[InfoCnt++],Buffer+10);
      strncpy(PdbIdent,Buffer+62,4);
      PdbIdent[4] = '\0';
    }
    else
    if(!strncmp(Buffer,"AUTHOR",6)) {
      Info[InfoCnt] = (char *)ckalloc(BUFSZ*sizeof(char));
      strcpy(Info[InfoCnt],"AUT  ");
      strcat(Info[InfoCnt++],Buffer+10);
    }
    else
    if(!strncmp(Buffer,"SOURCE",6)) {
      Info[InfoCnt] = (char *)ckalloc(BUFSZ*sizeof(char));
      strcpy(Info[InfoCnt],"SRC  ");
      strcat(Info[InfoCnt++],Buffer+10);
    }
    else
    if(!strncmp(Buffer,"COMPND",6)) {
      if( !Process_COMPND(Buffer,&Method) ) 
	return(FAILURE);
      else {
	Info[InfoCnt] = (char *)ckalloc(BUFSZ*sizeof(char));
	strcpy(Info[InfoCnt],"CMP  ");
	strcat(Info[InfoCnt++],Buffer+10);
      }
    }
    else if(!strncmp(Buffer,"JRNL",4) && !Process_JRNL(Buffer,&Published)) 
	return(FAILURE);
    else if(!strncmp(Buffer,"REMARK",6) && !Process_REMARK(Buffer,&Method,&Resolution,
							   &DsspAssigned)) 
      return(FAILURE);
    else if(!strncmp(Buffer,"EXPDTA",6) && !Process_EXPDTA(Buffer,&Method)) 
      return(FAILURE);
    else if(!strncmp(Buffer,"MODEL",5) && !Process_MODEL(&Method)) 
      return(FAILURE);
    else if(!strncmp(Buffer,"ENDMDL",6)) {
      Process_ENDMDL(Buffer,Chain,Cn);
      break;
    }
    else if(!strncmp(Buffer,"HELIX",5) && !Process_HELIX(Buffer,Chain,Cn,Cmd)) 
      return(FAILURE);
    else if(!strncmp(Buffer,"SHEET",5) && !Process_SHEET(Buffer,Chain,Cn,Cmd)) 
      return(FAILURE);
    else if(!strncmp(Buffer,"TURN",4) && !Process_TURN(Buffer,Chain,Cn,Cmd)) 
      return(FAILURE);
    else if(!strncmp(Buffer,"SSBOND",6) && !Process_SSBOND(Buffer,Chain,Cn,Cmd)) 
      return(FAILURE);
    else if(!strncmp(Buffer,"ATOM",4) && !Process_ATOM(Buffer,Chain,Cn,&First_ATOM,Cmd)) 
      return(FAILURE);
  }
  fclose(pdb);

  for( ChainCnt=0; ChainCnt< *Cn; ChainCnt++ ) {
    c = Chain[ChainCnt];
    if( c->NRes != 0  && !FindAtom(c,c->NRes,"CA",&i) )
      c->NRes--;
    strcpy(c->File,Cmd->InputFile);

    strcpy(c->PdbIdent,PdbIdent);
    if( c->NRes != 0 )  c->NRes++;
    if( c->NSheet != -1 ) c->NSheet++;
    c->Resolution = Resolution;
    c->Method = Method;
    c->Published = Published;
    c->DsspAssigned = DsspAssigned;
    c->NInfo = InfoCnt;
    for(i=0; i<InfoCnt; i++) {
      c->Info[i] = (char *)ckalloc(BUFSZ*sizeof(char));
      strcpy(c->Info[i],Info[i]);
      c->Info[i][71] = '\0';
    }
    for( i=0; i<c->NRes; i++ ) {
      r = c->Rsd[i];
      r->Inv =  (INVOLVED *)ckalloc(sizeof(INVOLVED));
      r->Prop = (PROPERTY *)ckalloc(sizeof(PROPERTY));
      r->Inv->NBondDnr = 0;
      r->Inv->NBondAcc = 0;
      r->Inv->InterchainHBonds = NO;
      r->Prop->Asn     = 'C';
      r->Prop->PdbAsn  = 'C';
      r->Prop->DsspAsn = 'C';
      r->Prop->Solv    = 0.0;
      r->Prop->Phi     = 360.0;
      r->Prop->Psi     = 360.0;
    }
  }
  
  for(i=0; i<InfoCnt; i++)
    free(Info[i]);
  
  return(SUCCESS);
}


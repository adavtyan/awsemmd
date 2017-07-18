#include "stride.h"

/* #include <console.h> */  /* For Macintosh only, see readme.mac */

int main(int argc, char **argv)
{

  CHAIN **Chain;
  HBOND **HBond;
  COMMAND *Cmd;
  int Cn, NChain=0, NHBond=0, ValidChain=0;
  float **PhiPsiMapHelix, **PhiPsiMapSheet;
  register int i;

  /* argc = ccommand(&argv); */ /* For Macintosh only, see readme.mac */



  Chain = (CHAIN  **)ckalloc(MAX_CHAIN*sizeof(CHAIN *));
  HBond = (HBOND  **)ckalloc(MAXHYDRBOND*sizeof(HBOND *));
  Cmd   = (COMMAND *)ckalloc(sizeof(COMMAND));

  ProcessStrideOptions(argv,argc,Cmd);

  if( !ReadPDBFile(Chain,&NChain,Cmd) || !NChain )
    die("Error reading PDB file %s\n",Cmd->InputFile);
  
  for( Cn=0; Cn<NChain; Cn++ )
    ValidChain += CheckChain(Chain[Cn],Cmd);

/*   if( Cmd->Stringent )
 *     exit(0);
 */
  if( !ValidChain ) 
    die("No valid chain in %s\n",Chain[0]->File);
  
  if( Cmd->BrookhavenAsn )
    GetPdbAsn(Chain,NChain);

  if( Cmd->DsspAsn )
    GetDsspAsn(Chain,NChain,Cmd);

  BackboneAngles(Chain,NChain);

  if( Cmd->OutSeq )
    OutSeq(Chain,NChain,Cmd);

  if( Cmd->ContactOrder )
    ContactOrder(Chain,NChain,Cmd);

  if( Cmd->ContactMap )
    ContactMap(Chain,NChain,Cmd);

  if( !strlen(Cmd->MapFileHelix) )
    PhiPsiMapHelix = DefaultHelixMap(Cmd);
  else
    ReadPhiPsiMap(Cmd->MapFileHelix,&PhiPsiMapHelix,Cmd);

  if( !strlen(Cmd->MapFileSheet) )
    PhiPsiMapSheet = DefaultSheetMap(Cmd);
  else
    ReadPhiPsiMap(Cmd->MapFileSheet,&PhiPsiMapSheet,Cmd);
  
  for( Cn=0; Cn<NChain; Cn++ )
    PlaceHydrogens(Chain[Cn]);
  
  if( (NHBond = FindHydrogenBonds(Chain,Cn,HBond,Cmd)) == 0 ) 
    die("No hydrogen bonds found in %s\n",Cmd->InputFile);
  
  NoDoubleHBond(HBond,NHBond);
  
  DiscrPhiPsi(Chain,NChain,Cmd);
  
  if(Cmd->ExposedArea)
    Area(Chain,NChain,Cmd);

  for( Cn=0; Cn<NChain; Cn++ ) {

    if( Chain[Cn]->Valid ) {
    
      Helix(Chain,Cn,HBond,Cmd,PhiPsiMapHelix);
      
      for( i=0; i<NChain; i++ ) 
	if( Chain[i]->Valid )
	  Sheet(Chain,Cn,i,HBond,Cmd,PhiPsiMapSheet);    
      
      BetaTurn(Chain,Cn);
      GammaTurn(Chain,Cn,HBond);
      
    }
  }
    
  Report(Chain,NChain,HBond,Cmd);

  if( Cmd->MolScript )
    MolScript(Chain,NChain,Cmd);

  for( i=0; i<Cn; i++ ) free(Chain[i]);
  for( i=0; i<NHBond; i++ ) free(HBond[i]);
  free(Cmd);

  return(0);
}

void ProcessStrideOptions(char **List, int ListLength, COMMAND *Cmd)
{
  
  int i, InpFile = 0;
  char OPTION;
  BOOLEAN Password = NO;

  if( Uniq(List,ListLength) == 0 ) { 
    fprintf(stderr,"All options must be unique\n"); 
    PrintStrideHelp(Cmd); 
  }
 
  DefaultCmd(Cmd);

  Password = Specified(List,ListLength,'$');

  for( i=1; i<ListLength; i++ ) {
    if( *List[i] == '-' ) {
      
      OPTION = toupper((*(List[i]+1)));

      /*********************** Process public options ************************/
      if(OPTION == 'M') {
	strcpy(Cmd->MolScriptFile,List[i]+2); 
	Cmd->MolScript = YES;
      }
      else if( OPTION == 'O' ) Cmd->ReportSummaryOnly = YES;
      else if( OPTION == 'H' ) Cmd->ReportBonds = YES;
      else if( OPTION == 'R' ) strcpy(Cmd->Active,List[i]+2);
      else if( OPTION == 'C' ) strcpy(Cmd->Processed,List[i]+2);
      else if( OPTION == 'F' ) strcpy(Cmd->OutFile,List[i]+2);
      else if( OPTION == 'Q' ) {
	strcpy(Cmd->SeqFile,List[i]+2);
	Cmd->OutSeq = YES;
      }
      /*********************** Process private options ************************/
      else if( OPTION == 'I' && Password ) Cmd->Info = YES;
      else if( OPTION == 'Z' && Password ) Cmd->Measure = YES;
      else if( OPTION == 'K' && Password ) Cmd->ContactOrder = YES;
      else if( OPTION == 'W' && Password ) Cmd->ContactMap = YES;
      else if( OPTION == 'X' && Password ) strcpy(Cmd->FirstResidue,List[i]+2);
      else if( OPTION == 'Y' && Password ) strcpy(Cmd->LastResidue,List[i]+2);
      else if( OPTION == 'D' && Password ) {
	strcpy(Cmd->DsspFile,List[i]+2);
	Cmd->DsspAsn = YES;
      }
      else if( OPTION == 'B' && Password ) Cmd->BrookhavenAsn = YES;
      else
      if(OPTION == 'P'&& Password ) {
	if( toupper(*(List[i]+2)) == 'H' ) 
	  strcpy(Cmd->MapFileHelix,List[i]+3); 
	else
	if( toupper(*(List[i]+2)) == 'E' ) 
	  strcpy(Cmd->MapFileSheet,List[i]+3); 
	else
	  PrintStrideHelp(Cmd);
      }
      else
      if( OPTION == 'T' && Password ) {
	if( toupper(*(List[i]+2)) == 'H' ) {
	  if( toupper(*(List[i]+3)) == 'A' )
	    Cmd->Treshold_H1 = atof(List[i]+4);
	  else
	  if( toupper(*(List[i]+3)) == 'C' ) 
	    Cmd->Treshold_H3 = atof(List[i]+4);
	  else
	  if( toupper(*(List[i]+3)) == 'D' ) 
	    Cmd->Treshold_H4 = atof(List[i]+4);
	  else
	    PrintStrideHelp(Cmd);
	}
	else
	if( toupper(*(List[i]+2)) == 'E' ) {
	  if( toupper(*(List[i]+3)) == 'A' ) 
	    Cmd->Treshold_E1 = atof(List[i]+4);
	  else
	  if( toupper(*(List[i]+3)) == 'B' ) 
	    Cmd->Treshold_E2 = atof(List[i]+4);
	  else
	  if( toupper(*(List[i]+3)) == 'C' ) 
	    Cmd->Treshold_E3 = atof(List[i]+4);
	  else
	  if( toupper(*(List[i]+3)) == 'D' ) 
	    Cmd->Treshold_E4 = atof(List[i]+4);
	  else
	    PrintStrideHelp(Cmd);
	}
      }
      else if( OPTION == 'S' && Password ) {
	Cmd->Stringent = YES;
	if( (int)strlen(List[i]+2) > MAXCONDITIONS-1 )
	  PrintStrideHelp(Cmd);	
	strcpy(Cmd->Cond,List[i]+2);
      }
      else if( OPTION == 'L' && Password ) {
	if( toupper(*(List[i]+2)) == 'L' ) Cmd->MaxLength = atoi(List[i]+3);
	if( toupper(*(List[i]+2)) == 'G' ) Cmd->MinLength = atoi(List[i]+3);
      }
      else if( OPTION == 'U' && Password ) {
	if( toupper(*(List[i]+2)) == 'L' ) Cmd->MaxResolution = atof(List[i]+3);
	if( toupper(*(List[i]+2)) == 'G' ) Cmd->MinResolution = atof(List[i]+3);
      }
      else
	if( OPTION != '$' )
	  PrintStrideHelp(Cmd);
    }
    else {
      strcpy(Cmd->InputFile,List[i]);
      InpFile++;
    }
  }

  if( InpFile > 1 ) {
    fprintf(stderr,"\nOnly one input file is allowed\n");
    PrintStrideHelp(Cmd); 
  }
  else
  if( !InpFile ) {
    fprintf(stderr,"\nYou must specify input file \n");
    PrintStrideHelp(Cmd); 
  }

  Cmd->NActive = (int)strlen(Cmd->Active);
  Cmd->NProcessed = (int)strlen(Cmd->Processed);

  if( Cmd->Measure ) {
    Cmd->BrookhavenAsn = YES;
    Cmd->DsspAsn = YES;
  }
}


void PrintStrideHelp(COMMAND *Cmd)
{

  fprintf(stderr,"\nAction: secondary structure assignment\n");
  fprintf(stderr,"Usage: stride [Options] InputFile [ > file ]\n");
  fprintf(stderr,"Options:  \n");
  fprintf(stderr,"  -fFile      Output file\n");
  fprintf(stderr,"  -mFile      MolScript file\n");
  fprintf(stderr,"  -o          Report secondary structure summary Only\n");
  fprintf(stderr,"  -h          Report Hydrogen bonds\n");
  fprintf(stderr,"  -rId1Id2..  Read only chains Id1, Id2 ...\n");
  fprintf(stderr,"  -cId1Id2..  Process only Chains Id1, Id2 ...\n");
  fprintf(stderr,"  -q[File]    Generate SeQuence file in FASTA format and die\n");
  fprintf(stderr,"\nOptions are position  and case insensitive\n");

  /*************** Private options - not for general use ****************/
  if( Cmd->Info ) {
    fprintf(stderr,"  -phFile        PhiPsiMap file for alpha helix\n");
    fprintf(stderr,"  -peFile        PhiPsiMap file for beta-sheet\n");
    fprintf(stderr,"  -pgFile        PhiPsiMap file for 3-1 helix\n");
    fprintf(stderr,"  -tha           treshold (Ener*Conf*Conf) for H\n");
    fprintf(stderr,"  -thc           treshold (Conf) for H\n");
    fprintf(stderr,"  -thd           treshold (Conf) for H\n");
    fprintf(stderr,"  -tea           treshold (Ener*Conf*Conf) for E\n");
    fprintf(stderr,"  -teb           treshold (Conf) for E\n");
    fprintf(stderr,"  -dFile         Dssp file\n");
    fprintf(stderr,"  -b             Include Brookhaven assignment\n");
    fprintf(stderr,"  -i             Print detailed information\n");
    fprintf(stderr,"  -llNumber      Sequence Length less then Number (with -s)\n");
    fprintf(stderr,"  -lgNumber      Sequence Length greater then Number (with -s)\n");
    fprintf(stderr,"  -ulNumber      ResolUtion less (better) then Number (with -s)\n");
    fprintf(stderr,"  -ugNumber      ResolUtion greater (worse) then Number (with -s)\n");
    fprintf(stderr,"  -s[xcnmapd]    Stringent: discard wrong or suspisious PDB files \n");
    fprintf(stderr,"       -x Process X-Ray structures\n");
    fprintf(stderr,"       -n Process NMR structures\n");
    fprintf(stderr,"       -m Process Model structures\n");
    fprintf(stderr,"       -c Exlude structures with only Calpha atoms\n");
    fprintf(stderr,"       -d Exculde structures with automatic secondary structure assignment\n");
    fprintf(stderr,"       -a Exclude structures without secondary structure assignment\n");
    fprintf(stderr,"       -p Exlude structures that are not published or in press\n");
    fprintf(stderr,"  -z             Measure differences with PDB and DSSP\n");
    
  }

  exit(0);

}

void DefaultCmd(COMMAND *Cmd)
{

  Cmd->SideChainHBond    = NO;
  Cmd->MainChainHBond    = YES;
  Cmd->MainChainPolarInt = YES;
  Cmd->Published         = NO;
  Cmd->DsspAssigned      = NO;
  Cmd->UseResolution     = NO;
  Cmd->Info              = NO;
  Cmd->Truncate          = YES;
  Cmd->ExposedArea       = YES;
  Cmd->ReportSummaryOnly = NO;
  Cmd->ReportBonds       = NO;
  Cmd->BrookhavenAsn     = NO;
  Cmd->DsspAsn           = NO;
  Cmd->MolScript         = NO;
  Cmd->OutSeq            = NO;
  Cmd->Stringent         = NO;
  Cmd->Measure           = NO;

  Cmd->EnergyType        = 'G';

  Cmd->DistCutOff        =  6.0;
  Cmd->PhiPsiStep        =  0.0;

  Cmd->C1_H              = -1.0;
  Cmd->C2_H              =  1.0;
  Cmd->C1_E              = -0.2;
  Cmd->C2_E              =  0.2;

  Cmd->Treshold_H1       = -230.0;
  Cmd->Treshold_H3       =  0.12;
  Cmd->Treshold_H4       =  0.06;
  Cmd->Treshold_E1       = -240.0;
  Cmd->Treshold_E2       = -310.0;

  Cmd->MinResolution     =  0.1;
  Cmd->MaxResolution     =  100.0;

  Cmd->MinLength         = 0;
  Cmd->MaxLength         = MAX_RES;
  
  Cmd->NPixel            = 0;
  Cmd->NActive           = 0;
  Cmd->NProcessed        = 0;

  strcpy(Cmd->FirstResidue,"");
  strcpy(Cmd->LastResidue,"");

  strcpy(Cmd->MapFileHelix,""); 
  strcpy(Cmd->MapFileSheet,""); 
  strcpy(Cmd->OutFile,"");
  strcpy(Cmd->Active,"");
  strcpy(Cmd->Processed,"");
  strcpy(Cmd->Cond,"");


}

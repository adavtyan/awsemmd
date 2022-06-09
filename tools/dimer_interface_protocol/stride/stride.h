#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <stdarg.h>

#define Pi                        3.1415927
#define Eps                       0.000001
#define Minimum(x,y)              ((x)<(y) ? x : y)
#define Maximum(x,y)              ((x)<(y) ? y : x)
#define Sign(x)                   ((x)<0 ? -1 : 1)
#define IN(x, target, range)      ( (x >= (target - range)) && (x <= (target + range)) )
#define RAD(x)                    (x)*Pi/180.0
#define DEG(x)                    (x)*180.0/Pi
#define RADDEG                    57.2958
#define BREAKDIST                 2.5
#define SSDIST                    3.0

#define SUCCESS                   1
#define FAILURE                   0
#define YES                       1
#define NO                        0
#define ERR                      -1

#define BUFSZ                     1024
#define MAX_FIELD                 50

#define MAX_AtomType              200
#define MAX_ResType               50
#define MAXNONSTAND               4.0
#define MAX_CHAIN                 100
#define MAX_RES                   20000
#define MAX_HETRES                20000
#define MAX_HET                   200
#define MAX_HELIX                 500
#define MAX_SHEET                 500
#define MAX_STRAND_IN_SHEET       20
#define MAX_TURN                  300
#define MAX_BOND                  100
#define MAX_ASSIGN                500
#define MAX_INFO                  1000
#define MAX_AT_IN_RES             75
#define MAX_AT_IN_HETERORES       200
#define MAXRESDNR                 6
#define MAXRESACC                 6 

#define RES_FIELD                 6
#define AT_FIELD                  5

#define MAX_X                180.000
#define MAX_Y                180.000
#define MAX_Z                180.000
#define MIN_X               -100.000
#define MIN_Y               -100.000
#define MIN_Z               -100.000
#define MAX_Occupancy           1.00
#define MIN_Occupancy           0.00
#define MAX_TempFactor       1000.00
#define MIN_TempFactor          0.00

#define OUTPUTWIDTH               80
#define MAXCONDITIONS             20

#define MAXHYDRBOND               50000
#define MAXDONOR                  MAX_RES
#define MAXACCEPTOR               MAX_RES

#define	MINPHIPSI     	         -180.0
#define MAXPHIPSI      	          180.0
#define DEFNUMPIXEL               18

#define DIST_N_H                  1.0
#define RmGRID                    3.0
#define EmGRID                   -2.8
#define CGRID                    -3.0*EmGRID*pow(RmGRID,8.0)
#define DGRID                    -4.0*EmGRID*pow(RmGRID,6.0)
#define K1GRID                    0.9/pow(cos(RAD(110.0)),6.0)
#define K2GRID                    pow(cos(RAD(110.0)),2.0)

#define MINACCANG_SP2   90.0
#define MAXACCANG_SP2   180.0
#define MINACCANG_SP3   60.0
#define MAXACCANG_SP3   180.0
#define MINDONANG_SP2   90.0
#define MAXDONANG_SP2   180.0
#define MINDONANG_SP3   90.0
#define MAXDONANG_SP3   180.0
#define ACCDONANG       60.0
#define DONACCANG       90.0

#define DSSPPATH                  "/data/dssp/"

#define NAcd            20

enum ASNSOURCE {Stride, Pdb, Dssp};
enum METHOD {XRay, NMR, Model};
enum HYBRID {Nsp2, Nsp3, Osp2, Osp3, Ssp3};
enum GROUP  {Peptide, Trp, Asn, Gln, Arg, His, Lys, Ser, Thr, Tyr, Asp, Glu, Met, Cys};
enum HBONDTYPE {MM, MS, SM, SS };

typedef char BUFFER[BUFSZ+1];

typedef int BOOLEAN;

#define MAX_FILE 500

typedef struct {
                 char *FileName[MAX_FILE];
		 char *ChainId[MAX_FILE];
	       } CHAINLIST;

typedef struct {
                 float Phi, Psi, Omega;
		 int PhiZn, PsiZn;
		 float Solv, DsspSolv;
		 char Asn, PdbAsn, DsspAsn;
	       } PROPERTY;

typedef struct {
                 int HBondDnr[MAXRESDNR];
                 int HBondAcc[MAXRESACC];
		 int NBondDnr, NBondAcc;
		 BOOLEAN InterchainHBonds;
	       } INVOLVED;

typedef struct {
                 int NAtom;
		 char PDB_ResNumb[RES_FIELD];
		 char ResType[RES_FIELD];
		 char AtomType[MAX_AT_IN_RES][AT_FIELD];
		 float Coord[MAX_AT_IN_RES][3];
		 float Occupancy[MAX_AT_IN_RES];
		 float TempFactor[MAX_AT_IN_RES];
		 PROPERTY *Prop;
		 INVOLVED *Inv;
	       } RESIDUE;

typedef struct {
                 int NAtom;
		 char PDB_ResNumb[RES_FIELD];
		 char ResType[RES_FIELD];
		 char AtomType[MAX_AT_IN_HETERORES][AT_FIELD];
		 char Mendeleev[MAX_AT_IN_HETERORES];
		 float Coord[MAX_AT_IN_HETERORES][3];
		 float Occupancy[MAX_AT_IN_HETERORES];
		 float TempFactor[MAX_AT_IN_HETERORES];
	       } HETERORESIDUE;

typedef struct {
                 char HetId[4];
		 char PDB_ResNumb[RES_FIELD];
		 char InsCode;
		 int AtomNumb;
	       } HET;

typedef struct {
                 char Res1[RES_FIELD];
		 char Res2[RES_FIELD];
		 char PDB_ResNumb1[RES_FIELD], PDB_ResNumb2[RES_FIELD];
		 char InsCode1, InsCode2;
		 int Class;
	       } HELIX;

typedef struct {
                 int  NStrand;
		 char SheetId[RES_FIELD];
		 char ResType1[MAX_STRAND_IN_SHEET][RES_FIELD];
		 char ResType2[MAX_STRAND_IN_SHEET][RES_FIELD];
		 char PDB_ResNumb1[MAX_STRAND_IN_SHEET][RES_FIELD];
		 char PDB_ResNumb2[MAX_STRAND_IN_SHEET][RES_FIELD];
		 char InsCode1[MAX_STRAND_IN_SHEET];
		 char InsCode2[MAX_STRAND_IN_SHEET];
		 int  Sence[MAX_STRAND_IN_SHEET];
		 int  RegYN[MAX_STRAND_IN_SHEET];
		 char AtomNameReg1[MAX_STRAND_IN_SHEET][AT_FIELD];
		 char AtomNameReg2[MAX_STRAND_IN_SHEET][AT_FIELD];
		 char ResTypeReg1[MAX_STRAND_IN_SHEET][RES_FIELD];
		 char ResTypeReg2[MAX_STRAND_IN_SHEET][RES_FIELD];
		 char PDB_ResNumbReg1[MAX_STRAND_IN_SHEET][RES_FIELD];
		 char PDB_ResNumbReg2[MAX_STRAND_IN_SHEET][RES_FIELD];
		 char InsCodeReg1[MAX_STRAND_IN_SHEET];
		 char InsCodeReg2[MAX_STRAND_IN_SHEET];
	       } SHEET;

typedef struct {
                 char Res1[RES_FIELD];
		 char Res2[RES_FIELD];
		 char PDB_ResNumb1[RES_FIELD], PDB_ResNumb2[RES_FIELD];
		 char InsCode1, InsCode2;
		 char TurnType;
	       } TURN;

typedef struct {
                 char PDB_ResNumb1[RES_FIELD], PDB_ResNumb2[RES_FIELD];
                 char InsCode1, InsCode2;
		 char ChainId1, ChainId2;
		 enum ASNSOURCE AsnSource;
	       } SSBOND;

typedef struct {
                 int NRes, NHetRes, NonStandRes, Ter;
		 int NHet, NAtom, NonStandAtom, NHelix, NSheet;
		 int NTurn, NAssignedTurn, NBond, NHydrBond, NHydrBondInterchain, NHydrBondTotal, NInfo;
		 char Id, *File;
		 float Resolution;
		 enum METHOD Method;
		 BOOLEAN Valid, Published, DsspAssigned;

		 RESIDUE **Rsd;
		 HETERORESIDUE **HetRsd;
		 HET     **Het;
		 HELIX   **Helix;
		 SHEET   **Sheet;
		 TURN    **Turn;
		 TURN    **AssignedTurn;
		 SSBOND  **SSbond;
		 char    **Info;
		 char    PdbIdent[5];

	       } CHAIN;


typedef struct {
                 CHAIN *Chain;
		 int D_Res, DD_Res, DDI_Res;
 		 int D_At, DD_At, DDI_At, H;
		 enum HYBRID Hybrid;
		 enum GROUP Group;
		 float HB_Radius;
	       } DONOR;

                 
typedef struct {
                 CHAIN *Chain;
                 int A_Res, AA_Res, AA2_Res;
		 int A_At, AA_At, AA2_At;
		 enum HYBRID Hybrid;
		 enum GROUP Group;
		 float HB_Radius;
	       } ACCEPTOR;

typedef struct { 
                 BUFFER InputFile, OutFile, SeqFile;
		 BUFFER MapFileHelix, MapFileSheet;
		 BUFFER MolScriptFile, DsspFile;
		 char EnergyType, Active[MAX_CHAIN+1]; 
		 char Processed[MAX_CHAIN+1], Cond[MAXCONDITIONS];
                 char FirstResidue[RES_FIELD], LastResidue[RES_FIELD];

                 int NPixel, NActive, NProcessed;
		 int MinLength, MaxLength;

		 float PhiPsiStep, DistCutOff;
		 float Treshold_H1, Treshold_H2, Treshold_H3, Treshold_H4;
		 float Treshold_E1, Treshold_E2, Treshold_E3, Treshold_E4;
		 float MinResolution, MaxResolution;
		 float C1_H, C2_H, C1_E, C2_E;

		 BOOLEAN SideChainHBond, MainChainHBond, MainChainPolarInt;
		 BOOLEAN Published, DsspAssigned, UseResolution, Info, Truncate;
		 BOOLEAN ExposedArea, ReportSummaryOnly, ReportBonds, BrookhavenAsn, DsspAsn;
		 BOOLEAN MolScript, OutSeq, Stringent, Measure, ContactOrder, ContactMap;

		 /* Not used by STRIDE */
		 BOOLEAN Shrink;
		 BUFFER CarteFile, MapFile, MathFile;
		 char AsnSource, Mode, SecStrType;
		 int FilterOrder;
		 int NStepA, NStepB, NStepC, NStepD;
		 float Treshold_From_A, Treshold_To_A, StepA;
		 float Treshold_From_B, Treshold_To_B, StepB;
		 float Treshold_From_C, Treshold_To_C, StepC;
		 float Treshold_From_D, Treshold_To_D, StepD;

	       } COMMAND;

typedef struct {
                 DONOR *Dnr;
		 ACCEPTOR *Acc;
		 BOOLEAN ExistPolarInter, ExistHydrBondRose, ExistHydrBondBaker;
		 float Energy, Er, Et, Ep, ti, to, p;
		 float AccDonDist, OHDist, AngNHO, AngCOH;
		 float AccAng, DonAng, AccDonAng, DonAccAng;
	       } HBOND;

typedef struct PAT {
                 HBOND *Hb1, *Hb2;
		 struct PAT *Nei1, *Nei2;
		 BOOLEAN ExistPattern;
		 BUFFER Type;
	       } PATTERN;

typedef struct {
		 BUFFER File;
		 char Id;
                 int NRes;
		 char **ResType;
		 char *SecondStr;
		 char **PDB_ResNumb;
		 float *Accessibility;
	       } DSSP;

typedef struct {
                 int TP, TN, FP, FN;
		 float Corr, Perc;
	       } QUALITY;

#include "protot.h"
#include "nsc.h"


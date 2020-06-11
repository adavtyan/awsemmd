
/*********************************************************************
**********************************************************************
****        <<<<< 3SPN.2 Configuration Generator    >>>>>         ****
****                                                              ****  
****        Dan Hinckley    `                                     ****
****        Instutite for Molecular Engineering                   ****
****        University of Chicago                                 ****
****        Chicago, IL 60637                                     ****
****                                                              ****
**********************************************************************
*********************************************************************/

#include "decl_ilib.h"
#include "decl_ivar.h"
#include "decl_ifun.h"

void wrte_lammps(char *dnme)
{
    long i, ncyc;
    char fnme[_TW_], angleStyle[_TW_];
    FILE *fptr;

    // Conversion factor for converting from kJ to kCal
    double kJ2kCal = 0.239005736; 
    double k2bond, k3bond, k4bond,kbend, ktors, bhfx, bhfy, bhfz;
    long stea, steb, stec, sted, typa, typb, typc, typd, match, j;

    // Getting the correct units for bonded interactions
    k2bond = 0.600000 * kJ2kCal;
    k3bond = 0.000000 * kJ2kCal;
    k4bond = k2bond * 100.0;
    kbend = 200.000000 * kJ2kCal;
    ktors = 6.000000 * kJ2kCal;

    sprintf(fnme, "%s/conf_lammps.in", dnme);
    fptr = fopen(fnme, "w");

    fprintf(fptr, "LAMMPS Description for 3SPN.2\n\n");
    ncyc = site.totl;

    // Writing down the number of atoms and bonded interactions
    fprintf(fptr, "\t%ld atoms\n",ncyc);
    fprintf(fptr, "\t%ld bonds\n",cbon);
    fprintf(fptr, "\t%ld angles\n",cben);
    fprintf(fptr, "\t%ld dihedrals\n\n",ctor);
    
    long ntype_atom, ntype_bond, ntype_angle, ntype_dihedral;
        ntype_atom = 14;
    if (dna_type != 1) {
        ntype_bond = 6;
        ntype_angle = 26;
        ntype_dihedral = 2;
    } else {
        ntype_bond = 1;
        ntype_angle = 17;
        ntype_dihedral = 1;
    }
    // Specifying the number of types
    fprintf(fptr, "\t%ld atom types\n",ntype_atom);
    fprintf(fptr, "\t%ld bond types\n",ntype_bond);
    fprintf(fptr, "\t%ld angle types\n",ntype_angle);
    fprintf(fptr, "\t%ld dihedral types\n\n",ntype_dihedral);

    // Specifying the box size
    bhfx = boxs.lenx / 2.0;
    bhfy = boxs.leny / 2.0;
    bhfz = boxs.lenz / 2.0;
    fprintf(fptr, "\t%lf\t%lf xlo xhi\n", -bhfx ,bhfx);
    fprintf(fptr, "\t%lf\t%lf ylo yhi\n", -bhfy ,bhfy);
    fprintf(fptr, "\t%lf\t%lf zlo zhi\n\n", -bhfz ,bhfz);

    // Specifying the masses (gm/mol)
    fprintf(fptr,"Masses\n\n");
    fprintf(fptr,"\t%d\t%lf\n", 1,94.9696); // P
    fprintf(fptr,"\t%d\t%lf\n", 2,83.1104); // S
    fprintf(fptr,"\t%d\t%lf\n", 3,134.1220); // A
    fprintf(fptr,"\t%d\t%lf\n", 4,125.1078); // T
    fprintf(fptr,"\t%d\t%lf\n", 5,150.1214); // G
    fprintf(fptr,"\t%d\t%lf\n", 6,110.0964); // C
    fprintf(fptr,"\t%d\t%lf\n", 7,134.1220); // 5A
    fprintf(fptr,"\t%d\t%lf\n", 8,125.1078); // 5T
    fprintf(fptr,"\t%d\t%lf\n", 9,150.1214); // 5G
    fprintf(fptr,"\t%d\t%lf\n", 10,110.0964); // 5C
    fprintf(fptr,"\t%d\t%lf\n", 11,134.1220); // 3A
    fprintf(fptr,"\t%d\t%lf\n", 12,125.1078); // 3T
    fprintf(fptr,"\t%d\t%lf\n", 13,150.1214); // 3G
    fprintf(fptr,"\t%d\t%lf\n\n", 14,110.0964); // 3C
    // fprintf(fptr,"\t%d\t%lf\n", 15,22.9898); // Na+
    // fprintf(fptr,"\t%d\t%lf\n", 16,24.305); // Mg2+
    // fprintf(fptr,"\t%d\t%lf\n", 17,35.453); // Cl-
    // fprintf(fptr,"\t%d\t%lf\n\n", 18,49.4); // N+

    // I use this same ordering of the indices for the rest of my forcefield. 
    // not consistent with the in-house 3SPN.2. I won't worry about it. 

    // epsilon sigma R_solv lambda Gsolve
    double eps_r = 1.000000 * kJ2kCal; 
    // Remember, the units are in kcal/mol and angtroms or radians

    // Specifying the pair coefficients
    fprintf(fptr,"Bond Coeffs\n\n");
    double blen[ntype_bond];
    if ((dna_type == 0) || (dna_type == 1)) {
        blen[0] = 3.899;   // PS
        blen[1] = 3.559;   // SP
        blen[2] = 4.670;   // SA
        blen[3] = 4.189;   // ST
        blen[4] = 4.829;   // SG
        blen[5] = 4.112;   // SC
    } else  {
        blen[0] = 4.157;   // PS
        blen[1] = 3.780;   // SP
        blen[2] = 4.697;   // SA
        blen[3] = 4.220;   // ST
        blen[4] = 4.852;   // SG
        blen[5] = 4.066;   // SC
    }
    if (dna_type != 1) {
        for (i = 0; i < ntype_bond; i++) {
            fprintf(fptr,"\t%ld\t%lf\t%lf\t%lf\t%lf\n",i+1,blen[i],k2bond,k3bond,k4bond);
        }
    } else {
        fprintf(fptr,"\t%ld\tlist\n",(long)1);
    }
    fprintf(fptr,"\n");

    // Specifying the angle coefficients
    fprintf(fptr,"Angle Coeffs\n\n");

    // Harmonic <K> <theta0>

    double bAngle[10];
    if  ((dna_type == 0) || (dna_type == 1)) {
        bAngle[0] = 94.49;   // SPS
        bAngle[1] = 120.15;  // PSP
        bAngle[2] = 112.07;  // ASP
        bAngle[3] = 116.68;  // TSP
        bAngle[4] = 110.12;  // GSP
        bAngle[5] = 114.34;  // CSP
        bAngle[6] = 103.53;  // PSA
        bAngle[7] = 92.06;   // PST
        bAngle[8] = 107.40;  // PSG
        bAngle[9] = 96.96;   // PSC
        if (dna_type == 0) {
            sprintf(angleStyle,"%s","harmonic");
        } else {
            sprintf(angleStyle,"%s","list");
        }
    } else  {
        bAngle[0] = 92.77;  // SPS
        bAngle[1] = 91.24;  // PSP
        bAngle[2] = 104.86; // ASP
        bAngle[3] = 110.58; // TSP
        bAngle[4] = 103.86; // GSP
        bAngle[5] = 106.94; // CSP
        bAngle[6] = 103.71; // PSA
        bAngle[7] = 93.27;  // PST
        bAngle[8] = 107.49; // PSG
        bAngle[9] = 97.58;  // PSC
        sprintf(angleStyle,"%s","harmonic");
    }
    if (dna_type != 1) {
        for (i = 0; i < 10; i++) {
            fprintf(fptr,"\t%ld\t%s\t%lf\t%lf\n",i+1,angleStyle,kbend,bAngle[i]); 
        }
    }
    else
    {
        fprintf(fptr,"\t%ld\tlist\n",(long)1); // List angle style takes no arguments
    }

    // A = 0; T = 1; G = 2; C = 3
    double bstk_sigm[16], bstk_thta[16], bstk_eps[16], bstk_alpha, bstk_range;
    bstk_alpha = 3.0;
    bstk_range = 6.0;

    if (dna_type == 0) { 
        // B- DNA without curvature
        bstk_sigm[0] = 3.716;
        bstk_sigm[1] = 3.675;
        bstk_sigm[2] = 3.827;
        bstk_sigm[3] = 3.744;
        bstk_sigm[4] = 4.238;
        bstk_sigm[5] = 3.984;
        bstk_sigm[6] = 4.416;
        bstk_sigm[7] = 4.141;
        bstk_sigm[8] = 3.576;
        bstk_sigm[9] = 3.598;
        bstk_sigm[10] = 3.664;
        bstk_sigm[11] = 3.635;
        bstk_sigm[12] = 4.038;
        bstk_sigm[13] = 3.798;
        bstk_sigm[14] = 4.208;
        bstk_sigm[15] = 3.935;

        bstk_thta[0] = 101.15;      // AA
        bstk_thta[1] = 85.94;       // AT
        bstk_thta[2] = 105.26;      // AG
        bstk_thta[3] = 89.00;       // AC
        bstk_thta[4] = 101.59;      // TA
        bstk_thta[5] = 89.5;        // TT
        bstk_thta[6] = 104.31;      // TG
        bstk_thta[7] = 91.28;       // TC
        bstk_thta[8] = 100.89;      // GA
        bstk_thta[9] = 84.83;       // GT
        bstk_thta[10] = 105.48;     // GG
        bstk_thta[11] = 88.28;      // GC
        bstk_thta[12] = 106.49;     // CA
        bstk_thta[13] =  93.31;     // CT
        bstk_thta[14] = 109.54;     // CG
        bstk_thta[15] =  95.46;     // CC

    } else if (dna_type == 1) {
        // With curvature
        bstk_sigm[0] = 3.58;
        bstk_sigm[1] = 3.56;
        bstk_sigm[2] = 3.85;
        bstk_sigm[3] = 3.45;
        bstk_sigm[4] = 4.15;
        bstk_sigm[5] = 3.93;
        bstk_sigm[6] = 4.32;
        bstk_sigm[7] = 3.87;
        bstk_sigm[8] = 3.51;
        bstk_sigm[9] = 3.47;
        bstk_sigm[10] = 3.67;
        bstk_sigm[11] = 3.42;
        bstk_sigm[12] = 4.15;
        bstk_sigm[13] = 3.99;
        bstk_sigm[14] = 4.34;
        bstk_sigm[15] = 3.84;

        bstk_thta[0] = 100.13;
        bstk_thta[1] = 90.48;
        bstk_thta[2] = 104.39;
        bstk_thta[3] = 93.23;
        bstk_thta[4] = 102.59;
        bstk_thta[5] = 93.32;
        bstk_thta[6] = 103.70;
        bstk_thta[7] = 94.55;
        bstk_thta[8] = 95.45;
        bstk_thta[9] = 87.63;
        bstk_thta[10] = 106.36;
        bstk_thta[11] = 83.12;
        bstk_thta[12] = 102.69;
        bstk_thta[13] = 96.05;
        bstk_thta[14] = 100.46;
        bstk_thta[15] = 100.68;
    } else { 
        // A-DNA
        bstk_sigm[0] = 4.022;
        bstk_sigm[1] = 3.344;
        bstk_sigm[2] = 4.261;
        bstk_sigm[3] = 3.737;
        bstk_sigm[4] = 4.794;
        bstk_sigm[5] = 4.031;
        bstk_sigm[6] = 5.064;
        bstk_sigm[7] = 4.445;
        bstk_sigm[8] = 3.855;
        bstk_sigm[9] = 3.217;
        bstk_sigm[10] = 4.077;
        bstk_sigm[11] = 3.592;
        bstk_sigm[12] = 4.499;
        bstk_sigm[13] = 3.708;
        bstk_sigm[14] = 4.772;
        bstk_sigm[15] = 4.116;

        bstk_thta[0] = 108.32;
        bstk_thta[1] = 96.74;
        bstk_thta[2] = 111.32;
        bstk_thta[3] = 97.36;
        bstk_thta[4] = 103.33;
        bstk_thta[5] = 94.85;
        bstk_thta[6] = 105.36;
        bstk_thta[7] = 94.51;
        bstk_thta[8] = 108.25;
        bstk_thta[9] = 95.59;
        bstk_thta[10] = 111.66;
        bstk_thta[11] = 96.71;
        bstk_thta[12] = 111.39;
        bstk_thta[13] = 102.73;
        bstk_thta[14] = 113.47;
        bstk_thta[15] = 102.14;
    }


    if ((dna_type == 0) || (dna_type == 2)) {

        bstk_eps[0] = 14.39;  // AA
        bstk_eps[1] = 14.34;  // AT
        bstk_eps[2] = 13.25;  // AG
        bstk_eps[3] = 14.51;  // AC
        bstk_eps[4] = 10.37;  // TA
        bstk_eps[5] = 13.36;  // TT
        bstk_eps[6] = 10.34;  // TG
        bstk_eps[7] = 12.89;  // TC
        bstk_eps[8] = 14.81;  // GA
        bstk_eps[9] = 15.57;  // GT
        bstk_eps[10] = 14.93; // GG
        bstk_eps[11] = 15.39; // GC
        bstk_eps[12] = 11.42; // CA
        bstk_eps[13] = 12.79; // CT
        bstk_eps[14] = 10.52; // CG
        bstk_eps[15] = 13.24; // CC
    }
    else if (dna_type == 1 ){
       bstk_eps[0] =  13.820 ;// AA
       bstk_eps[1] =  15.050 ;// AT
       bstk_eps[2] =  13.320 ;// AG
       bstk_eps[3] =  15.820 ;// AC
       bstk_eps[4] =   9.150 ;// TA
       bstk_eps[5] =  12.440 ;// TT
       bstk_eps[6] =   9.580 ;// TG
       bstk_eps[7] =  13.110 ;// TC
       bstk_eps[8] =  13.760 ;// GA
       bstk_eps[9] =  14.590 ;// GT
       bstk_eps[10] = 14.770 ;// GG
       bstk_eps[11] = 15.170 ;// GC
       bstk_eps[12] =  9.250 ;// CA
       bstk_eps[13] = 12.420 ;// CT
       bstk_eps[14] =  8.830 ;// CG
       bstk_eps[15] = 14.010 ;// CC
    }



    for (i = 0; i < 16; i++)
    {   
        if (dna_type != 1) fprintf(fptr,"\t%ld\tstacking/3spn2\t%lf\t%lf\t%lf\t%lf\t%lf\n",i+11,bstk_eps[i] * kJ2kCal,bstk_sigm[i], bstk_thta[i], bstk_alpha, bstk_range);
        else fprintf(fptr,"\t%ld\tstacking/3spn2\t%lf\t%lf\t%lf\t%lf\t%lf\n",i+2,bstk_eps[i] * kJ2kCal,bstk_sigm[i], bstk_thta[i], bstk_alpha, bstk_range);
    }
    fprintf(fptr,"\n");

    // Specifying the dihedral coefficients
    // <k> <phi> <sigm>
    fprintf(fptr,"Dihedral Coeffs\n\n");
    if (dna_type == 0) {
        fprintf(fptr,"\t%d\t%lf\t%lf\t%lf\n",1,ktors, -359.17,0.3000); // SPSP
        fprintf(fptr,"\t%d\t%lf\t%lf\t%lf\n",2,ktors, -334.79,0.3000); // PSPS
    } else if (dna_type == 1) {
        fprintf(fptr,"\t%ld\tlist\n",(long)1); // List style does not take any parameters
    } else {
        fprintf(fptr,"\t%d\t%lf\t%lf\t%lf\n",1,ktors, -9.58 ,0.3000); // SPSP
        fprintf(fptr,"\t%d\t%lf\t%lf\t%lf\n",2,ktors,-328.40,0.3000); // PSPS
    }
    fprintf(fptr,"\n");

    // First, we identify the 5' and 3' sites...
    int prime_shift[ncyc];
    long ndna = site.dna_totl;
    for (i = 0; i < ndna; i++)
    {
        prime_shift[i] = 0;
        if (i == 1) prime_shift[i] = 4;
        if (i == ndna - 1) prime_shift[i] = 8;
        if (comp)
        {
            if (i == ndna/2 - 1) prime_shift[i] = 8;
            if (i == ndna/2 + 1) prime_shift[i] = 4;
        }
    }

    // Specifying the Atoms
    fprintf(fptr,"Atoms\n\n");
    long moln = 1;
    for (i = 0; i < ndna; i++)
    {
        fprintf(fptr, "\t%ld",i+1); // Atom number
        if (comp)
        {
            if (i == ndna/ 2)
            {
                moln++;
            }
        }
        fprintf(fptr, "\t%ld",moln); // Molecule number
        fprintf(fptr, "\t%ld",atom[i].stid2 + prime_shift[i]); // Atom type
        fprintf(fptr, "\t%lf",atom[i].chrg); // Atom charge
        fprintf(fptr, "\t%lf",atom[i].xval); // Atom x coord
        fprintf(fptr, "\t%lf",atom[i].yval); // Atom y coord
        fprintf(fptr, "\t%lf\n",atom[i].zval); // Atom z coord
    }

    // Add ions
    if (ions_flag)
    {
        for (i = ndna; i < site.totl; i++)
        {
            moln++;
            fprintf(fptr, "\t%ld",i + 1); // Atom number
            fprintf(fptr, "\t%ld",moln); // Molecule number
            fprintf(fptr, "\t%ld",atom[i].stid2); // Atom type
            fprintf(fptr, "\t%lf",atom[i].chrg); // Atom charge
            fprintf(fptr, "\t%lf",atom[i].xval); // Atom x coord
            fprintf(fptr, "\t%lf",atom[i].yval); // Atom y coord
            fprintf(fptr, "\t%lf\n",atom[i].zval); // Atom z coord
        }
    }
    fprintf(fptr,"\n");

    // Specifying the bonds

    // First I classify each of the bonds
    long bt[6][2];
    bt[0][0] = 1; bt[0][1] = 2; // PS
    bt[1][0] = 2; bt[1][1] = 1; // SP
    bt[2][0] = 2; bt[2][1] = 3; // SA
    bt[3][0] = 2; bt[3][1] = 4; // ST
    bt[4][0] = 2; bt[4][1] = 5; // SG
    bt[5][0] = 2; bt[5][1] = 6; // SC

    long b[cbon];
	for(i=0; i < cbon; i++)
	{	
        stea = bond[i].aste;
		steb = bond[i].bste;
		typa = atom[stea].stid2;
		typb = atom[steb].stid2;
		match = 0;
        if (dna_type != 1) {
            for(j = 0; j < ntype_bond; j++)
            {	
                if( (typa == bt[j][0]) && (typb == bt[j][1]) )
                {	
                    match = 1;
                    b[i] = j+1;
                    break;
                }
            }
            if(!match)
            {	
                printf("WARNING: classify(): could not classify bond %ld\n", i+1);
            }
        } else b[i] = 1;
	}

    if (cbon > 0) fprintf(fptr,"Bonds\n\n");
    for (i = 0; i < cbon; i++)
    {
        fprintf(fptr,"\t%ld\t%ld\t%ld\t%ld\n",i+1,b[i],bond[i].aste + 1, bond[i].bste + 1);
    }
    fprintf(fptr,"\n");

    // Specifying the angles
    // First I classify each of the bonds

    long a[cben];
    if (dna_type != 1) {
        long at[26][3];
        at[0][0] = 2; at[0][1] = 1; at[0][2] = 2;// SPS
        at[1][0] = 1; at[1][1] = 2; at[1][2] = 1;// PSP
        at[2][0] = 3; at[2][1] = 2; at[2][2] = 1;// ASP
        at[3][0] = 4; at[3][1] = 2; at[3][2] = 1;// TSP
        at[4][0] = 5; at[4][1] = 2; at[4][2] = 1;// GSP
        at[5][0] = 6; at[5][1] = 2; at[5][2] = 1;// CSP
        at[6][0] = 1; at[6][1] = 2; at[6][2] = 3;// PSA
        at[7][0] = 1; at[7][1] = 2; at[7][2] = 4;// PST
        at[8][0] = 1; at[8][1] = 2; at[8][2] = 5;// PSG
        at[9][0] = 1; at[9][1] = 2; at[9][2] = 6;// PSC

        for (i = 0; i < 16; i++)
        {
            at[10 + i][0] = 2;
            at[10 + i][1] = i / 4 + 3;
            at[10 + i][2] = i % 4 + 3;
        }

        for(i=0; i < cben; i++)
        {	
            stea = bend[i].aste;
            steb = bend[i].bste;
            stec = bend[i].cste;
            typa = atom[stea].stid2;
            typb = atom[steb].stid2;
            typc = atom[stec].stid2;
            match = 0;
            for(j = 0; j < ntype_angle; j++)
            {	
                if( (typa == at[j][0]) && (typb == at[j][1]) && (typc == at[j][2]) )
                {	
                    match = 1;
                    a[i] = j+1;
                    break;
                }
            }
            if(!match)
            {	
                printf("WARNING: classify(): could not classify bend %ld\n", i+1);
            }
        }
    } else {
        long at[16][3];
        for (i = 0; i < 16; i++)
        {
            at[i][0] = 2;
            at[i][1] = i / 4 + 3;
            at[i][2] = i % 4 + 3;
        }

        for(i=0; i < cben; i++)
        {	
            stea = bend[i].aste;
            steb = bend[i].bste;
            stec = bend[i].cste;
            typa = atom[stea].stid2;
            typb = atom[steb].stid2;
            typc = atom[stec].stid2;
            match = 0;
            for(j = 0; j < ntype_angle-1; j++)
            {	
                if( (typa == at[j][0]) && (typb == at[j][1]) && (typc == at[j][2]) )
                {	
                    match = 1;
                    a[i] = j+2;
                    break;
                }
            }
            if(!match) a[i] = 1;
        }
    }
    
    if (cben > 0) fprintf(fptr,"Angles\n\n");
    for (i = 0; i < cben; i++)
    {
        fprintf(fptr,"\t%ld\t%ld\t%ld\t%ld\t%ld\n",(long)i+1,a[i],bend[i].aste + 1, bend[i].bste + 1, bend[i].cste + 1); 
    }
    fprintf(fptr,"\n");

    // Specifying the dihedarls

    // First I classify each of the bonds
    long dt[2][4];
    dt[0][0] = 2; dt[0][1] = 1; dt[0][2] = 2; dt[0][3] = 1; // SPSP
    dt[1][0] = 1; dt[1][1] = 2; dt[1][2] = 1; dt[1][3] = 2; // PSPS

    long d[ctor];
	for(i=0; i < ctor; i++)
	{	
        stea = tors[i].aste;
		steb = tors[i].bste;
		stec = tors[i].cste;
		sted = tors[i].dste;
		typa = atom[stea].stid2;
		typb = atom[steb].stid2;
		typc = atom[stec].stid2;
		typd = atom[sted].stid2;
		match = 0;
        if (dna_type != 1) {
            for(j = 0; j < ntype_dihedral; j++)
            {	
                if( (typa == dt[j][0]) && (typb == dt[j][1]) && (typc == dt[j][2]) && (typd == dt[j][3]) )
                {	
                    match = 1;
                    d[i] = j+1;
                    break;
                }
            }
            if(!match)
            {	
                printf("WARNING: classify(): could not classify dihedral %ld\n", i+1);
            }
        } else d[i] = 1;
	}

    if (ctor > 0) fprintf(fptr,"Dihedrals\n\n");
    for (i = 0; i < ctor; i++)
    {
        fprintf(fptr,"\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\n",i+1,d[i],tors[i].aste + 1, tors[i].bste + 1, tors[i].cste + 1, tors[i].dste +1);
    }
    fclose(fptr);

    // B-DNA
    double rA,rT,rG,rC;
    char dnaString[_TW_];
    if ((dna_type == 0) || (dna_type == 1)) {
        rA = 2.7;
        rT = 3.55;
        rG = 2.45;
        rC = 3.2;
        sprintf(dnaString,"B-DNA");
    } else {
        rA = 2.23;
        rT = 2.75;
        rG = 2.1;
        rC = 2.85;
        sprintf(dnaString,"A-DNA");
    }

    // Write to screen the values for containing 3SPN.2 parameters
    fprintf(stdout,"################################################################\n");
    fprintf(stdout,"Pair coefficients for 3SPN.2 representation of %s\n",dnaString);
    fprintf(stdout,"%s\t%d\t%d\t%s\t%lf\t%lf\n","pair_coeff",1,1,"3spn2",eps_r, 2.0 * 2.25);
    fprintf(stdout,"%s\t%d\t%d\t%s\t%lf\t%lf\n","pair_coeff",2,2,"3spn2",eps_r, 2.0 * 3.10); 
    fprintf(stdout,"%s\t%d\t%d\t%s\t%lf\t%lf\n","pair_coeff",3,3,"3spn2",eps_r, 2.0 * rA); 
    fprintf(stdout,"%s\t%d\t%d\t%s\t%lf\t%lf\n","pair_coeff",4,4,"3spn2",eps_r, 2.0 * rT); 
    fprintf(stdout,"%s\t%d\t%d\t%s\t%lf\t%lf\n","pair_coeff",5,5,"3spn2",eps_r, 2.0 * rG); 
    fprintf(stdout,"%s\t%d\t%d\t%s\t%lf\t%lf\n","pair_coeff",6,6,"3spn2",eps_r, 2.0 * rC); 
    fprintf(stdout,"%s\t%d\t%d\t%s\t%lf\t%lf\n","pair_coeff",7,7,"3spn2",eps_r, 2.0 * rA); 
    fprintf(stdout,"%s\t%d\t%d\t%s\t%lf\t%lf\n","pair_coeff",8,8,"3spn2",eps_r, 2.0 * rT); 
    fprintf(stdout,"%s\t%d\t%d\t%s\t%lf\t%lf\n","pair_coeff",9,9,"3spn2",eps_r, 2.0 * rG); 
    fprintf(stdout,"%s\t%d\t%d\t%s\t%lf\t%lf\n","pair_coeff",10,10,"3spn2",eps_r, 2.0 * rC);
    fprintf(stdout,"%s\t%d\t%d\t%s\t%lf\t%lf\n","pair_coeff",11,11,"3spn2",eps_r, 2.0 * rA);
    fprintf(stdout,"%s\t%d\t%d\t%s\t%lf\t%lf\n","pair_coeff",12,12,"3spn2",eps_r, 2.0 * rT);
    fprintf(stdout,"%s\t%d\t%d\t%s\t%lf\t%lf\n","pair_coeff",13,13,"3spn2",eps_r, 2.0 * rG);
    fprintf(stdout,"%s\t%d\t%d\t%s\t%lf\t%lf\n","pair_coeff",14,14,"3spn2",eps_r, 2.0 * rC);
    fprintf(stdout,"################################################################\n");
}

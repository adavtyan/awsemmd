#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


// Modify this function to remap atom types if necessary
inline int remap_arom_type(int type)
{
/*	if (type<=400) return 1;
	if (type==401) return 2;
	if (type==402) return 3;*/

	return type;
}

class Atom
{
public:
	Atom(): id(0), type(0) {
		x[0] = x[1] = x[2] = 0.0;
	}
	~Atom() {}
	void set(int id, int ty, double x, double y, double z) {
		this->id = id;
		this->type = ty;
		this->x[0] = x;
		this->x[1] = y;
		this->x[2] = z;
	}

	int id, type;
	double x[3];
};

class Histograms
{
public:
	Histograms(int nty, double rmax, double dr) {
		int i, j;

		this->ntypes = nty;
		this->rmax = rmax;
		this->dr = dr;

		nbins = int(round(rmax/dr));
		nhists = ntypes*(ntypes+1)/2;

		hist_map = new int*[ntypes];
		for (i=0;i<ntypes;i++)
			hist_map[i] = new int[ntypes];

		type_counts = new double[ntypes];
		hist_sums = new double*[nhists];
		for (i=0;i<nhists;i++) 
			hist_sums[i] = new double[nbins];

		for (i=0;i<ntypes;i++) {
			type_counts[i] = 0.0;
			for (j=0;j<=i;j++) {
				hist_map[i][j] = i*(i+1)/2+j;
				hist_map[j][i] = i*(i+1)/2+j;
			}
		}

		for (i=0;i<nhists;i++) {
			for(j=0;j<nbins;j++) {
				hist_sums[i][j] = 0.0;
			}
		}
	}

	~Histograms() {
		int i;

		for (i=0;i<ntypes;i++) delete [] hist_map[i];
		for (i=0;i<nhists;i++) delete [] hist_sums[i];

		delete [] hist_map;
		delete [] type_counts;
		delete [] hist_sums;
	}

	inline void hist_add(int ity, int jty, int bin) {
		hist_sums[hist_map[ity][jty]][bin] += 1.0;
	}

	inline void type_add(int ity) {
		type_counts[ity] += 1.0;
	}

	inline void ihist2types(int ih, int &ity, int &jty) {
		ity = int(0.5*(sqrt(1+8*(double)ih)-1));
		jty = ih - ity*(ity+1)/2;
	}

	void calc_rdfs(double V_avg, int n_frame) {	
		double dV;
		int ity, jty, i, j;
		for (i=0;i<nhists;i++) {
			ihist2types(i, ity, jty);
			for (j=0;j<nbins;j++) {
				dV = 4*M_PI*dr*dr*dr*(3*j*j+3*j+1)/3; // 4pi/3 * [r(i+1)^3 - r(i)^3] = 4pi/3 * dr^3 * [(i+1)^3 - i^3]
				hist_sums[i][j] *= n_frame*V_avg/(type_counts[ity]*type_counts[jty]*dV);
				if (ity==jty) hist_sums[i][j] *=2;
			}
		}
	}

	int nhists, **hist_map;
	int nbins, ntypes;
	double rmax, dr;
	double *type_counts, **hist_sums;
};

double distance(const double *x1, const double *x2, const double *prd);
void calc_rdfs(Atom *atoms, int n_atoms, double *prd, Histograms *hists);

int main(int argc, const char* argv[])
{
	if (argc<5) {
		printf("\n> %s Input_file Output_file dr rmax [start] [end] [freq]\n\n", argv[0]);
		return -1;
	}

	// Read input/output file names
	const char *infile=argv[1];
	const char *outfile=argv[2];

	double dr = 0.5;
	double rmax = 20.0;

	dr = atof(argv[3]);
	rmax = atof(argv[4]);

	// Read the rest of the arguments
	int start=0, end=-1, freq=1;
	if (argc>5) start = atoi(argv[5]);
	if (argc>6) end = atoi(argv[6]);
	if (argc>7) freq = atoi(argv[7]);
	
	// Open input file
	FILE *input = fopen(infile, "r");
	if (!input) {
		printf("\nERROR: Input file not found!\n\n");
		return -1;
	}

	// Read input file and analyze frames one by one
	int i, atom_id, atom_type, nstr;
	double x, y, z;
	Atom *atoms=NULL;
	Histograms *hists=NULL;
	int n_frame = 0;
	int i_frame = 0;
	int step = 0;
	int n_atoms = 0;
	int n_max_type = 0;
	double V_avg=0.0;
	double A[3][2], prd[3]; // Box boundaries and sizes
	char line[128], buff[128], item[32], *str[10];
	while ( fgets(line, sizeof(line), input) != NULL ) {
		if (strncmp(line, "ITEM:", 5)==0) {
			strcpy(item, line+6);

			if (strncmp(item, "TIMESTEP", 8)==0) {
				fgets(buff, sizeof(buff), input);
				step = atoi(buff);
			} else if (strncmp(item, "NUMBER OF ATOMS", 15)==0) {
				fgets(buff, sizeof(buff), input);
				n_atoms = atoi(buff);
			} else if (strncmp(item, "BOX BOUNDS", 10)==0) {
				for (i=0;i<3;i++) {
					fgets(buff, sizeof(buff), input);
					A[i][0] = atof(strtok(buff," \t\n"));
					A[i][1] = atof(strtok(NULL," \t\n"));
				}
			} else if (strncmp(item, "ATOMS", 5)==0) {
//				if (step==1000000) {
//					printf(step, start, end, i_frame, freq);
//				}
				if (step>=start && (end==-1 || step<=end) && i_frame%freq==0) {
					if (n_frame%100==0)
						printf("Processing frame # %d\n", step);

					if (atoms==NULL) atoms = new Atom[n_atoms];

					// Read coordinates
					for (i=0;i<n_atoms;++i) {
						fgets(buff, sizeof(buff), input);
						
						// Split the line string and read the first 5 values 
						nstr = 0;
						str[nstr] = strtok(buff," \t\n");
						while ( str[nstr]!=NULL ) {
							nstr++;
							if (nstr>=5) break;
							str[nstr] = strtok(NULL," \t\n");
						}
						if (nstr!=5) printf("ERROR: Error reading ATOMS section!\n");
						
						atom_id = atoi(str[0]);
						atom_type = atoi(str[1]);
						x = atof(str[2]);
						y = atof(str[3]);
						z = atof(str[4]);

						atom_type = remap_arom_type(atom_type);

						if (hists==NULL && atom_type>n_max_type) n_max_type = atom_type;

						x = (A[0][1] - A[0][0])*x + A[0][0];
						y = (A[1][1] - A[1][0])*y + A[1][0];
						z = (A[2][1] - A[2][0])*z + A[2][0];
						
						atoms[i].set(atom_id, atom_type, x, y, z);
					}

					if (hists==NULL) hists = new Histograms(n_max_type, rmax, dr);

					prd[0] = (A[0][1] - A[0][0]);
					prd[1] = (A[1][1] - A[1][0]);
					prd[2] = (A[2][1] - A[2][0]);

					V_avg += prd[0]*prd[1]*prd[2];
					calc_rdfs(atoms, n_atoms, prd, hists);
					n_frame++;
				}
				if (step>=start) i_frame++;
			}
		}
	}

	// Final calculation step for RDFs
	if (n_frame>0) {
		V_avg /= n_frame;
		hists->calc_rdfs(V_avg, n_frame);
	}

	printf("Number of frames %d\n", n_frame);
	printf("Number of atoms %d\n", n_atoms);
	printf("Number of atom types %d\n", n_max_type);
	printf("Avarage volume %f\n", V_avg);
	
	fclose(input);

	// Open output file
	FILE *output = fopen(outfile, "w");
	if (!output) {
		printf("\nERROR: Output file not found!\n");
		return -1;
	}

	// Write header header
	int ity, jty;
	fprintf(output, "# r");
	for (i=0;i<hists->nhists;i++) {
		hists->ihist2types(i,ity,jty);
		fprintf(output, "\t%d-%d", ity, jty);
	}
	fprintf(output, "\n");
	
	// Write histograms
	int j;
	double r;
	for (i=0;i<hists->nbins;i++) {
		r = i*dr+dr/2;
		fprintf(output, "%8.4f", r);
		for (j=0;j<hists->nhists;j++) {
			fprintf(output, "\t%8.6f", hists->hist_sums[j][i]);
		}
		fprintf(output, "\n");
	}

	fclose(output);

	return 0;
}

void calc_rdfs(Atom *atoms, int n_atoms, double *prd, Histograms *hists)
{
	int i, j, ity, jty, bin;
	double d;
	for (i=0;i<n_atoms;++i) {
		ity = atoms[i].type-1;
		hists->type_add(ity);
		for (j=i+1;j<n_atoms;++j) {
			jty = atoms[j].type-1;
			d = distance(atoms[i].x, atoms[j].x, prd);
			if (d<=hists->rmax) {
				bin = int(d/hists->dr);
				hists->hist_add(ity, jty, bin);
			}
		}
	}
}

double distance(const double *x1, const double *x2, const double *prd)
{
	double half_prd[3], v[3];

	half_prd[0] = prd[0]/2;	
	half_prd[1] = prd[1]/2;	
	half_prd[2] = prd[2]/2;	

	v[0] = x2[0] - x1[0];
	v[1] = x2[1] - x1[1];
	v[2] = x2[2] - x1[2];

	if (v[0]<=-half_prd[0]) v[0] += prd[0];
	else if (v[0]>half_prd[0]) v[0] -= prd[0];
	if (v[1]<=-half_prd[1]) v[1] += prd[1];
	else if (v[1]>half_prd[1]) v[1] -= prd[1];
	if (v[2]<=-half_prd[2]) v[2] += prd[2];
	else if (v[2]>half_prd[2]) v[2] -= prd[2];

	return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

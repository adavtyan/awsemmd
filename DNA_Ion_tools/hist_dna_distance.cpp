#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

inline int remap_arom_type(int type)
{
	if (type<=400) return 1;
	if (type==401) return 2;
	if (type==402) return 3;
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

class Point
{
public:
	Point() {
		this->x=0.0;
		this->y=0.0;
		this->z=0.0;
	}
	Point(double xx, double yy, double zz) {
		this->x=xx;
		this->y=yy;
		this->z=zz;
	}
	Point(double *xx) {
		this->x=xx[0];
		this->y=xx[1];
		this->z=xx[2];
	}
	~Point() {}

	double x,y,z;
};

double distance(const double *x1, const double *x2, const double *prd);
double *unit_vector(const double *x1, const double *x2, const double *prd);
double *midpoint(const double *x1, const double *x2, const double *prd);
double dot(const double *x1, const double *x2, double a);
double *cross(const double *x1, const double *x2);
void calc_midpoints(Atom *atoms, int n_atoms, double *prd, Point *mid_points);
void normalize(double *x);
double prd_coord_add(double *x, double *xc, const double *x0, const double *prd);
inline void prd_gather(const double *x1, const double *x2, const double *x3, const double *x4, const double *prd, double **x);
inline void move_point_to_box(Point p, double A[3][2], double *prd);
void output_coordinates(FILE *outfile, Atom *atoms, int natoms, double A[3][2]);

int main(int argc, const char* argv[])
{
	if (argc<3) {
		printf("\n> %s Input_file Output_file [start] [end] [freq]\n\n", argv[0]);
		return -1;
	}

	// Read input/output file names
	const char *infile=argv[1];
	const char *outfile=argv[2];

	// Read the rest of the arguments
	int start=0, end=-1, freq=1;
	if (argc>3) start = atoi(argv[3]);
	if (argc>4) end = atoi(argv[4]);
	if (argc>5) freq = atoi(argv[5]);
	
	// Open input file
	FILE *input = fopen(infile, "r");
	if (!input) {
		printf("\nERROR: Input file not found!\n\n");
		return -1;
	}
	
	// Open output file
	FILE *output = fopen(outfile, "w");
	if (!output) {
		printf("\nERROR: Output file not found!\n");
		return -1;
	}
	fprintf(output, "# Instantaneous distances between DNA base-pairs \n");

	// Read input file and analyze frames one by one
	int i,j, atom_id, atom_type, nstr;
	double x, y, z;
	Atom *atoms=NULL;
	Point *mid_points=NULL;
	double *sum_dis=NULL;
	double *sum_sq_dis=NULL;
	double dis;
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
				if (step>=start && (end==-1 || step<=end) && i_frame%freq==0) {
					if (i_frame%1000==0)
						printf("Processing frame # %d\n", step);

					if (atoms==NULL) atoms = (Atom*) malloc(n_atoms*sizeof(Atom));

					int i_atom=0;
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
						atom_type = remap_arom_type(atom_type);

						if (atom_type!=1) break;

						x = atof(str[2]);
						y = atof(str[3]);
						z = atof(str[4]);


						x = (A[0][1] - A[0][0])*x + A[0][0];
						y = (A[1][1] - A[1][0])*y + A[1][0];
						z = (A[2][1] - A[2][0])*z + A[2][0];
						
						i_atom++;
						atoms[i].set(i_atom, atom_type, x, y, z);
					}
					n_atoms = i_atom;
					

					prd[0] = (A[0][1] - A[0][0]);
					prd[1] = (A[1][1] - A[1][0]);
					prd[2] = (A[2][1] - A[2][0]);

					if (mid_points==NULL) {
						mid_points = new Point[n_atoms/2];
						sum_dis = (double*)calloc(n_atoms/4,sizeof(double));
						sum_sq_dis = (double*)calloc(n_atoms/4,sizeof(double));
//						atoms = (Atom*)realloc(atoms, 3*n_atoms*sizeof(Atom)/2);
					}

					calc_midpoints(atoms, n_atoms, prd, mid_points);

/*					for (i=0;i<n_atoms/2;i++) {
						i_atom++;
						atom_type = 2;
						move_point_to_box(mid_points[i], A, prd);
						atoms[i+n_atoms].set(i_atom, atom_type, mid_points[i].x, mid_points[i].y, mid_points[i].z);
					}
					n_atoms = i_atom;*/

//					output_coordinates(output, atoms, n_atoms, A);

					for (i=0;i<n_atoms/4;i++) {
						double x1[]={mid_points[i].x, mid_points[i].y, mid_points[i].z};
						double x2[]={mid_points[i+n_atoms/4].x, mid_points[i+n_atoms/4].y, mid_points[i+n_atoms/4].z};
						dis = distance(x1, x2, prd);
						fprintf(output, "%.4f ", dis);
						sum_dis[i] += dis;
						sum_sq_dis[i] += dis*dis;
					}
					fprintf(output, "\n");

					n_frame++;
				}
				if (step>=start) i_frame++;
			}
		}
	}

	char outfileavg[256];
	strcpy(outfileavg, "avg_");
	strcat(outfileavg, outfile);
	FILE *output_avg = fopen(outfileavg,"w");
	fprintf(output_avg, "# Average distances between DNA base-pairs \n");
	for (i=0;i<n_atoms/4;i++) {
		sum_dis[i] /= n_frame;
		fprintf(output_avg, "%.4f ", sum_dis[i]);
	}
	fprintf(output_avg, "\n");

	double dev;
	fprintf(output_avg, "# Standard Deviation\n");
	for (i=0;i<n_atoms/4;i++) {
                sum_sq_dis[i] /= n_frame;
		dev = sqrt(sum_sq_dis[i] - sum_dis[i]*sum_dis[i]);
                fprintf(output_avg, "%.4f ", dev);
        }
        fprintf(output_avg, "\n");

	printf("Number of frames %d\n", n_frame);
	printf("Number of atoms %d\n", n_atoms);
	
	fclose(input);
	fclose(output);
	fclose(output_avg);

	if (atoms!=NULL) free(atoms);
	if (mid_points!=NULL) delete [] mid_points;
	if (sum_dis!=NULL) free(sum_dis);

	return 0;
}

inline void find_intersect(int i1, int i2, int i3, int i4, const Atom *atoms, const double *prd, const double *n1, const double *n2, double *x, double *x0=NULL)
{
	int i;
	double *n3, *xm1, *xm2, xc[3], p1, p2, p3, d, **xa;
	n3 = cross(n1, n2); // n3 = [n1 x n2]
	normalize(n3); // n3/|n3|

	xa = new double*[4];
	for (i=0;i<4;i++) xa[i] = new double[3];

	prd_gather(atoms[i1].x, atoms[i2].x, atoms[i3].x, atoms[i4].x, prd, xa);

	xm1 = midpoint(xa[0], xa[1], prd); // xm1 = (x1+x2)/2
	xm2 = midpoint(xa[2], xa[3], prd); // xm2 = (x3+x4)/2

	p1 = dot(n1, xm1, -1.0); // p1 = -(n1.xm1)
	p2 = dot(n2, xm2, -1.0); // p2 = -(n2.xm2)
	p3 = dot(n3, xa[0], -1.0); // p3 = -(n3.x1)

	d = (n1[2]*n2[1]-n1[1]*n2[2])*n3[0]+(n1[0]*n2[2]-n1[2]*n2[0])*n3[1]+(n1[1]*n2[0]-n1[0]*n2[1])*n3[2];
	xc[0] = -(p1*(n2[2]*n3[1]-n2[1]*n3[2])+p2*(n1[1]*n3[2]-n1[2]*n3[1])+p3*(n1[2]*n2[1]-n1[1]*n2[2]))/d;
	xc[1] = -(p1*(n2[0]*n3[2]-n2[2]*n3[0])+p2*(n1[2]*n3[0]-n1[0]*n3[2])+p3*(n1[0]*n2[2]-n1[2]*n2[0]))/d;
	xc[2] = -(p1*(n2[1]*n3[0]-n2[0]*n3[1])+p2*(n1[0]*n3[1]-n1[1]*n3[0])+p3*(n1[1]*n2[0]-n1[0]*n2[1]))/d;
	
	if (x0==NULL) {
		x[0] = xc[0];
		x[1] = xc[1];
		x[2] = xc[2];
	} else {
		prd_coord_add(x, xc, x0, prd);
	}

	for (i=0;i<4;i++) delete [] xa[i];
	delete [] xa;
}

void calc_midpoints1(Atom *atoms, int n_atoms, double *prd, Point *mid_points)
{
	int i,j,i1,i2,im;
        int n_dna=2;
        int nbp=n_atoms/(2*n_dna);
        double *xm;

        for (i=0;i<n_dna;i++) {
                for (j=0;j<nbp;j++) {
                        i1 = 2*nbp*i+j;
                        i2 = 2*nbp*i+nbp+j;
			im = i*nbp+j;
			
			xm = midpoint(atoms[i1].x, atoms[i2].x, prd); // xm1 = (x1+x2)/2
			
			mid_points[im].x = xm[0];
                        mid_points[im].y = xm[1];
                        mid_points[im].z = xm[2];
		}
	}			
}

void calc_midpoints(Atom *atoms, int n_atoms, double *prd, Point *mid_points)
{
	int i,j,k,i1,i2,i3,i4,im;
	int n_dna=2;
	int nbp=n_atoms/(2*n_dna);
	double *n1, *n2; 
	double x[3], x0[3], dn, max;
	int nx=0,k_max;

	for (i=0;i<n_dna;i++) {
		for (j=0;j<nbp;j++) {
			i1 = 2*nbp*i+j;
			i2 = 2*nbp*i+nbp+j;
			im = i*nbp+j;
			nx = 0; max = 0.0;
			
			n1 = unit_vector(atoms[i1].x, atoms[i2].x, prd); // n1 = (x2 - x1)/|x2-x1|

			x[0] = x[1] = x[2] = 0.0;
			for (k=-3;k<=3;++k) {
				if (k==0 || j+k<0 || j+k>=nbp) continue;
				
				i3 = 2*nbp*i+j+k;
                                i4 = 2*nbp*i+nbp+j+k;
				
				n2 = unit_vector(atoms[i3].x, atoms[i4].x, prd); // n2 = (x4 - x3)/|x4-x3|
				dn = 1-dot(n1,n2,1.0);
				if (dn>max) { max=dn; k_max=k; }

				if (dn>0.05) {
					find_intersect(i1, i2, i3, i4, atoms, prd, n1, n2, x, nx==0 ? NULL : x0);
					if (nx==0) { x0[0] = x[0]; x0[1] = x[1]; x0[2] = x[2]; }
					nx++;
				}
			}

			if (nx==0) {
				i3 = 2*nbp*i+j+k_max;
                                i4 = 2*nbp*i+nbp+j+k_max;

                                n2 = unit_vector(atoms[i3].x, atoms[i4].x, prd); // n2 = (x4 - x3)/|x4-x3|

				find_intersect(i1, i2, i3, i4, atoms, prd, n1, n2, x, nx==0 ? NULL : x0);
				nx++;
			}
			mid_points[im].x = x[0]/nx;
			mid_points[im].y = x[1]/nx;
			mid_points[im].z = x[2]/nx;
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

double *unit_vector(const double *x1, const double *x2, const double *prd)
{
	double half_prd[3],  vv;
	double *v = new double[3];

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

	vv = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);

	v[0] /= vv;
	v[1] /= vv;
	v[2] /= vv;
	
	return v;
}

double *midpoint(const double *x1, const double *x2, const double *prd)
{
	double half_prd[3], v[3]; 
	double *vm = new double[3];

	half_prd[0] = prd[0]/2;	
	half_prd[1] = prd[1]/2;	
	half_prd[2] = prd[2]/2;

	v[0] = x2[0] - x1[0];
	v[1] = x2[1] - x1[1];
	v[2] = x2[2] - x1[2];

	vm[0] = (x1[0] + x2[0])/2;
	vm[1] = (x1[1] + x2[1])/2;
	vm[2] = (x1[2] + x2[2])/2;

	if (v[0]<=-half_prd[0]) vm[0] += half_prd[0];
	else if (v[0]>half_prd[0]) vm[0] -= half_prd[0];
	if (v[1]<=-half_prd[1]) vm[1] += half_prd[1];
	else if (v[1]>half_prd[1]) vm[1] -= half_prd[1];
	if (v[2]<=-half_prd[2]) vm[2] += half_prd[2];
	else if (v[2]>half_prd[2]) vm[2] -= half_prd[2];

	return vm;
}

double prd_coord_add(double *x, double *xc, const double *x0, const double *prd)
{
	double half_prd[3], v[3];

	half_prd[0] = prd[0]/2;
	half_prd[1] = prd[1]/2;
	half_prd[2] = prd[2]/2;

	v[0] = xc[0] - x0[0];
	v[1] = xc[1] - x0[1];
	v[2] = xc[2] - x0[2];

	if (v[0]<=-half_prd[0]) xc[0] += prd[0];
	else if (v[0]>half_prd[0]) xc[0] -= prd[0];
	if (v[1]<=-half_prd[1]) xc[1] += prd[1];
	else if (v[1]>half_prd[1]) xc[1] -= prd[1];
	if (v[2]<=-half_prd[2]) xc[2] += prd[2];
	else if (v[2]>half_prd[2]) xc[2] -= prd[2];

	x[0] += xc[0];
	x[1] += xc[1];
	x[2] += xc[2];
}

double dot(const double *x1, const double *x2, double a=1.0)
{
	return a*(x1[0]*x2[0]+x1[1]*x2[1]+x1[2]*x2[2]);
}

double *cross(const double *x1, const double *x2)
{
	double *v = new double[3];

	v[0] = x1[1]*x2[2]-x1[2]*x2[1];
	v[1] = -x1[0]*x2[2]+x1[2]*x2[0];
	v[2] = x1[0]*x2[1]-x1[1]*x2[0];

	return v;
}

void normalize(double *x)
{
	double xm=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
	x[0] /= xm;
	x[1] /= xm;
	x[2] /= xm;
}

void output_coordinates(FILE *output, Atom *atoms, int natoms, double A[3][2])
{
	int i;
	double x,y,z;
	fprintf(output, "ITEM: TIMESTEP\n");
	fprintf(output, "0\n");
	fprintf(output, "ITEM: NUMBER OF ATOMS\n");
	fprintf(output, "%d\n", natoms);
	fprintf(output, "ITEM: BOX BOUNDS pp pp pp\n");
	fprintf(output, "%f %f\n", A[0][0], A[0][1]);
	fprintf(output, "%f %f\n", A[1][0], A[1][1]);
	fprintf(output, "%f %f\n", A[2][0], A[2][1]);
	fprintf(output, "ITEM: ATOMS id type xs ys zs\n");
	for (i=0;i<natoms;i++) {
		x = (atoms[i].x[0]-A[0][0])/(A[0][1]-A[0][0]);
		y = (atoms[i].x[1]-A[1][0])/(A[1][1]-A[1][0]);
		z = (atoms[i].x[2]-A[2][0])/(A[2][1]-A[2][0]);
		fprintf(output, "%d %d %f %f %f\n", atoms[i].id, atoms[i].type, x, y, z);
	}	
}

inline void move_point_to_box(Point p, double A[3][2], double *prd)
{
	if (p.x<A[0][0]) p.x += prd[0];
	else if (p.x>A[0][1]) p.x -= prd[0];
	if (p.y<A[1][0]) p.y += prd[1];
	else if (p.y>A[1][1]) p.y -= prd[1];
	if (p.z<A[2][0]) p.z += prd[2];
	else if (p.z>A[2][1]) p.z -= prd[2];
}

inline void prd_gather(const double *x1, const double *x2, const double *x3, const double *x4, const double *prd, double **x)
{
	double half_prd[3], v[3];

        half_prd[0] = prd[0]/2;
        half_prd[1] = prd[1]/2;
        half_prd[2] = prd[2]/2;

	x[0][0] = x1[0]; x[0][1] = x1[1]; x[0][2] = x1[2];
	x[1][0] = x2[0]; x[1][1] = x2[1]; x[1][2] = x2[2];
	x[2][0] = x3[0]; x[2][1] = x3[1]; x[2][2] = x3[2];
	x[3][0] = x4[0]; x[3][1] = x4[1]; x[3][2] = x4[2];
	
	for (int i=1;i<4;++i) {
		v[0] = x[i][0] - x[0][0];
		v[1] = x[i][1] - x[0][1];
		v[2] = x[i][2] - x[0][2];

		if (v[0]<=-half_prd[0]) x[i][0] += prd[0];
		else if (v[0]>half_prd[0]) x[i][0] -= prd[0];
		if (v[1]<=-half_prd[1]) x[i][1] += prd[1];
		else if (v[1]>half_prd[1]) x[i][1] -= prd[1];
		if (v[2]<=-half_prd[2]) x[i][2] += prd[2]; 
		else if (v[2]>half_prd[2]) x[i][2] -= prd[2];
	}
}

/* ----------------------------------------------------------------------
Copyright (2010) Aram Davtyan and Garegin Papoian

Papoian's Group, University of Maryland at Collage Park
http://papoian.chem.umd.edu/

The modified density depending on zim file is included

Last Update: 3/9/2012
------------------------------------------------------------------------- */

typedef struct WPV {
  double kappa;
  double kappa_sigma;
  double treshold;
  int n_wells, well_flag[5];
  double well_r_min[5], well_r_max[5];
  WPV() {}
  WPV(double k, double ks, double t, int nw, int *wflag, double *wr_min, double *wr_max): 
  	  kappa(k), kappa_sigma(ks), treshold(t), n_wells(nw)
  {
  	for (int i=0; i<n_wells; ++i) {
  		well_flag[i] = wflag[i];
  		well_r_min[i] = wr_min[i];
  		well_r_max[i] = wr_max[i];
  	}
  }
  WPV &operator = (const WPV &x) { 
  	kappa=x.kappa; kappa_sigma=x.kappa_sigma; treshold=x.treshold; n_wells=x.n_wells;
  	
  	for (int i=0; i<n_wells; ++i) {
  		well_flag[i] = x.well_flag[i];
  		well_r_min[i] = x.well_r_min[i];
  		well_r_max[i] = x.well_r_max[i];
  	}
  	
  	return *this; 
  }
} _WPV;

typedef struct TBV {
  double energy;
  double force;
  TBV() : energy(0.0), force(0.0) {}
  TBV(double ee, double ff): energy(ee), force(ff) {}
  TBV &operator = (const TBV &x) {
    energy=x.energy;
    force=x.force;
  
    return *this;
  }
} _TBV;

//=============================================================================================//

template <typename T, typename U>
class cP_AP {
public:
	cP_AP(int n, int m, int *indicator, U *lclass);
	~cP_AP();
	
	inline T &nu(int i, int j);
	inline T &prd_nu(int i, int j);
	void compute(int i, int j);
private:
	int n, m;
	T **v_nu;
	T **v_prd_nu;
	int **g;
	int *ind;
	U *lc;
	int drmax, drmin;
};

template <typename T, typename U>
cP_AP<T, U>::cP_AP(int nn, int mm, int *indicator, U *lclass)
{
	n = nn;
	m = mm;
	ind = indicator;
	lc = lclass;

	drmax = lc->P_AP_cut + int(8.0*2.302585/lc->P_AP_pref) + 1;
	drmin = lc->P_AP_cut - int(8.0*2.302585/lc->P_AP_pref) - 1;

	v_nu = new T*[n];
	v_prd_nu = new T*[n];
	g = new int*[n];

	for (int i=0;i<n;++i) {
		v_nu[i] = new T[m];
		v_prd_nu[i] = new T[m];
		g[i] = new int[m];
		for (int j=0;j<m;++j) g[i][j] = -1;
	}
}

template <typename T, typename U>
cP_AP<T, U>::~cP_AP()
{
	for (int i=0;i<n;++i) {
		delete [] v_nu[i];
		delete [] v_prd_nu[i];
		delete [] g[i];
	}

	delete [] v_nu;
	delete [] v_prd_nu;
	delete [] g;
}

template <typename T, typename U>
inline T &cP_AP<T, U>::nu(int i, int j)
{
	if ( (ind && g[i][j]!=*ind) || (!ind && g[i][j]!=1) ) {
		compute(i, j);
		if (ind) g[i][j] = *ind;
		else g[i][j] = 1;
	}

	return v_nu[i][j];
}

template <typename T, typename U>
inline T &cP_AP<T, U>::prd_nu(int i, int j)
{
	if ( (ind && g[i][j]!=*ind) || (!ind && g[i][j]!=1) ) {
		compute(i, j);
		if (ind) g[i][j] = *ind;
		else g[i][j] = 1;
	}

	return v_prd_nu[i][j];
}

template <typename T, typename U>
void cP_AP<T, U>::compute(int i, int j)
{
	T dx[3], dr, th;

	dx[0] = lc->xca[i][0] - lc->xca[j][0];
	dx[1] = lc->xca[i][1] - lc->xca[j][1];
	dx[2] = lc->xca[i][2] - lc->xca[j][2];
	
	dr = sqrt(pow(dx[0],2) + pow(dx[1],2) + pow(dx[2],2));

	if (dr>drmax) {
		v_nu[i][j] = 0.0;

		v_prd_nu[i][j] = 0.0;
	} else if(dr<drmin) {
		v_nu[i][j] = 1.0;

		v_prd_nu[i][j] = 0.0;
	} else {
		th = tanh(lc->P_AP_pref*(lc->P_AP_cut - dr));

		v_nu[i][j] = 0.5*(1+th);

		v_prd_nu[i][j] = 0.5*lc->P_AP_pref*(1-pow(th,2))/dr;
	}
}

//=============================================================================================//

template <typename T, typename U>
class cR {
public:
	cR(int n, int m, int *indicator, U *lclass);
	~cR();
	
	inline T &rNO(int i, int j);
	inline T &rHO(int i, int j);
private:
	int n, m;
	T **v_rNO;
	T **v_rHO;
	int **gNO;
	int **gHO;
	int *ind;
	U *lc;
};

template <typename T, typename U>
cR<T, U>::cR(int nn, int mm, int *indicator, U *lclass)
{
	n = nn;
	m = mm;
	ind = indicator;
	lc = lclass;

	v_rNO = new T*[n];
	v_rHO = new T*[n];
	gNO = new int*[n];
	gHO = new int*[n];

	for (int i=0;i<n;++i) {
		v_rNO[i] = new T[m];
		v_rHO[i] = new T[m];
		gNO[i] = new int[m];
		gHO[i] = new int[m];
		for (int j=0;j<m;++j) { gNO[i][j] = -1; gHO[i][j] = -1; }
	}
}

template <typename T, typename U>
cR<T, U>::~cR()
{
	for (int i=0;i<n;++i) {
		delete [] v_rNO[i];
		delete [] v_rHO[i];
		delete [] gNO[i];
		delete [] gHO[i];
	}

	delete [] v_rNO;
	delete [] v_rHO;
	delete [] gNO;
	delete [] gHO;
}

template <typename T, typename U>
inline T &cR<T, U>::rNO(int i, int j)
{
	if ( (ind && gNO[i][j]!=*ind) || (!ind && gNO[i][j]!=1) ) {
		v_rNO[i][j] = sqrt( pow(lc->xo[i][0] - lc->xn[j][0], 2) +
				    pow(lc->xo[i][1] - lc->xn[j][1], 2) +
				    pow(lc->xo[i][2] - lc->xn[j][2], 2) );

		if (ind) gNO[i][j] = *ind;
		else gNO[i][j] = 1;
	}

	return v_rNO[i][j];
}

template <typename T, typename U>
inline T &cR<T, U>::rHO(int i, int j)
{
	if ( (ind && gHO[i][j]!=*ind) || (!ind && gHO[i][j]!=1) ) {
		v_rHO[i][j] = sqrt( pow(lc->xo[i][0] - lc->xh[j][0], 2) +
			  	    pow(lc->xo[i][1] - lc->xh[j][1], 2) +
			  	    pow(lc->xo[i][2] - lc->xh[j][2], 2) );

		if (ind) gHO[i][j] = *ind;
		else gHO[i][j] = 1;
	}

	return v_rHO[i][j];
}

//=============================================================================================//

template <typename T, typename U>
class cWell {
public:
	cWell(int n, int m, int nw, const WPV &par, int *indicator, U *lclass);
	~cWell();
	
	inline T &theta(int i, int j, int i_well);
	inline T &prd_theta(int i, int j, int i_well);
	inline T &sigma(int i, int j);
	inline T &H(int i);
	inline T &prd_H(int i);
	inline T &ro(int i);
	
	void compute_theta(int i, int j, int i_well);
	void compute_sigma(int i, int j);
	void compute_H(int i);
	void compute_ro(int i);
	WPV par;

private:
	int n, m, nw;
	T ***v_theta;
	T ***v_prd_theta;
	T **v_sigma;
	T *v_H;
	T *v_prd_H;
	T *v_ro;
	int ***gTheta;
	int **gSigma;
	int *gH;
	int *gRo;
	int *ind;
	U *lc;
};

template <typename T, typename U>
cWell<T, U>::cWell(int nn, int mm, int ww, const WPV &p, int *indicator, U *lclass)
{
	n = nn;
	m = mm;
	nw = ww;
	ind = indicator;
	par = p;
	lc = lclass;

	v_sigma = new T*[n];
	v_H = new T[n];
	v_prd_H = new T[n];
	v_ro = new T[n];
	gSigma = new int*[n];
	gH = new int[n];
	gRo = new int[n];

	for (int i=0;i<n;++i) {
		gH[i] = -1;
		gRo[i] = -1;
		
		v_sigma[i] = new T[m];
		gSigma[i] = new int[m];

		for (int j=0;j<m;++j) {
			gSigma[i][j] = -1;
		}
	}

	v_theta = new T**[nw];
	v_prd_theta = new T**[nw];
	gTheta = new int**[nw];

	for (int k=0;k<nw;++k) {		
		v_theta[k] = new T*[n];
		v_prd_theta[k] = new T*[n];
		gTheta[k] = new int*[n];

		for (int i=0;i<n;++i) {
			v_theta[k][i] = new T[m];
			v_prd_theta[k][i] = new T[m];
			gTheta[k][i] = new int[m];
			for (int j=0;j<m;++j) gTheta[k][i][j] = -1;
		}
	}
}

template <typename T, typename U>
cWell<T, U>::~cWell()
{
	for (int k=0;k<nw;++k) {
		for (int i=0;i<n;++i) {
			delete [] v_theta[k][i];
			delete [] v_prd_theta[k][i];
			delete [] gTheta[k][i];
		}

		delete [] v_theta[k];
		delete [] v_prd_theta[k];
		delete [] gTheta[k];
	}

	for (int i=0;i<n;++i) {
		delete [] v_sigma[i];
		delete [] gRo[i];
	}

	delete [] v_theta;
	delete [] v_prd_theta;
	delete [] v_sigma;
	delete [] v_H;
	delete [] v_prd_H;
	delete [] v_ro;
	delete [] gTheta;
	delete [] gSigma;
	delete [] gH;
	delete [] gRo;
}

template <typename T, typename U>
inline T &cWell<T, U>::theta(int i, int j, int i_well)
{
	if ( (ind && gTheta[i_well][i][j]!=*ind) || (!ind && gTheta[i_well][i][j]!=1) ) {
		compute_theta(i, j, i_well);
		if (ind) gTheta[i_well][i][j] = gTheta[i_well][j][i] = *ind;
		else gTheta[i_well][i][j] = gTheta[i_well][j][i] = 1;
	}

	return v_theta[i_well][i][j];
}

template <typename T, typename U>
inline T &cWell<T, U>::prd_theta(int i, int j, int i_well)
{
	if ( (ind && gTheta[i_well][i][j]!=*ind) || (!ind && gTheta[i_well][i][j]!=1) ) {
		compute_theta(i, j, i_well);
		if (ind) gTheta[i_well][i][j] = gTheta[i_well][j][i] = *ind;
		else gTheta[i_well][i][j] = gTheta[i_well][j][i] = 1;
	}

	return v_prd_theta[i_well][i][j];
}

template <typename T, typename U>
inline T &cWell<T, U>::sigma(int i, int j)
{
	if ( (ind && gSigma[i][j]!=*ind) || (!ind && gSigma[i][j]!=1) ) {
		compute_sigma(i, j);
		if (ind) gSigma[i][j] = gSigma[j][i] = *ind;
		else gSigma[i][j] = gSigma[j][i] = 1;
	}

	return v_sigma[i][j];
}

template <typename T, typename U>
inline T &cWell<T, U>::H(int i)
{
	if ( (ind && gH[i]!=*ind) || (!ind && gH[i]!=1) ) {
		compute_H(i);
		if (ind) gH[i] = *ind;
		else gH[i] = 1;
	}

	return v_H[i];
}

template <typename T, typename U>
inline T &cWell<T, U>::prd_H(int i)
{
	if ( (ind && gH[i]!=*ind) || (!ind && gH[i]!=1) ) {
		compute_H(i);
		if (ind) gH[i] = *ind;
		else gH[i] = 1;
	}

	return v_prd_H[i];
}

template <typename T, typename U>
inline T &cWell<T, U>::ro(int i)
{
	if ( (ind && gRo[i]!=*ind) || (!ind && gRo[i]!=1) ) {
		compute_ro(i);
		if (ind) gRo[i] = *ind;
		else gRo[i] = 1;
	}

	return v_ro[i];
}

template <typename T, typename U>
void cWell<T, U>::compute_theta(int i, int j, int i_well)
{
	T dx[3], rij, t_min, t_max;
	T *xi, *xj;
	
	int i_resno = lc->res_no[i]-1;
	int j_resno = lc->res_no[j]-1;
	
	if (lc->se[i_resno]=='G') xi = lc->xca[i];
	else xi = lc->xcb[i];
	if (lc->se[j_resno]=='G') xj = lc->xca[j];
	else xj = lc->xcb[j];

	dx[0] = xi[0] - xj[0];
	dx[1] = xi[1] - xj[1];
	dx[2] = xi[2] - xj[2];
	
	rij = sqrt(pow(dx[0],2) + pow(dx[1],2) + pow(dx[2],2));
	
	t_min = tanh( par.kappa*(rij - par.well_r_min[i_well]) );
	t_max = tanh( par.kappa*(par.well_r_max[i_well] - rij) );
	v_theta[i_well][i][j] = 0.25*(1.0 + t_min)*(1.0 + t_max);
	v_theta[i_well][j][i] = v_theta[i_well][i][j];
	v_prd_theta[i_well][i][j] = par.kappa*v_theta[i_well][i][j]*(t_max - t_min)/rij;
	v_prd_theta[i_well][j][i] = v_prd_theta[i_well][i][j];
}

template <typename T, typename U>
void cWell<T, U>::compute_sigma(int i, int j)
{	
	v_sigma[i][j] = v_sigma[j][i] = H(i)*H(j);
}

template <typename T, typename U>
void cWell<T, U>::compute_H(int i)
{	
	T g, th;
	
	g = par.kappa_sigma*(ro(i) - par.treshold);
	th = tanh(g);
	v_H[i] = 0.5*( 1.0 -  th);
	v_prd_H[i] = -0.5*par.kappa_sigma*( 1.0 - th*th );
}

template <typename T, typename U>
void cWell<T, U>::compute_ro(int i)
{
  int j;
  
  v_ro[i] = 0;
  for (j=0;j<lc->nn;++j) {
  	if (lc->res_info[j]==lc->OFF) continue;
  	
  	if ( lc->chain_no[i]!=lc->chain_no[j] || abs(lc->res_no[j] - lc->res_no[i])>1 )
  		v_ro[i] += theta(i, j, 0);
  }
  
// add new density which depend on z (if memb potential is on)
  if ( lc->memb_flag){
   if (lc->z_res[i]==2){
     v_ro[i] +=lc->memb_dens_offset;
  }
}


//  for (j=0;j<i-1;++j) v_ro[i] += theta(i, j, 0);
//  for (j=i+2;j<n;++j) v_ro[i] += theta(i, j, 0);
}

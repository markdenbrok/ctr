#include<fftw3.h>
#include"gaussfit.h"

#define True 1
#define False 0

#define AU  1.4959787068E8 //AU in km
#define G 0.00430237 // G in (km/s)^2 pc/Msun
#define RADEG 57.29578 // degrees per radian
#define RA2DEG 57.29578 // degrees per radian
#define pc2km  3.0856776e+13 // (km per parsec)


struct multigaussexp {
  double *area;
  double *sigma;
  double *q;
  double *pa;
  int ntotal;
};



struct mdb_conv_struc
{
  double *a; 
  double *psf;  
  double *c; 
  double *f;
  //double *fs;
  double **awrap;
  double **psfwrap;
  double **cwrap;
  double **fwrap;

  fftw_complex *A;  // model
  fftw_complex *PSF; 
  fftw_complex *C;  // PSF convolved model

  fftw_plan pa;     // model FFT
  fftw_plan ppsf; 
  fftw_plan pcinv;  // PSF convolved model inverse FFT

  int M;            // the size of the thing
  int N;
  int initialized;
  int psf_xoffset;
  int psf_yoffset;
  double **vfield, **sfield, ***Z;
  double *vdata, *sdata, *Zdata;
  int nvel;
  double vsample;
  //double Z[40][68][69];

  double *fluxarr;
  double **outsigma;
  double **outvel;
  double *outsigmadata;
  double *outveldata;
  struct gauss_data gdat;
  double  *gpars;
};




struct image {
  double *pixels;
  double **pixelwrap;
  int nx, ny;
  int zeropoint;
  int extcoeff;
};



struct funcparams {
  struct multigaussexp potential;
  double R;
};


// In kepler.c **********************

//galaxy distance is in pc
//inclination is in radians? maybe better in degrees
//
void keplerian_disc(double, struct multigaussexp *, double *, int,double *,double,double , double,int ,int,double,double, double,double *, double **);


//********************************

// In potential.c

// The pot mge should be spherical and with sigma in pc, dens = mass of gauss in solar mass
double *dphi(double *, struct multigaussexp *, double, int);
double qpint1d(double, double, double, double,struct funcparams *);
double qpfunc(double, void *);
double potential_integrand(double, struct multigaussexp *, double);
// this is for the hot disc
//everything in pc and km/s
double *dsigma2drhorho (struct multigaussexp *,struct multigaussexp *, double *, int);


// r should be in pc, pot->dens=mass in Msun, sig=[pc],
double *v_azi2(double *, struct multigaussexp *, double, int);


//********************************
// in line_model.c

double lsf(double,double, double, double);


//
//mbh in solar mass
//galaxy distance in Mpc

// the psf should be in I0, arcsec.
// Normalization will be taken care of.
// all the mges in solarmass/pc^2-arcsec
// rad should be in arcsec
// pa should be in radians
// xc, yc in pix
// platescale in arcsec
// cs has to be made with initialized FALSE the first time.
// outimage has to be present as well.
void line_modeling(double, struct multigaussexp *,struct multigaussexp *,double *, double *, int, double, double, double,int,int,double ,double ,struct multigaussexp *, double,double *, struct mdb_conv_struc *, double **, int, double, struct image *,double, double, double, double);


void initialize_convolution(struct mdb_conv_struc *, int, double);
void initialize_velocity_arrays(int, int, int, double ****, double ***, double ***); /* This module declares the arrays in that will contain the velocity field, sigma field and Z (the spectral cube) */



void setup_psf_stuff(struct multigaussexp *, struct mdb_conv_struc *, double );
void setup_flux_stuff(struct multigaussexp *, struct mdb_conv_struc *, double , double , double);
void setup_flux_stuff_inputimage(struct image *, struct mdb_conv_struc *, double ,double, double, double, double);
void clean_up_convolution(struct mdb_conv_struc *);
void mdb_fftw_convolve_2d(struct mdb_conv_struc *);
void rebin(struct mdb_conv_struc *,double **,int,int);




void setup_conv_stuff(struct mdb_conv_struc *,int,int,struct multigaussexp *,double); 
/* This module calculates the size of the convolution arrays (M,N) in the convolution structure 
 */


//****************************8
//read_mge.c

void read_mge(char *, int , struct multigaussexp *);
void free_mge(struct multigaussexp *);


//*****************************
// write_fits.c

int fits_write_wrapper(double **,int, int,char *);
int write_fits(char *,struct image *);

//***************************
// read_fits.c
struct image read_fits(char *);

//*****************************
// print_mge.c
void print_mge(struct multigaussexp *);



//***************************
// calc_sigma_field.c
void calc_sigma_field(struct multigaussexp *,double,double, double, int, int, double **, double);



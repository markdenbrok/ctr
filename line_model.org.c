/************************************************************************
 * This code was written by Mark den Brok (Utah 2014)                   *
 * based on the IDL code of Cappellari                                  *
 * however, the convolution has changed (bigger size)                   *
 * denbrok@physics.utah.edu                                             *
 ************************************************************************/
//#define DEBUG

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <fftw3.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_erf.h>
#include "tr.h"
#define fftw_real double
#define Pi 3.14159265



#define oversample 10


double lsf(double v,double v_pix, double sigma, double dv){
  double res;
  res=gsl_sf_erf( (v+0.5*dv-v_pix)/(sqrt(2.0)*sigma)) - gsl_sf_erf( (v-0.5*dv-v_pix)/(sqrt(2.0)*sigma));  
  return 0.5/dv* res;
}



//mbh in solar mass
//galaxy distance in Mpc

// the psf should be in I0, arcsec.
// Normalization will be taken care of.
void line_modeling(double mbh, struct multigaussexp *pot, struct multigaussexp *disc,struct multigaussexp *sigma,double *rad, double *pa_varia, int nrad, double xc, double yc, double plateScale,int nxPix,int nyPix,double vSist,double sigma_inst,struct multigaussexp *psf, double galaxy_distance,double *incl, struct mdb_conv_struc *cs, double **outimage){

  double dx, distance, xc_os, yc_os;
  double *rad_pc;
  int nx, ny, i,j;
  struct multigaussexp pot_pc, disc_pc, sigma_pc; 
  // We take oversampling of factor oversample
  nx=nxPix*oversample;
  ny=nyPix*oversample;
  dx=plateScale/oversample; // This is still in arcsecond/pixel
  distance = galaxy_distance*(1.0e6); //to pc

  
 
  // prepare the image to fill with values
  if (cs->initialized==False){
    
    setup_conv_stuff(cs,nx,ny, psf,dx);
  
    initialize_convolution(cs);
  
    setup_psf_stuff(psf, cs, dx);
  
    
    }
  else {
    //printf("Convolution already initialized\n");
    // clean the arrays from last time
    for (i=0; i<cs->M; i++){
      for (j=0; j<cs->N; j++){
	cs->awrap[i][j]=0.;
	cs->cwrap[i][j]=0.;
      }
    }
  }
  
  // The pixel coordinate falls at the centre of the pixel
  xc_os = xc*((double) oversample) + ((double) oversample-1.0)/(2.0)+ cs->psf_xoffset;
  yc_os = yc*((double) oversample) + ((double) oversample-1.0)/(2.0)+ cs->psf_yoffset;
 

  setup_flux_stuff(disc,cs,dx,xc_os,yc_os);
#ifdef DEBUG
  fits_write_wrapper(cs->fwrap,cs->M,cs->N,"f.fits");
#endif
  //convert distance to pc
  distance=galaxy_distance*(1.e6);
  rad_pc=(double *) malloc((nrad)*sizeof(double));
  // if the outer raidus is smaller than the image, we extrapolate everything
  
  for (i=0; i<nrad; i++){ 
    rad_pc[i]=rad[i]*distance*M_PI/180.0/3600.0; // convert radii in pc
  } 
 


  //create new MGEs that are in pc and total mass
  pot_pc.area=(double *)malloc(pot->ntotal*sizeof(double));
  pot_pc.sigma=(double *)malloc(pot->ntotal*sizeof(double));
  pot_pc.ntotal=pot->ntotal;



  for(i=0; i<pot_pc.ntotal; i++){
    pot_pc.sigma[i]=pot->sigma[i]*distance*M_PI/180.0/3600.0; // sigma from arcsec to pc
    pot_pc.area[i]=pot->area[i]*2.0*M_PI*pot_pc.sigma[i]*pot_pc.sigma[i]*pot->q[i];
  }
  pot_pc.q=pot->q;

  disc_pc.area=(double *)malloc(disc->ntotal*sizeof(double));
  disc_pc.sigma=(double *)malloc(disc->ntotal*sizeof(double));
  
  disc_pc.ntotal=disc->ntotal;
 
  for(i=0; i<disc_pc.ntotal; i++){
    disc_pc.sigma[i]=disc->sigma[i]*distance*M_PI/180.0/3600.0; // sigma from arcsec to pc
    disc_pc.area[i]=disc->area[i]*2.0*M_PI*disc_pc.sigma[i]*disc_pc.sigma[i]*disc->q[i];
  
  }
  disc_pc.q=disc->q;

  sigma_pc.area=(double *)malloc(sigma->ntotal*sizeof(double));
  sigma_pc.sigma=(double *)malloc(sigma->ntotal*sizeof(double));
  sigma_pc.ntotal=sigma->ntotal;
 
  for(i=0; i<sigma_pc.ntotal; i++){
    sigma_pc.sigma[i]=sigma->sigma[i]*distance*M_PI/180.0/3600.0; // sigma from arcsec to pc
    sigma_pc.area[i]=sigma->area[i]*2.0*M_PI*sigma_pc.sigma[i]*sigma_pc.sigma[i]*sigma->q[i];
  }
  sigma_pc.q=sigma->q;
  //fill the image
  
  keplerian_disc(mbh, &pot_pc, &disc_pc, rad_pc, nrad, pa_varia, xc_os, yc_os, dx,cs->M,cs->N, vSist, sigma_inst, distance, incl, cs->awrap);
  
#ifdef DEBUG 
 fits_write_wrapper(cs->awrap,cs->M,cs->N,"a1.fits");
#endif
  // we weigh the image by the flux
  for (i=0; i<cs->M; i++){
    for (j=0; j<cs->N; j++){
      cs->awrap[i][j]*=cs->fwrap[i][j];
      //cs->awrap[i][j]=cs->fwrap[i][j]; // as a test
    }
  }
#ifdef DEBUG
  fits_write_wrapper(cs->awrap,cs->M,cs->N,"a2.fits");
#endif  
  
  //do the convolution
 
  mdb_fftw_convolve_2d(cs);

#ifdef DEBUG
  fits_write_wrapper(cs->cwrap,cs->M,cs->N,"c.fits");
#endif
  //rebin the model to match the observations
  
  // this is a test for debuggin purposes
#ifdef DEBUG 
  double **testim;
  testim=(double **) malloc(cs->M*sizeof(double *));
  for (i=0; i<cs->M; i++){
    testim[i]=(double *) malloc(cs->N*sizeof(double));
    for (j=0; j<cs->N; j++){
      testim[i][j]=cs->cwrap[i][j]/cs->fswrap[i][j];
    }
  } 
  fits_write_wrapper(testim,cs->M,cs->N,"testim.fits");
  for (i=0; i<cs->M; i++){
    free(testim[i]);
  }
#endif
  
 
  rebin(cs,outimage,nxPix,nyPix);
  
  

  //clean up

  free(rad_pc);

  free(sigma_pc.sigma);
  free(disc_pc.sigma);
  free(pot_pc.sigma);
  
  free(sigma_pc.area);
  free(disc_pc.area);
  free(pot_pc.area);
  

  //return with the new model.

  return;
}



void initialize_convolution(struct mdb_conv_struc *cs){
  int i, M,N;


 
  M=cs->M;
  N=cs->N;
  // First we initialize the pfs. We only use one block for this
  cs->psf = (double *) fftw_malloc(M*N * sizeof(double));
  cs->PSF = (fftw_complex*) fftw_malloc(M*(N/2+1) * sizeof(fftw_complex));
  cs->ppsf = fftw_plan_dft_r2c_2d(M, N, cs->psf, cs->PSF, FFTW_MEASURE);
  
 
  
  cs->psfwrap=(double **) malloc(M*sizeof(double *));

  //for (i=0; i<M; i++){
  //  cs->psfwrap[i]=&cs->psf[N*i+0];
  //}  
  // Note that the PSF has been screwed up by the plans


  // a denotes the model
  // c denotes the psf convolved model
  cs->a=(double *)fftw_malloc(M*N * sizeof(double));
  cs->c=(double *)fftw_malloc(M*N * sizeof(double));
  cs->f=(double *)fftw_malloc(M*N * sizeof(double));
  cs->fs=(double *)fftw_malloc(M*N * sizeof(double));
  cs->A=(fftw_complex *)fftw_malloc(M*(N/2+1) * sizeof(fftw_complex));
  cs->C=(fftw_complex *)fftw_malloc(M*(N/2+1) * sizeof(fftw_complex));
  cs->F=(fftw_complex *)fftw_malloc(M*(N/2+1) * sizeof(fftw_complex));
  cs->FS=(fftw_complex *)fftw_malloc(M*(N/2+1) * sizeof(fftw_complex));
  cs->awrap=(double **)malloc(M* sizeof(double *));
  cs->cwrap=(double **)malloc(M* sizeof(double *));
  cs->fwrap=(double **)malloc(M* sizeof(double *));
  cs->fswrap=(double **)malloc(M* sizeof(double *));

  for (i=0; i<M; i++){
    cs->psfwrap[i]=&cs->psf[N*i+0];
    cs->awrap[i]=&cs->a[N*i+0];
    cs->cwrap[i]=&cs->c[N*i+0];
    cs->fwrap[i]=&cs->f[N*i+0];
    cs->fswrap[i]=&cs->fs[N*i+0];
  } // this allows us to write psfwrap[M][N] & cetera; 
 
  cs->pa=fftw_plan_dft_r2c_2d(M, N, cs->a, cs->A, FFTW_MEASURE);
  cs->pcinv=fftw_plan_dft_c2r_2d(M, N, cs->C, cs->c, FFTW_MEASURE);
  cs->pf=fftw_plan_dft_r2c_2d(M, N, cs->f, cs->F, FFTW_MEASURE);
  cs->pfsinv=fftw_plan_dft_c2r_2d(M, N, cs->FS, cs->fs, FFTW_MEASURE);
  
  cs->initialized=True;


  return;
}



void setup_psf_stuff(struct multigaussexp *psf, struct mdb_conv_struc *cs, double dx){
  int i,j,k,M,N;
  double total_flux, dist;
  M=cs->M;
  N=cs->N;
 
  total_flux=0.;
  
  for (i=0; i<psf->ntotal; i++){
    total_flux+=psf->area[i]*2.*M_PI*pow(psf->sigma[i]/dx,2.0);
  }// the total flux
  //printf("Total flux: %lf\n",total_flux);

  //make sure that the array is emtpy
  //since fftw uses it for initializing the convolution
  for (i=0; i<M; i++){
    for (j=0; j<N; j++){
      cs->psfwrap[i][j]=0.;
    }
  }

  //we put the center at (0,0)
  //so that the psf is distributed over the four corners
  for (i=-M/4; i<M/4; i++){
    for (j=-N/4; j<N/4; j++){
      for (k=0; k<psf->ntotal; k++)
	{
	  dist=sqrt(pow((double) i,2.)+pow((double) j,2.));
	  cs->psfwrap[(i+M)%M][(j+N)%N]+=psf->area[k]/total_flux*exp(-0.5*pow( dist/(psf->sigma[k]/dx),2.0)); // this is an approximation, 
	  
	  // maybe should change to gsl_sf_erf (although with the oversampling
	  // there is probaly fairly little difference between it)
	}
      
    }
  }
  //cs->psfwrap[0][0]=1.; // for testing
#ifdef DEBUG
  fits_write_wrapper(cs->psfwrap,M,N,"p.fits");  
#endif
 
  // actual fft
  fftw_execute(cs->ppsf);
 
  return;
}


void setup_flux_stuff(struct multigaussexp *disc, struct mdb_conv_struc *cs, double dx, double xc, double yc){
  int i,j,k,M,N;
  //double total_flux;
  double dist;
  double xacc,yacc;
  M=cs->M;
  N=cs->N;

  //printf("Setting up flux stuff\n");
  //printf("%f \n",disc->q[0]);

  //printf("=========\n");
  // normalizing is not necessary
  //total_flux=0.;
  //for (i=0; i<psf->ntotal; i++){
  //  total_flux+=psf->area[i]*2.*M_PI*pow(psf->sigma[i]/dx,2.0);
  //}// the total flux
  

  //make sure that the array is emtpy
  //since fftw uses it for initializing the convolution
  //for (i=0; i<M; i++){
  //  for (j=0; j<N; j++){
  //    cs->fwrap[i][j]=0.;
  //  }
  //}
  //we put the center at xc,yc which should be the oversampled centers
  //

  for (i=0; i<M; i++){
    for (j=0; j<N; j++){
      for (k=0; k<disc->ntotal; k++)
	{
	  //printf("PA: %f\n",disc->pa[k]);
	  // first define xacc, yacc
	  // xacc = (i-xc) cos(pa) - (j-yc) sin(pa)
	  // yacc = (i-xc) sin(pa) + (j-yc) cos(pa)
	  // pa should be radians, defined in the mathematical way (x through y)
	  xacc = (i-xc)*cos(disc->pa[k]) - (j-yc)*sin(disc->pa[k]);
	  yacc = (i-xc)*sin(disc->pa[k]) + (j-yc)*cos(disc->pa[k]);
	  
	  dist=sqrt(pow(xacc,2.)+pow(yacc,2.)/pow(disc->q[k],2.0));
	  cs->fwrap[i][j]+=disc->area[k]*exp(-0.5*pow( dist/(disc->sigma[k]/dx),2.0)); // this is an approximation, 
	  // maybe should change to gsl_sf_erf (although with the oversampling
	  // there is probaly fairly little difference between it)
	}
      
    }
  }
  //printf("Writing line model\n");
  //fits_write_wrapper(cs->fwrap,M,N,"f2.fits");


  // we do no actual fft, since it is done somewhere else
  return;
  

}




void clean_up_convolution(struct mdb_conv_struc *cs){
  // destroy your plans!
  fftw_destroy_plan(cs->pa);
  fftw_destroy_plan(cs->ppsf);
  fftw_destroy_plan(cs->pcinv);
  fftw_destroy_plan(cs->pf);
  fftw_destroy_plan(cs->pfsinv);
  fftw_free(cs->a);
  fftw_free(cs->psf);
  fftw_free(cs->c);
  fftw_free(cs->f);
  fftw_free(cs->fs);
  
  free(cs->awrap);
  free(cs->psfwrap);
  free(cs->cwrap);
  free(cs->fwrap);
  free(cs->fswrap);
 
  
    
  fftw_free(cs->A);
  fftw_free(cs->C);
  fftw_free(cs->PSF);
  fftw_free(cs->F);
  fftw_free(cs->FS);
  
  cs->initialized=False;
  
  
  return;
  
}



void mdb_fftw_convolve_2d(struct mdb_conv_struc *cs){
  int M,N;
  int i, j, ij;
  double omega;
  fftw_complex *A, *B, *C, *F, *FS;
  fftw_real scale;
  M=cs->M;
  N=cs->N;
  scale = 1.0 / (M * N);
  omega=2.0*M_PI;
  A=cs->A;
  B=cs->PSF;
  C=cs->C;
  F=cs->F;
  FS=cs->FS;

  

  fftw_execute(cs->pa);
  fftw_execute(cs->pf);

  for (i=0; i<M; i++){
    for (j=0; j<(N/2+1); j++){
      ij = i*(N/2+1) + j;
      C[ij][0] = (A[ij][0] * B[ij][0] - A[ij][1] * B[ij][1]) * scale;
      C[ij][1] = (A[ij][0] * B[ij][1] + A[ij][1] * B[ij][0]) * scale;
      FS[ij][0] = (F[ij][0] * B[ij][0] - F[ij][1] * B[ij][1]) * scale;
      FS[ij][1] = (F[ij][0] * B[ij][1] + F[ij][1] * B[ij][0]) * scale;
    }
  }
  fftw_execute(cs->pcinv);  // the smoothed surface density weighted velocity
  fftw_execute(cs->pfsinv); // the smoothed surface density
  
#ifdef DEBUG  
  fits_write_wrapper(cs->fswrap,M,N,"fsmooth.fits");
#endif
  // now we have the smooth images, which we first have to rebin and then divide
  // in reality the multiplication with scale and the normalization of the PSF
  // are not necessary -- though I left it here because
  // more elegant and not too time-consuming.

  return;
}


void rebin(struct mdb_conv_struc *cs,double **outimage,int nx,int ny){
  int i,j,k,l;
  double q,r;
 
  for (i=0; i<nx; i++){
    for (j=0; j<ny; j++){
      q=0.;
      r=0.;
      for (k=0; k<oversample; k++){
	for (l=0; l<oversample; l++){
	  r+=cs->fswrap[i*oversample+k+cs->psf_xoffset][j*oversample+l+cs->psf_yoffset];
	  q+=cs->cwrap[i*oversample+k+cs->psf_xoffset][j*oversample+l+cs->psf_yoffset];
	}
      }
    
      outimage[i][j]=q/r;
      
    }
  }
  return;

}


void setup_conv_stuff(struct mdb_conv_struc *cs,int nx,int ny,struct multigaussexp *psf,double dx){
  int i, M,N;
  double maxsigma;
  

  maxsigma=psf->sigma[0];
  for (i=0; i<psf->ntotal; i++){
    if (psf->sigma[i]>maxsigma)
      maxsigma=psf->sigma[i];
  }
  maxsigma/=dx;
  // we use a width of 4 sigma or the image width as the maximum width of the psf.
  maxsigma*=4;
  if (maxsigma > nx) {
    M=2*nx;
    cs->psf_xoffset=nx/2;
  }
  else {M=nx+maxsigma;
    cs->psf_xoffset=maxsigma/2;
  }
  if (maxsigma > ny) {
    N=2*ny;
    cs->psf_yoffset=ny/2;
  }
  else {
    N=ny+maxsigma;
    cs->psf_yoffset=maxsigma/2;
  }

  
  // next, we see if this is a power of 2
  // otherwise, we find the nearest larger power of 2 
  // some bitwise stuff here
  // it seems that zero size arrays will be a problem
  // not checking for it though. Only idiots want to produce
  // zero-sized images.

  //cs->M=M;
  //cs->N=N;
  
  if (!(((M & (M - 1)) == 0) && (M>0)))
    {
      cs->M=2;
      while (M >>=1) cs->M*=2;
    }
  else{cs->M=M;}
  if (!(((N & (N - 1)) == 0) && (N>0)))
    {
      cs->N=2;
      while (N >>=1) cs->N*=2;
      
    }
  else{cs->N=N;}
  // that should do the trick.

  return;
  



}

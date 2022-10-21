/************************************************************************
 * This code was written by Mark den Brok (Utah 2014)                   *
 * based on the IDL code of Cappellari                                  *
 * however, the convolution has changed (bigger size)                   *
 * also, now complete line modeling for first moment                     *
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



#define oversample 2


double lsf(double v,double v_pix, double sigma, double dv){
  double res;
  res=gsl_sf_erf( (v+0.5*dv-v_pix)/(sqrt(2.0)*sigma)) - gsl_sf_erf( (v-0.5*dv-v_pix)/(sqrt(2.0)*sigma));  
  return 0.5/dv* res;
}



//mbh in solar mass
//galaxy distance in Mpc

// the psf should be in I0, arcsec.
// Normalization will be taken care of.

// hole_rad should be in arcsec
void line_modeling(double mbh, struct multigaussexp *pot,struct multigaussexp *sigma,double *rad, double *pa_varia, int nrad, double xc, double yc, double plateScale,int nxPix,int nyPix,double vSist,double sigma_inst,struct multigaussexp *psf, double galaxy_distance,double *incl, struct mdb_conv_struc *cs, double **outimage, int nvel, double dv, struct image* fluxmap, double inc_gal, double hole_rad, double hole_PA, double hole_q){
  
  double dx, distance, xc_os, yc_os;
  double *rad_pc;
  int nx, ny, i,j, k,l, m;
  double ***Z; // Z contains the  PSF convolved line model, size (M x N * nv)
  double q;
  double *velocities; // Vector containing the velocity axis (z axis of Z)
  struct multigaussexp pot_pc, sigma_pc; 
  double **outflux;
  
  // We take oversampling of factor oversample
  nx=nxPix*oversample;
  ny=nyPix*oversample;
  dx=plateScale/oversample; // This is still in arcsecond/pixel
  distance = galaxy_distance*(1.0e6); //to pc
  
  
  //for (i=0; i<nrad; i++) printf("%lf  \n",rad[i]);  
  
   
  // prepare the image to fill with values
  if (cs->initialized==False){
    
    setup_conv_stuff(cs,nx,ny, psf,dx);
  
    initialize_convolution(cs,nvel,dv);
  
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
  Z=cs->Z;
  //initialize_velocity_arrays(nvel,cs->M,cs->N,&Z, &vfield, &sfield);
  
  velocities=cs->gdat.x;
  //printdarr(velocities,nvel,"velocities");
  // The pixel coordinate falls at the centre of the pixel
  //xc_os = xc*((double) oversample) + ((double) oversample-1.0)/(2.0)+ cs->psf_xoffset;
  //yc_os = yc*((double) oversample) + ((double) oversample-1.0)/(2.0)+ cs->psf_yoffset;
  xc_os = (xc-1.0)*((double) oversample) + cs->psf_xoffset;
  yc_os = (yc-1.0)*((double) oversample) + cs->psf_yoffset;

  //setup_flux_stuff(disc,cs,dx,xc_os,yc_os);
  setup_flux_stuff_inputimage(fluxmap,cs,hole_rad/dx,xc_os,yc_os,hole_PA,hole_q);
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
  pot_pc.q=(double *)malloc(pot->ntotal*sizeof(double));
  pot_pc.ntotal=pot->ntotal;


  
  for(i=0; i<pot_pc.ntotal; i++){
    pot_pc.sigma[i]=pot->sigma[i]*distance*M_PI/180.0/3600.0; // sigma from arcsec to pc
    pot_pc.area[i]=pot->area[i]*2.0*M_PI*pot_pc.sigma[i]*pot_pc.sigma[i]*pot->q[i]; // For clarification: this is the mass in each MGE component, no longer I0
    pot_pc.q[i]=sqrt(pot->q[i]*pot->q[i]-cos(inc_gal)*cos(inc_gal))/sin(inc_gal);
  }
  // This is hardcoded: don't deproject these two. This is a test to see if 
  // using spherical symmetry for those 2 components makes a difference
  //pot_pc.q[0]=pot->q[0];
  //pot_pc.q[1]=pot->q[1];
  //pot_pc.q[2]=pot->q[2];



  sigma_pc.area=(double *)malloc(sigma->ntotal*sizeof(double));
  sigma_pc.sigma=(double *)malloc(sigma->ntotal*sizeof(double));
  sigma_pc.ntotal=sigma->ntotal;
 
  for(i=0; i<sigma_pc.ntotal; i++){
    sigma_pc.sigma[i]=sigma->sigma[i]*distance*M_PI/180.0/3600.0; // sigma from arcsec to pc
    sigma_pc.area[i]=sigma->area[i]*2.0*M_PI*sigma_pc.sigma[i]*sigma_pc.sigma[i]*sigma->q[i];
  }
  sigma_pc.q=sigma->q;
  //fill the image
  
  keplerian_disc(mbh, &pot_pc, rad_pc, nrad, pa_varia, xc_os, yc_os, dx,cs->M,cs->N, vSist, sigma_inst, distance, incl, cs->vfield);
  
  //printf("Calculated velocity model\n");

  calc_sigma_field(&sigma_pc,xc_os,yc_os,dx,cs->M,cs->N,cs->sfield,sigma_inst);

  //printf("Calculated dispersion model\n");

#ifdef DEBUG 
  fits_write_wrapper(cs->vfield,cs->M,cs->N,"a1.fits");
#endif

  for (k=0; k<nvel; k++){

 
   // we copy the vfield image to do convolution
   for (i=0; i<cs->M; i++){
     for (j=0; j<cs->N; j++){
       cs->awrap[i][j]=cs->fwrap[i][j]*lsf(cs->vfield[i][j], velocities[k],cs->sfield[i][j], cs->vsample);
       //cs->awrap[i][j]=lsf(cs->vfield[i][j], velocities[k],cs->sfield[i][j], vsample);
     }
   }
   //if (k==(nvel/2)){
   //  fits_write_wrapper(cs->awrap,cs->M,cs->N,"a2.fits");
   //}
   //if (k==(nvel-2)){
   //  fits_write_wrapper(cs->awrap,cs->M,cs->N,"a3.fits");
   // }
#ifdef DEBUG
  fits_write_wrapper(cs->awrap,cs->M,cs->N,"a2.fits");
#endif  
  
  
  
  //do the convolution
 
  mdb_fftw_convolve_2d(cs);

  //if (k==(nvel/2)){
  //   fits_write_wrapper(cs->cwrap,cs->M,cs->N,"c2.fits");
  // }
  // if (k==(nvel-2)){
  //fits_write_wrapper(cs->cwrap,cs->M,cs->N,"c3.fits");
  // }

  //printf("Done convolution\n");
  for (i=0; i<cs->M; i++)
    {
     for (j=0; j<cs->N; j++)
       {

	 
	 Z[k][i][j]=cs->cwrap[i][j];
	
       }
    }
  }

  //printf("Done making spectral cube. Fitting Gaussians!\n");

  // Z[k][i][j] is the oversampled convolved data cube with flux values along
  // the velocity x and y directions. 

  // Next, for each spatial bin along the velocity axis we fit a gaussian
  
  //i=25;
  //j=34;
  //printf("mean velocity: %f\n ",cs->vfield[i*oversample+oversample/2+cs->psf_xoffset][j*oversample+oversample/2+cs->psf_yoffset]);
  //printf("mean sigma: %f\n ",cs->sfield[i*oversample+oversample/2+cs->psf_xoffset][j*oversample+oversample/2+cs->psf_yoffset]);
  
  //for (m=0; m<nvel; m++){
  //  q=0.;
  //  for (k=0; k<oversample; k++){
  //    for (l=0; l<oversample; l++){	
  //  	q+=Z[m][i*oversample+k+cs->psf_xoffset][j*oversample+l+cs->psf_yoffset];
  //    }
  //  }
  //  cs->gdat.y[m]=q;
  // }
  //printdarr(cs->gdat.x,cs->gdat.n,"Xdata");
  //printdarr(cs->gdat.y,cs->gdat.n,"Ydata");
  //printdarr(cs->gdat.dy,cs->gdat.n,"dY data");
  
  //gaussfit(cs->gdat,cs->gpars,3);
  //printf("Done!\n");
  //printdarr(cs->gpars,3,"gpars");
  //exit(0);
  
  outflux=(double **)malloc(nxPix*sizeof(double *));
  for (i=0; i<nxPix; i++){
    outflux[i]=(double *)malloc(nyPix*sizeof(double));
  }
  

  //
  //FILE *f;	
  //f=fopen("outcube.txt","w");
  //fprintf(f,"# i j v f\n"); 
  for (i=0; i<nxPix; i++){
    for (j=0; j<nyPix; j++){
      //printf("i ,j= %i %i \n",i,j);
      for (m=0; m<nvel; m++){
	q=0.;
	for (k=0; k<oversample; k++){
	  for (l=0; l<oversample; l++){
	    q+=Z[m][i*oversample+k+cs->psf_xoffset][j*oversample+l+cs->psf_yoffset];
	  }
	}
	cs->gdat.y[m]=q;
	//fprintf(f,"%i %i %lf %lf\n",i,j,cs->gdat.x[m],cs->gdat.y[m]);
      }
      gaussfit(cs->gdat,cs->gpars,3);
      //printf("Gauss params: %f %f %f\n",cs->gpars[0],cs->gpars[1],cs->gpars[2]);
      outflux[i][j]=cs->gpars[0]*cs->gpars[2]*sqrt(2.0*M_PI);
      cs->outvel[i][j]=cs->gpars[1];
      cs->outsigma[i][j]=cs->gpars[2]*cs->gpars[2]-sigma_inst*sigma_inst;
      if (cs->outsigma[i][j] < 0.){
	cs->outsigma[i][j]=0.;
      } 
      cs->outsigma[i][j]=sqrt(cs->outsigma[i][j]);
      //cs->outsigma[i][j]=cs->gpars[2];
    }
  }
  //fclose(f);	
  #ifdef DEBUG
  fits_write_wrapper(outflux,nxPix,nyPix,"outflux.fits");
  #endif
  //clean up

  free(rad_pc);

  free(sigma_pc.sigma);
 
  free(pot_pc.sigma);
  
  free(sigma_pc.area);
 
  free(pot_pc.area);
  free(pot_pc.q);
 
  for (i=nxPix-1; i>=0; i--){
    free(outflux[i]);
  }
  free(outflux);
  

  //return with the new model.

  return;
}



void initialize_convolution(struct mdb_conv_struc *cs, int nvelo, double vsample){
  int i, j, M,N;
 
  M=cs->M;
  N=cs->N;
  cs->nvel=nvelo;
  cs->vsample=vsample;
  
  // First we initialize the pfs. We only use one block for this
  cs->psf = (double *) fftw_malloc(M*N * sizeof(double));
  cs->PSF = (fftw_complex*) fftw_malloc(M*(N/2+1) * sizeof(fftw_complex));
  cs->ppsf = fftw_plan_dft_r2c_2d(M, N, cs->psf, cs->PSF, FFTW_MEASURE);
    
  cs->psfwrap=(double **) malloc(M*sizeof(double *));




  cs->a=(double *)fftw_malloc(M*N * sizeof(double));
  cs->c=(double *)fftw_malloc(M*N * sizeof(double));
  cs->f=(double *)fftw_malloc(M*N * sizeof(double));

  cs->A=(fftw_complex *)fftw_malloc(M*(N/2+1) * sizeof(fftw_complex));
  cs->C=(fftw_complex *)fftw_malloc(M*(N/2+1) * sizeof(fftw_complex));

  cs->awrap=(double **)malloc(M* sizeof(double *));
  cs->cwrap=(double **)malloc(M* sizeof(double *));
  cs->fwrap=(double **)malloc(M* sizeof(double *));
  cs->vfield=(double **)malloc(M * sizeof(double *));
  cs->sfield=(double **)malloc(M * sizeof(double *));
  cs->vdata=(double *)malloc(M*N*sizeof(double));
  cs->sdata=(double *)malloc(M*N*sizeof(double));
  cs->Zdata=(double *)malloc(M*N*cs->nvel*sizeof(double));
    
  
  //printf("Size of Z: %i %i %i",M,N,cs->nvel);

  for (i=0; i<M; i++){
    cs->psfwrap[i]=&cs->psf[N*i+0];
    cs->awrap[i]=&cs->a[N*i+0];
    cs->cwrap[i]=&cs->c[N*i+0];
    cs->fwrap[i]=&cs->f[N*i+0];
    cs->vfield[i]=&cs->vdata[N*i+0];
    cs->sfield[i]=&cs->sdata[N*i+0];
    
  } // this allows us to write psfwrap[M][N] & cetera; 
 
  cs->Z=(double ***)malloc(cs->nvel*sizeof(double **));
  for (j=0; j<cs->nvel; j++){
    cs->Z[j]=(double **)malloc(M*sizeof(double *));
    for (i=0; i<M; i++){
      cs->Z[j][i]=&cs->Zdata[M*N*j + N*i+0];
    }
  }



  cs->pa=fftw_plan_dft_r2c_2d(M, N, cs->a, cs->A, FFTW_MEASURE);
  cs->pcinv=fftw_plan_dft_c2r_2d(M, N, cs->C, cs->c, FFTW_MEASURE);
  
  cs->gdat.dy=(double *)malloc(cs->nvel*sizeof(double));
  cs->gdat.y=(double *)malloc(cs->nvel*sizeof(double));
  cs->gdat.x=(double *)malloc(cs->nvel*sizeof(double));
  for (i=0; i<cs->nvel; i++){
    cs->gdat.x[i]=(-0.5*cs->nvel + i) * cs->vsample;
    cs->gdat.dy[i]=1.;

  }
  
  cs->gpars=(double *)malloc(3*sizeof(double));

  cs->gdat.n=cs->nvel;


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
	  cs->psfwrap[(i+M)%M][(j+N)%N]+=psf->area[k]/total_flux*exp(-0.5*pow( dist/(psf->sigma[k]/dx),2.0));
	  // maybe should change to gsl_sf_erf (although with the high oversampling
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
  M=cs->M; // This is the size of the oversampled image
  N=cs->N;

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
  //fits_write_wrapper(cs->fwrap,M,N,"f1.fits");


  // we do no actual fft, since it is done somewhere else
  return;
  

}




void clean_up_convolution(struct mdb_conv_struc *cs){
  // destroy your plans!
  int i;
  fftw_destroy_plan(cs->pa);
  fftw_destroy_plan(cs->ppsf);
  fftw_destroy_plan(cs->pcinv);


  fftw_free(cs->a);
  fftw_free(cs->psf);
  fftw_free(cs->c);
  fftw_free(cs->f);

  
  free(cs->awrap);
  free(cs->psfwrap);
  free(cs->cwrap);
  free(cs->fwrap);

  for (i=cs->nvel-1; i>0; i--){
    free(cs->Z[i]);
  }
  //  
  free(cs->Z);
  free(cs->Zdata);

  free(cs->sdata);
  free(cs->vdata);
  free(cs->vfield);
  free(cs->sfield);

  free(cs->gdat.x);
  free(cs->gdat.y);
  free(cs->gdat.dy);
  
  free(cs->outsigma);
  free(cs->outvel);
  
  fftw_free(cs->A);
  fftw_free(cs->C);
  fftw_free(cs->PSF);

  
  
  cs->initialized=False;
  
  
  return;
  
}



void mdb_fftw_convolve_2d(struct mdb_conv_struc *cs){
  int M,N;
  int i, j, ij;
  //double omega;
  fftw_complex *A, *B, *C;
  fftw_real scale;
  M=cs->M;
  N=cs->N;
  scale = 1.0 / (M * N);
  //omega=2.0*M_PI;
  A=cs->A;
  B=cs->PSF;
  C=cs->C;



  

  fftw_execute(cs->pa);


  for (i=0; i<M; i++){
    for (j=0; j<(N/2+1); j++){
      ij = i*(N/2+1) + j;
      C[ij][0] = (A[ij][0] * B[ij][0] - A[ij][1] * B[ij][1]) * scale;
      C[ij][1] = (A[ij][0] * B[ij][1] + A[ij][1] * B[ij][0]) * scale;
      

    }
  }
  fftw_execute(cs->pcinv);  // the smoothed surface density weighted velocity

  
 

  return;
}


void setup_conv_stuff(struct mdb_conv_struc *cs,int nx,int ny,struct multigaussexp *psf,double dx){
  int i, M,N;
  double maxsigma;
  

  cs->outsigma=(double **)malloc(nx*sizeof(double *));
  cs->outvel=(double **)malloc(nx*sizeof(double *));
  cs->outveldata=(double *)malloc(nx*ny*sizeof(double));
  cs->outsigmadata=(double *)malloc(nx*ny*sizeof(double));

  for (i=0; i<nx; i++){
    cs->outvel[i]=&cs->outveldata[ny*i+0];
    cs->outsigma[i]=&cs->outsigmadata[ny*i+0];
  }

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



  
// minradius should be in units of oversampled pixels, same for xc, yc
// PA should be in rad
  
void setup_flux_stuff_inputimage(struct image *fluxmap, struct mdb_conv_struc *cs,  double minradius, double xc, double yc, double PA, double q ){
  int i,j,k,M,N;
  double dist, xacc, yacc, minval;
  //double total_flux;
  int xorig, yorig;
  M=cs->M; // This is the size of the oversampled image
  N=cs->N;

  minval=0.;// not sure if this works!
  for (i=0; i<M; i++){
    for (j=0; j<N; j++){
      //printf("i ,j= %i %i \n",i,j);
      cs->fwrap[i][j]=minval;
      xorig=(int)(((double) i - cs->psf_xoffset)/((double) oversample));
      yorig=(int)(((double) j - cs->psf_yoffset)/((double) oversample)) ;
      if ((xorig< fluxmap->nx) && (xorig > 0) && (yorig > 0) && (yorig <fluxmap->ny))
	{
	  cs->fwrap[i][j]=fluxmap->pixelwrap[xorig][yorig];
	  xacc = (i-xc)*cos(PA) - (j-yc)*sin(PA);
	  yacc = (i-xc)*sin(PA) + (j-yc)*cos(PA);
	  
	  dist=sqrt(pow(xacc,2.)+pow(yacc,2.)/pow(q,2.0));


	  if (dist < minradius){
	    cs->fwrap[i][j]=0.;
	  }
	}
    }} 
			  

  
  //printf("Writing line model\n");
  //fits_write_wrapper(cs->fwrap,M,N,"f2.fits");


  // we do no actual fft, since it is done somewhere else
  return;
  

}


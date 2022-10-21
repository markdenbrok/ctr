/************************************************************************
 * This code was written by Mark den Brok (Utah 2015)                   *
 * to model NGC 404 for Dieu Nguyen                                     * 
 * denbrok@physics.utah.edu                                             *
 * Test model for the tilted ring stuff                                 *
 ************************************************************************/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"tr.h"
#include"readcol.h"
#include "memcee_structs.h"
#include "memcee.h"

#define oversample 2


struct n3862_ll_params {
  struct multigaussexp pot_orig, *pot_var; // the second one is an array for threads
 
  struct multigaussexp psf;
  struct multigaussexp sigma;
  struct image fluxmap;
  struct image velmap, velerr, dispmap, disperr;
  struct mdb_conv_struc *cs;
  double *rad;
  int nrad;
  int nvel;
  double plateScale;
  int nx;
  int ny;
  double dv;
};


double n3862_likelihood(double *arr,struct n3862_ll_params* p, int walker){
  int i, j,nx, ny, nborder;
  double *incl;
  double *pa;
  double mbh;
  double ml;
  double xc, yc;
  double galaxy_distance;
  double sigma_inst;
  double dx, dy, r_hole;
  double vSist, vsys;
  double chi2vel, chi2sig;
  incl=(double *)malloc(p->nrad*sizeof(double));
  pa=(double *)malloc(p->nrad*sizeof(double));

  nborder=22;
  nx=p->nx;
  ny=p->ny;
  xc=(24.0 + 0.0); //32.3, 28.43
  yc=(25.0 + 0.0); // this is the stellar photocenter

  galaxy_distance=10.; //Mpc
  sigma_inst=10.; // should be even lower
  vsys=0.;
  vSist=0.;
  printf("Thread: %i\n",walker);


  mbh=pow(10.,5.+3*arr[0]); // Black hole mass between 10^5 and 10^8
  ml=0.1+arr[1]*2.; // M/L between 0.1 and 2.1
  for (i=0; i<p->pot_orig.ntotal; i++){
    p->pot_var[walker].area[i]=p->pot_orig.area[i]*ml;
  }

  for (i=0; i<p->nrad; i++){
    incl[i]=(74.+arr[2]*16.)*M_PI/180.;
    pa[i]=(175.51 + (arr[3]*30.-15) - 90.)*M_PI/180.;
  }

  r_hole=arr[4]*0.08;
  dx=arr[5]*12.0-6.0;
  dy=arr[6]*12.0-6.0;
  
 
  line_modeling(mbh, &p->pot_var[walker], &p->sigma,p->rad, pa, p->nrad, xc+dx, yc+dy, p->plateScale, nx,ny,vSist,sigma_inst, &p->psf, galaxy_distance,incl, &p->cs[walker], NULL, p->nvel, p->dv, &p->fluxmap,incl[0],r_hole,pa[0],cos(incl[0]));

  vsys=arr[7]*100.-50.;
  

  chi2vel=0.;
  chi2sig=0.;
  for (i=nborder; i<nx-nborder; i++){
    for (j=nborder; j<ny-nborder; j++){	  
      if (p->velmap.pixelwrap[i][j] > -500.){ // only look at valid pixels
      chi2vel+=pow( (p->velmap.pixelwrap[i][j]-p->cs[walker].outvel[i][j]-vsys)/p->velerr.pixelwrap[i][j] , 2.);
      chi2sig+=pow( (p->cs[walker].outsigma[i][j]-p->dispmap.pixelwrap[i][j])/p->disperr.pixelwrap[i][j] , 2.);
      }
    } //j 
  } //i
  free(incl);
  free(pa);
  return -0.5*chi2vel;
}










int main(){
  
  int i, j, nthread;
  double *rad, *pa, *errpa, *q, *errq;
  int *nfit;
  
  struct n3862_ll_params llparams;
  struct multigaussexp mge_co;
  struct McParams mpa={200000,100,"mcmc/out_all_ngc3862.txt","mcmc/out_select_ngc3862.txt",0};

  char *fmt="%lf %lf %lf %lf %lf %i";
  char *rcflags=" bufsize=1024 filelen=39 skipsym=# ";
  char *filename="input/kinemetry_output.txt";
  
  nthread=18;

  //  llparams.fluxmap=read_fits("input/co_fluxmap.fits");
  llparams.fluxmap=read_fits("input/model_flux.fits");
  llparams.velmap=read_fits("input/co_velmap.fits");
  llparams.dispmap=read_fits("input/co_dispmap.fits");
  llparams.velerr=read_fits("input/co_velerr.fits");
  llparams.disperr=read_fits("input/co_disperr.fits");
  llparams.plateScale=0.059;
  llparams.nvel=100;
  llparams.dv=10.;
  llparams.nx=101;
  llparams.ny=101;

  llparams.pot_var=(struct multigaussexp *)malloc(nthread*sizeof(struct multigaussexp));
  printf("Reading potential MGE\n");
  read_mge("input/ngc3862_mge.txt",20, &llparams.pot_orig); 
  printf("Reading CO mge\n");
  read_mge("input/CO_ellip.txt",20, &mge_co); 
  printf("Read CO mge\n");
  for (i=0; i<nthread; i++){
    read_mge("input/ngc3862_mge.txt",20, &llparams.pot_var[i]);
    for (j=0; j<mge_co.ntotal; j++){
      llparams.pot_var[i].area[ llparams.pot_var[i].ntotal+j]=mge_co.area[j];
      llparams.pot_var[i].sigma[ llparams.pot_var[i].ntotal+j]=mge_co.sigma[j];
      llparams.pot_var[i].q[ llparams.pot_var[i].ntotal+j]=mge_co.q[j];
    }
    llparams.pot_var[i].ntotal+=mge_co.ntotal;
  }
  printf("Reading psf MGE\n");
  read_mge("input/psf.txt",20, &llparams.psf); // PSF MGE
  printf("Reading sigma MGE\n");
  read_mge("input/psf.txt",20, &llparams.sigma); // temporary until I make one 
 
  llparams.rad=(double *)malloc(39*sizeof(double));
  pa=(double *)malloc(39*sizeof(double));
  errpa=(double *)malloc(39*sizeof(double));
  q=(double *)malloc(39*sizeof(double));
  errq=(double *)malloc(39*sizeof(double));
  nfit=(int *)malloc(39*sizeof(int));
  
  llparams.nrad=readcol(filename,rcflags,fmt,llparams.rad,pa,q,errpa,errq,nfit);
  printf("Read the kinemetry file\n");

  llparams.cs=(struct  mdb_conv_struc *)malloc(nthread*sizeof(struct mdb_conv_struc));
  for (i=0; i<nthread; i++){
    llparams.cs[i].initialized=False;
  }


  printf("Initializing convolution structure in threadsafe way\n");
  for (i=0; i<nthread; i++){
    printf("Thread: %i \n",i);
    
    //n3862_likelihood(arr, &llparams, i);
    setup_conv_stuff(&llparams.cs[i],llparams.nx*oversample,llparams.ny*oversample, &llparams.psf,llparams.plateScale/oversample);
  
    initialize_convolution(&llparams.cs[i],llparams.nvel,llparams.dv);
  
    setup_psf_stuff(&llparams.psf, &llparams.cs[i], llparams.plateScale/oversample);
    
  }

  printf("Starting MCMC!\n");
  memcee(mpa, 8, &n3862_likelihood , &llparams,nthread);
 
 
  printf("Cleaning up!\n");
 
  for (i=0; i< nthread; i++){
    clean_up_convolution(&llparams.cs[i]);
  }
  free(llparams.cs);
  free_mge(&llparams.psf);
  free_mge(&llparams.pot_orig);
 
  free_mge(&llparams.sigma);

  free(llparams.rad);
  free(pa);
  free(errpa);
  free(q);
  free(errq);
  

  //  fits_write_wrapper(image,nx,ny,"out.fits");
  
  return 0;
}



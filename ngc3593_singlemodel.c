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


int main(){
  int i,j;
  int nx, ny;
  double xc,yc;
  double plateScale, mbh;
  double ml;
  struct multigaussexp psf, sigma, pot, potml, co;
  double **image;
  double vSist,vsys,sigma_inst;
  double galaxy_distance;
  struct mdb_conv_struc cs;
  double *rad, *pa, *errpa, *q, *errq, *pa_varia, *incl;
  char *fmt="%lf %lf %lf %lf %lf %i";
  char *rcflags=" bufsize=1024 filelen=50 skipsym=# ";
  struct image fluxmap;
  struct image velmap, velerr, dispmap, disperr;
  char *filename="input/kinemetry_output.txt";
  int nrad, nvalid;
  int *nfit;
  double mass;
  double chi2vel, chi2sig;
  int iinc, nborder;
  rad=(double *)malloc(50*sizeof(double));
  pa=(double *)malloc(50*sizeof(double));
  errpa=(double *)malloc(50*sizeof(double));
  q=(double *)malloc(50*sizeof(double));
  errq=(double *)malloc(50*sizeof(double));
  nfit=(int *)malloc(50*sizeof(int));
  pa_varia=(double *)malloc(50*sizeof(double));
  incl=(double *)malloc(50*sizeof(double));
  //fluxmap=read_fits("input/h2_intmap_real.fits");
  fluxmap=read_fits("input/co_fluxmap.fits");
  velmap=read_fits("input/co_velmap.fits");
  dispmap=read_fits("input/co_dispmap.fits");
  velerr=read_fits("input/co_velerr.fits");
  disperr=read_fits("input/co_disperr.fits");
  


  nrad=readcol(filename,rcflags,fmt,rad,pa,q,errpa,errq,nfit);
  
  printf("Read the kinemetry file\n");
  
  cs.initialized=False;
  nx=301;
  ny=301;
  xc=(150.0 + 0.) ; //32.3, 28.43
  yc=(150.0 + 0.0); // this is the stellar photocenter
  plateScale=0.059;
  mbh=pow(10.,6.93);//(4.e5);//1.0e4;//0.0;//5.0e5;
  ml=0.62;//pow(10.,0.006*0.7);//1.250;
  vSist=0.0;
  sigma_inst=20.;
  galaxy_distance=10.; //Mpc

  image=(double **)malloc(nx*sizeof(double *));
  for (i=0; i<nx; i++) image[i]=(double *)malloc(ny*sizeof(double));


  read_mge("input/ngc3593_mge.txt",20, &pot); //in M_sun/pc^2,arcmin
  read_mge("input/ngc3593_mge.txt",20, &potml); 
  //read_mge("input/h2_int_mge_nifs.txt",20, &disc);
  //disc.pa=(double *)malloc(20*sizeof(double));
  //for (i=0; i<20; i++) disc.pa[i]=(30.)*M_PI/180.; //According toDieu, but not checked 
  read_mge("input/CO_ellip.txt",20,&co);
  read_mge("input/psf.txt",20, &psf); // PSF Mge
  read_mge("input/psf.txt",20, &sigma); // temporary until I make one 



  mass=0.;
  for (i=0; i<pot.ntotal; i++){
    mass += 2.*M_PI*pot.area[i]*pot.q[i]*pow(pot.sigma[i]*galaxy_distance*(1.E6)*M_PI/180./3600.,2.);
  }
  printf("Luminosity [x10^6 Lo] in MGE expansion, %lf \n",mass/(1.E6));
  
  
  for (i=0; i<pot.ntotal; i++)
    { pot.area[i]=potml.area[i]*ml;}
  for (i=pot.ntotal; i<pot.ntotal+co.ntotal; i++)
    { pot.area[i]=co.area[i-pot.ntotal];
      pot.sigma[i]=co.sigma[i-pot.ntotal];
      pot.q[i]=co.q[i-pot.ntotal];
      
    }
  pot.ntotal+=co.ntotal;
  iinc=0;
  for (i=0; i<nrad; i++){
    pa_varia[i]=(175.51 +5. -90.)*M_PI/180.;
    incl[i]=(73.8)*M_PI/180.; // instead of 62!
    //incl[i]=(101.)*M_PI/180.; // instead of 62!
  }
  line_modeling(mbh, &pot, &sigma,rad, pa_varia, nrad, xc, yc, plateScale, nx,ny,vSist,sigma_inst, &psf, galaxy_distance,incl, &cs, image, 100, 10.0, &fluxmap,incl[0], 0.08*0.0507, pa_varia[0],cos(incl[i]));
  chi2vel=0.;
  chi2sig=0.;
  nborder=19;
  vsys=12.25;
  nvalid=0;
  for (i=nborder; i<nx-nborder; i++){
    for (j=nborder; j<ny-nborder; j++){
      if (velmap.pixelwrap[i][j] > -500.){ 
	nvalid++;
	chi2vel+=pow( (cs.outvel[i][j]+vsys-velmap.pixelwrap[i][j])/velerr.pixelwrap[i][j] , 2.);

     
	chi2sig+=pow( (cs.outsigma[i][j]-dispmap.pixelwrap[i][j])/disperr.pixelwrap[i][j] , 2.);
      }
      
    }
  }
  printf("%lf %lf %lf %lf %lf  %lf %i\n",incl[0],ml, mbh,-0.5*chi2vel,-0.5*chi2sig,chi2vel/nvalid,nvalid);
  fits_write_wrapper(cs.outvel,nx,ny,"single_model_velocity_newmge.fits");
  fits_write_wrapper(cs.outsigma,nx,ny,"single_model_sigma_newmge.fits");


  for (i=0; i<nx; i++) {
    free(image[i]);}
  free(image);
  clean_up_convolution(&cs);

  free_mge(&psf);
  free_mge(&pot);
  free_mge(&potml);
  free_mge(&sigma);

  free(rad);
  free(pa);
  free(errpa);
  free(q);
  free(errq);
 
  free(pa_varia);
  free(incl);
  

  //  fits_write_wrapper(image,nx,ny,"out.fits");
  
  return 0;
}



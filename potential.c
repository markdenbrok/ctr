
/************************************************************************
 * This code was written by Mark den Brok (Utah 2014)                   *
 * based on the IDL code of Cappellari                                  *
 * denbrok@physics.utah.edu                                             *
 ************************************************************************/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<gsl/gsl_sf_erf.h>
#include <gsl/gsl_integration.h>

#include "tr.h"


// The pot mge should be spherical and with sigma in pc, dens = mass of gauss in solar mass




double qpfunc(double x, void *params){
  double result;
  //convert params back to structs
  struct funcparams* fparams;
 
  struct multigaussexp potential;
  double R; 

  fparams=params;
  potential=fparams->potential;
  R=fparams->R;
  
  result = potential_integrand(x, &potential, R);
  return result;
  
}

double qpint1d(double a, double b, double epsrel, double epsabs,struct funcparams *fparams){
  gsl_integration_workspace *w  = gsl_integration_workspace_alloc (1000);
  double result, error;
  gsl_function F;
  F.function = &qpfunc;
  F.params = fparams;

  epsrel=1.E-5;

  gsl_integration_qag (&F, a, b, epsabs, epsrel, 1000 ,2,w,&result, &error); 

  gsl_integration_workspace_free (w);
  
  return result;
}

double potential_integrand(double u, struct multigaussexp *pot, double r){
  int i;
  double  result;
  result=0.;
  for (i=0; i<pot->ntotal; i++){
    result+=pot->area[i]/pow(pot->sigma[i],3.)*r*(u*u)*exp(-0.5*pow(u*r/pot->sigma[i],2.))/sqrt(1.-(1-pot->q[i]*pot->q[i])*u*u);
  }
  return result;
}



double *dphi(double *r, struct multigaussexp *pot, double mbh, int n){
  double *dphis;
  struct funcparams fparams;
 
  fparams.potential.area=pot->area;
  fparams.potential.sigma=pot->sigma;
  fparams.potential.q=pot->q;
  fparams.potential.ntotal=pot->ntotal;
  int i,j;
  dphis=(double *)malloc(n*sizeof(double));
  for (i=0; i<n; i++){
    dphis[i]=G*mbh/pow(r[i],2.0);
    fparams.R=r[i];
    dphis[i]+=G*sqrt(2./M_PI)*qpint1d(0.,1.,1E-5,0.0,&fparams);
    
    //for (j=0; j<pot->ntotal; j++){
    //  dphis[i]+=-1.0*(G*pot->area[j]/r[i])*
    //	(exp(-0.5*pow(r[i]/pot->sigma[j],2.)) * sqrt(2.0/M_PI)/pot->sigma[j] /
    //	 - gsl_sf_erf(r[i]/(sqrt(2.0)*pot->sigma[j]))/r[i]);
    //}
  }
  
  // here the hot disc can be implemented.


  

  return dphis;

}



// this is for a hot disk
//everything in pc and km/s
double *dsigma2drhorho (struct multigaussexp *mge_rho,struct multigaussexp *mge_sig, double *r, int n){
  // sig = sum_i sig_i * exp(-0.5*(r/s)^2)
  // dsig = sum_i sig_i * exp(-0.5*(r/s)^2) * r/(s^2)
  // disg2 = sig*disg*2.0
  int i, j;
  double *sigmar, *dsigmar;
  double *rhor, *drhor;
  sigmar= (double *)malloc(n*sizeof(double));
  dsigmar= (double *)malloc(n*sizeof(double));

  rhor= (double *)malloc(n*sizeof(double));
  drhor= (double *)malloc(n*sizeof(double));

  for (j=0; j<n; j++) {
    sigmar[j]=0.;
    dsigmar[j]=0.;
    rhor[j]=0.;
    drhor[j]=0.;
    for(i=0; i<mge_sig->ntotal; i++){
      sigmar[j]+=mge_sig->area[i]*exp(-0.5*pow(r[j]/mge_sig->sigma[i],2.0));
      dsigmar[j]+=mge_sig->area[i]*exp(-0.5*pow(r[j]/mge_sig->sigma[i],2.0))*r[j]/(mge_sig->sigma[i]*mge_sig->sigma[i]);
      rhor[j]+=mge_rho->area[i]*exp(-0.5*pow(r[j]/mge_rho->sigma[i],2.0));
      drhor[j]+=mge_rho->area[i]*exp(-0.5*pow(r[j]/mge_rho->sigma[i],2.0))*r[j]/(mge_rho->sigma[i]*mge_rho->sigma[i]);

    }
    sigmar[j]=sigmar[j]*dsigmar[j]*2.0;
  }
  free(dsigmar);




  return sigmar;
}


// return the squared azimuthal velocity. This is calculated as the local gradient in the potential (force)
// it is possible to put v_azi and dphi in one function, but this is clearer
// r should be in pc, pot->dens=mass in Msun, sig=[pc],
double *v_azi2(double *r, struct multigaussexp *pot, double mbh, int n){
  double *vazi;
  int i;


  vazi=dphi(r,pot,mbh,n); // cold disc


  for (i=0; i<n; i++){
    vazi[i]*=r[i];
  }

  // For the hot disk, see above
  
  return vazi;
 
}

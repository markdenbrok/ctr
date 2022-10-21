#include<stdio.h>
#include<math.h>
#include <gsl/gsl_integration.h>

struct funcparams {
  struct multigaussexp potential;
  double R;
};

double qpfunc(double, void *);

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


double qpfunc(double x, void *params){
  double result;
  //convert params back to structs
  struct funcparams* fparams;
  struct multigaussexp luminosity;
  struct multigaussexp potential;
  double X1; double Y1; 
  double inc; 
  double *beta;

  fparams=params;
  luminosity=fparams->luminosity;
  potential=fparams->potential;
  X1=fparams->X1;
  Y1=fparams->Y1;
  inc=fparams->inc;
  beta=fparams->beta;
  //printf("Doing integral %i\n",fparams->which_integral);
 
    result = janis2_jeans_mge_integrand(x, &luminosity, &potential, X1, Y1, inc, beta);
    return result;
  
}

  


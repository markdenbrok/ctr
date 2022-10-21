#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>

struct gauss_data{
  int n;
  double *x;
  double *y;
  double *dy;
};


int gauss_f (const gsl_vector *, void *, gsl_vector *);

int gauss_df (const gsl_vector *, struct gauss_data *,  gsl_matrix *);


int gauss_fdf (const gsl_vector *, void *,gsl_vector *, gsl_matrix *);


int gaussfit(struct gauss_data , double *, int );


void print_state (int , gsl_multifit_fdfsolver * );

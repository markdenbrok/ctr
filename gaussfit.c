#include "gaussfit.h"
#include "stuff.h" 


// this is based on the GSL example to fit an exponential+background
// except lambda->mu , b->sigma




int gauss_f (const gsl_vector *par, void *data, gsl_vector *f)
{
  int n = ((struct gauss_data *)data)->n;
  double *x = ((struct gauss_data *)data)->x;
  double *y = ((struct gauss_data *)data)->y;
  double *dy = ((struct gauss_data *) data)->dy;

  double A = gsl_vector_get (par, 0);
  double mu = gsl_vector_get (par, 1);
  double sigma = gsl_vector_get (par, 2);

  double Yi;
  
  size_t i;

  for (i = 0; i < n; i++)
    {
      /* Model Yi = A * exp(-lambda * i) + b */
      
      Yi = A * exp( -0.5*pow((x[i]-mu)/sigma,2.));
      gsl_vector_set (f, i, (Yi - y[i])/dy[i]);
    }

  return GSL_SUCCESS;
}


int gauss_df (const gsl_vector * par, struct gauss_data *data, 
         gsl_matrix * J)
{
  int n = data->n;
  double *x = data->x;
  double *dy = data->dy;

  double A = gsl_vector_get (par, 0);
  double mu = gsl_vector_get (par, 1);
  double sigma = gsl_vector_get (par, 2);
  double ei;

  // 


  size_t i;

  for (i = 0; i < n; i++)
    {
      /* Jacobian matrix J(i,j) = dfi / dxj, */
      /* where fi = (Yi - yi)/dy[i],      */
      /*       Yi = A * exp(-0.5 * ((x-mu)/sigma)^2)  */
      /* and the parj are the parameters (A,mu,sigma) */
      /* dfi/dA =  exp(-0.5 * ((x-mu)/sigma)^2)   */
      /* dfi/dmu = A exp(-0.5 * ((x-mu)/sigma)^2) * ((x-mu)/(sigma^2)) */
      /* dfi/dsigma =  A exp(-0.5 * ((x-mu)/sigma)^2)* (x-mu)/sigma * (x-mu)/sigma^2  */
      
	
      ei = exp(-0.5* pow( (x[i]-mu)/sigma,2.));      


      gsl_matrix_set (J, i, 0, ei/dy[i]); 
      gsl_matrix_set (J, i, 1, A*ei*(x[i]-mu)/(sigma*sigma)/dy[i]);
      gsl_matrix_set (J, i, 2, A*ei*pow((x[i]-mu)/sigma,2.)/sigma/dy[i]);
    }
  return GSL_SUCCESS;
}


int gauss_fdf (const gsl_vector * par, void *data,gsl_vector * f, gsl_matrix * J)
{
  gauss_f (par, data, f);
  gauss_df (par, data, J);
  return GSL_SUCCESS;
}




int gaussfit(struct gauss_data data, double *params, int nparams){

  const gsl_multifit_fdfsolver_type *T;
  gsl_multifit_fdfsolver *s;
  int status;
  unsigned int i, j, iter = 0;
  int n, imax;
  const size_t p = 3;
  const gsl_rng_type *type;
  gsl_rng * r;
  double par_init[3];
  //double *y, *dy;
  gsl_matrix *covar = gsl_matrix_alloc (p, p);
  gsl_matrix *J;
  
  gsl_multifit_function_fdf f;

  for (j=0; j<3; j++) {par_init[j]=params[j];}
  
  gsl_vector_view par = gsl_vector_view_array (par_init, p);
  
  n=data.n;
  //y=data.y;
  //dy=data.dy;

  imax=max_in_array_index(data.y, data.n);
  par_init[0]=data.y[imax];
  par_init[2]=20.;
  par_init[1]=data.x[imax];

  //printdarr(data.x,data.n,"xdata");
  //printdarr(data.y,data.n,"ydata");
  //printdarr(data.dy,data.n,"dydata");

  //printf("n: %i p: %i \n",n, p);

  gsl_rng_env_setup();

  type = gsl_rng_default;
  r = gsl_rng_alloc (type);

  f.f = &gauss_f;
  f.df = &gauss_df;
  f.fdf = &gauss_fdf;
  f.n = n;
  f.p = p;
  f.params = &data;

  /* This is the data to be fitted */

  T = gsl_multifit_fdfsolver_lmsder;
  s = gsl_multifit_fdfsolver_alloc (T, n, p);
  gsl_multifit_fdfsolver_set (s, &f, &par.vector);

  //print_state (iter, s);

  do
    {
      iter++;
      status = gsl_multifit_fdfsolver_iterate (s);

      //printf ("status = %s\n", gsl_strerror (status));

      //print_state (iter, s);

      if (status)
        break;

      status = gsl_multifit_test_delta (s->dx, s->x,
                                        1e-4, 1e-4);
    }
  while (status == GSL_CONTINUE && iter < 500);

  // MdB: 20160930 replacing this:  gsl_multifit_covar (s->J, 0.0, covar);
  J = gsl_matrix_alloc(s->fdf->n,s->fdf->p);
  gsl_multifit_fdfsolver_jac(s,J);
  gsl_multifit_covar(J,0.0,covar);
    

  
#define FIT(i) gsl_vector_get(s->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))

  //{ 
  //  double chi = gsl_blas_dnrm2(s->f);
  //  double dof = n - p;
  //  double c = GSL_MAX_DBL(1, chi / sqrt(dof)); 

    //printf("chisq/dof = %g\n",  pow(chi, 2.0) / dof);

    //printf ("A      = %.5f +/- %.5f\n", FIT(0), c*ERR(0));
    //printf ("mu     = %.5f +/- %.5f\n", FIT(1), c*ERR(1));
    //printf ("sigma  = %.5f +/- %.5f\n", FIT(2), c*ERR(2));
  //}

  //double chi = gsl_blas_dnrm2(s->f);
  //c=GSL_MAX_DBL(1, chi / sqrt(dof));
  for (i=0; i<3; i++){
    params[i]=FIT(i);
    //dparams[i]=c*ERR(i)
  }

  //printf ("status = %s\n", gsl_strerror (status));

  gsl_multifit_fdfsolver_free (s);
  gsl_matrix_free (covar);
  gsl_matrix_free (J);
  gsl_rng_free (r);
  return 0;
}

void print_state (int iter, gsl_multifit_fdfsolver * s)
{
  printf ("iter: %3i par = % 15.8f % 15.8f % 15.8f "
          "|f(x)| = %g\n",
          iter,
          gsl_vector_get (s->x, 0), 
          gsl_vector_get (s->x, 1),
          gsl_vector_get (s->x, 2), 
          gsl_blas_dnrm2 (s->f));

}

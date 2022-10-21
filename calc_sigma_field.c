#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"stuff.h"
#include "gsl/gsl_errno.h"
#include "gsl/gsl_spline.h"
#include"tr.h"


//sigma in pc, dx in pix/pc

void calc_sigma_field(struct multigaussexp *sigma,double xc_os,double yc_os,double dx, int M, int N,double **sigma_field, double sig_inst){
  int i,j, k;
  double xacc, yacc,xc, yc;
  double dist;
  xc=xc_os;
  yc=yc_os;
  for (i=0; i<M; i++){
    for (j=0; j<N; j++){
      //sigma_field[i][j]=2.;
      for (k=0; k<sigma->ntotal; k++)
	{
	  //printf("PA: %f\n",disc->pa[k]);
	  // first define xacc, yacc
	  // xacc = (i-xc) cos(pa) - (j-yc) sin(pa)
	  // yacc = (i-xc) sin(pa) + (j-yc) cos(pa)
	  // pa should be radians, defined in the mathematical way (x through y)
	  
	  //xacc = (i-xc)*cos(sigma->pa[k]) - (j-yc)*sin(sigma->pa[k]);
	  //yacc = (i-xc)*sin(sigma->pa[k]) + (j-yc)*cos(sigma->pa[k]);
	  
	  xacc = (i-xc);
	  yacc = (j-yc); // assuming spherical symmetry throughout

	  dist=sqrt(pow(xacc,2.)+pow(yacc,2.));
	  sigma_field[i][j]+=sigma->area[k]*exp(-0.5*pow( dist/(sigma->sigma[k]/dx),2.0)); 
	}
      sigma_field[i][j]=20.;
      sigma_field[i][j]=sqrt(pow(sigma_field[i][j],2.) + sig_inst*sig_inst);// add the instrumental dispersion sig_inst to the sigma field by squaring/sqrt
     
      
    }
  }
  


}

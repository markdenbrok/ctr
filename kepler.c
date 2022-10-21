/************************************************************************
 * This code was written by Mark den Brok (Utah 2014)                   *
 * based on the IDL code of Cappellari/Neumayer                         *
 * however, the coordinate transformation's taken care of in another way*
 * denbrok@physics.utah.edu                                             *
 ************************************************************************/



#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"stuff.h"
#include "gsl/gsl_errno.h"
#include "gsl/gsl_spline.h"
#include"tr.h"


//galaxy distance is in pc
//inclination is in radians? maybe better in degrees
// rad in pc
//Mbh in solar mass
void keplerian_disc(double Mbh, struct multigaussexp *pot, double *rad, int nrad,double *pa_varia,double xc,double yc, double highres_platescale,int nx,int ny,double vSist,double sigma_strum, double galaxy_distance,double *incl, double **image){
  int i,j,k;
  double dx_pc;
  double r, v2, medinc, *rad_pc, medpa;
  double *vr1_tot;
  gsl_spline *rspline;
  gsl_interp_accel *acc_r;
  gsl_spline *incspline;
  gsl_interp_accel *acc_inc;
  gsl_spline *paspline;
  gsl_interp_accel *acc_pa;
  double inc_pref, pa_pref;
  double x1,x2,y1,y2;
  double *incl_copy, *pa_varia_copy;

  // We will turn all distances into pc

  dx_pc=highres_platescale*galaxy_distance*M_PI/180.0/3600.0; // pc per high-res model pixel
  
  pa_varia_copy=(double *)malloc(nrad*sizeof(double));
  incl_copy=(double *)malloc(nrad*sizeof(double));
  
  for (i=0; i<nrad; i++){
    incl_copy[i]=incl[i];
    pa_varia_copy[i]=pa_varia[i];
  }

  medinc=median2(incl_copy,nrad);
  medpa=median2(pa_varia_copy,nrad);
  //printf("median pa: %lf\n",medpa);
  // rad should be in pc
  
  free(incl_copy);
  free(pa_varia_copy);

  rad_pc=(double *)malloc(nrad*sizeof(double));
  //printf("radii: ");
  for (i=0;i<nrad; i++){
    rad_pc[i]=rad[i];
    //printf("%f ",rad_pc[i]);
  }
  //printf("\n");
  
  vr1_tot = v_azi2(rad_pc, pot,Mbh,  nrad);
  // this for testing:
  //print_mge(pot);
  //printf("Mbh: %lf\n",Mbh);
  //printdarr(rad_pc,nrad,"rad");
  //printdarr(vr1_tot,nrad,"vr1_tot");

  
  //printf("Here!\n");
  acc_r = gsl_interp_accel_alloc();
  rspline = gsl_spline_alloc(gsl_interp_cspline,nrad);
  gsl_spline_init(rspline,rad_pc,vr1_tot,nrad);


  //printf("Here!!\n");
  acc_pa = gsl_interp_accel_alloc();
  paspline = gsl_spline_alloc(gsl_interp_cspline,nrad);
  gsl_spline_init(paspline,rad_pc,pa_varia,nrad);

  //printf("Here!!!\n");
  acc_inc = gsl_interp_accel_alloc();
  incspline = gsl_spline_alloc(gsl_interp_cspline,nrad);
  gsl_spline_init(incspline,rad_pc,incl,nrad);
  
  //printf("Here!!!!\n");
  // next we have to interpolate over this.

  
  
  
  for (i=0; i<nx; i++){
    for (j=0; j<ny; j++){
      pa_pref=medpa;
      inc_pref=medinc;
      x1=(i-xc); // this is the pixel position before rotation
      y1=(j-yc);
      x2=-1.0*x1*sin(medpa)+1.0*y1*cos(medpa);
      y2=-1.0*x1*cos(medpa)-y1*sin(medpa); // this makes it consistent with kinemetry PAs
      
      r=sqrt(x2*x2 + y2*y2/(cos(inc_pref)*cos(inc_pref)))*dx_pc; // this is a quick fix.
      for(k=0; k<3; k++){
	if (r>rad[nrad-1]){
	  v2=vr1_tot[nrad-1]*rad_pc[nrad-1]/r; // extrapolate keplerain prof.
	  pa_pref=pa_varia[nrad-1];
	  inc_pref=incl[nrad-1];
	 
	}
	else{
	  if (r<rad[0])
	    {//printf("r smaller than rad[0]: %f\n",r);
	      v2=G*Mbh/r;
	    }
	  else {
	    pa_pref=gsl_spline_eval(paspline,r,acc_pa);
	    inc_pref=gsl_spline_eval(incspline,r,acc_inc);
	    v2=gsl_spline_eval(rspline,r,acc_r); 
	  }
	}
	x2=-1.0*x1*sin(pa_pref)+1.0*y1*cos(pa_pref);
	y2=-1.0*x1*cos(pa_pref)-1.0*y1*sin(pa_pref);
	r=sqrt(x2*x2 + y2*y2/(cos(inc_pref)*cos(inc_pref)))*dx_pc;
      }
      //printf("%lf %lf\n",r,v2);
      if (v2 < 0)
	v2=0.; // to avoid NANs
      if (isnan(v2))
	v2=0.;
      if (isinf(v2))
	v2=0.;
      if (r==0.)
	r=0.001*dx_pc; // MdB fix to avoid Nan
      image[i][j]=sqrt(v2)*sin(inc_pref)*(x2)*dx_pc/r + vSist;
    }
  }
   
#ifdef DEBUG   
fits_write_wrapper(image,nx,ny,"out_kp.fits");
#endif

  gsl_spline_free(rspline);
  gsl_interp_accel_free(acc_r);
  gsl_spline_free(paspline);
  gsl_interp_accel_free(acc_pa);
  gsl_spline_free(incspline);
  gsl_interp_accel_free(acc_inc);
  free(vr1_tot);
  free(rad_pc);
  return;


}

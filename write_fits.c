/************************************************************************
 * This code was written by Mark den Brok (Utah 2012)                   *
 * based on early code written by MdB (Groningen 2010)                  *
 * denbrok@physics.utah.edu                                             *
 ************************************************************************/




#include<stdio.h>
#include<stdlib.h>
#include<fitsio.h>
#include"tr.h"

// writes a simple image to a fitsfile. 
// mainly for diagnosis purposes

int fits_write_wrapper(double **im,int nx, int ny,char *filename){
  int i,j;
  struct image m;
  m.nx=nx;
  m.ny=ny;
  m.pixels=(double *)malloc(m.nx*m.ny*sizeof(double));
  m.pixelwrap=(double **)malloc(m.nx*sizeof(double *));
  for (i=0; i<m.nx; i++){
    m.pixelwrap[i]=&m.pixels[i*m.ny];
  }
  for (i=0; i<m.nx; i++){
    for (j=0; j<m.ny; j++){
      m.pixelwrap[i][j]=im[i][j];
    }
  }
  // MdB 20180829 fixed free here
  
  i= write_fits(filename, &m);
  free(m.pixels);
  free(m.pixelwrap);
  return i;

}


int write_fits(char *filename,struct image *im){
  fitsfile *fptr;
  int status =0;
  long naxes[2];
  long fpixel[2]={1l,1l};
  int nelements=im->nx*im->ny;
  int i,j;
  float *pixels;
  //char *filename="out.fits";
  FILE *out;

  pixels=(float *)malloc(im->nx*im->ny*sizeof(float));

  for (i=0; i<im->nx; i++){
    for (j=0; j<im->ny; j++) {
      pixels[j*im->nx + i] =(float) im->pixelwrap[i][j];
    } // for j
  }//for i


  naxes[0]=im->nx;
  naxes[1]=im->ny;
  status=0;
  // If the output file exists remove it
  if ((out = fopen (filename, "r")) != (FILE *) 0)
    {
      remove(filename);
    }
  

  fits_create_file(&fptr,filename, &status);
  status=0;

  fits_create_img(fptr,FLOAT_IMG,2,naxes,&status);
  status=0;
  fits_write_comment(fptr, "MdB writefits", &status);
  status=0;
  fits_write_pix(fptr,TFLOAT,fpixel,nelements,pixels,&status);
  status=0;
  fits_close_file(fptr,&status);
  
  free(pixels);
  return 0;

}//write_image

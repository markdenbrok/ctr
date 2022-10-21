#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<fitsio.h>
#include"tr.h"

/**********************************************

This reads an image into a struct image
form of image is:
image.pixels[y][x]
because this is the way that fits is stored.

 *********************************************/



struct image read_fits(char *filename){
  struct image im;
  fitsfile *fptr;
  int status, nf, anynull, i,j;
  long naxes[3];				    
  long fpixel,ntot;
  float *imasrow;
  float nullval;

  anynull=0;
  fpixel=1;
  nullval=0.;
  status=0;
#ifdef DEBUG
  printf("Reading file %s\n",filename);
#endif
  if (fits_open_file(&fptr, filename, READONLY, &status)) 
    {
      printf("Cannot open file: %s %i\n",filename,status);
        exit (1);
    } 
  if (fits_read_keys_lng(fptr, "NAXIS", 1, 2, &naxes[1], &nf, &status)) 
    {
      printf("Failed reading fits keys in %s, errorcode: %i\n",filename,status );
      exit(1);
    }
  //im=(struct image*)malloc(sizeof(struct image));
  im.nx=naxes[1];
  im.ny=naxes[2];
  ntot=(long)im.nx*im.ny;
  //we read the full image as one single row, later we will set pointers to this
  imasrow=(float *)malloc(im.nx*im.ny*sizeof(float));
  if ( fits_read_img(fptr, TFLOAT, fpixel, ntot, &nullval,imasrow, &anynull, &status) ) 
    {
      free(imasrow);
      printf("Failed reading img %s, errorcode: %i\n",filename,status );  
      exit(1);
    }

  // this is important: the x-coordinate is the fast coordinate
  // which we're going to change, since we're converting to doubles
  // anyhow.
  im.pixels=(double *)malloc(im.nx*im.ny*sizeof(double));
  im.pixelwrap=(double **)malloc(im.nx*sizeof(double *));
  for (i=0; i<im.nx; i++){
    im.pixelwrap[i]=&im.pixels[i*im.ny];}
  for (i=0; i<im.nx; i++){
    for (j=0; j<im.ny; j++) {
      im.pixelwrap[i][j]=(double) imasrow[j*im.nx + i];
    } // for j
  }//for i
  // close stuff
  fits_close_file(fptr, &status);
  free(imasrow);
  // can something really go wrong with closing?
  // and if it does do I care?

  return im;
}

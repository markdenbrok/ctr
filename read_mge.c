/************************************************************************
 * This code was written by Mark den Brok (Heidelberg 2011)             *
 * denbrok@astro.rug.nl                                                 *
 ************************************************************************/




#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"tr.h"
#include"readcol.h"


void read_mge(char *file, int length, struct multigaussexp *mge){
  int ntot;
  int *intarr;
  char *flags;
  char *fmt="%i %lf %lf %lf";
  flags=(char *)malloc((39+strlen(file))*sizeof(char));
  sprintf(flags,"quiet bufsize=1024 filelen=%i skipsym=# ",length);
  
  intarr=(int *)malloc(length*sizeof(int));
  mge->area=(double *)malloc(length*sizeof(double));
  mge->sigma=(double *)malloc(length*sizeof(double));
  mge->q=(double *)malloc(length*sizeof(double));
  ntot=readcol(file,flags,fmt,intarr,mge->area,mge->sigma,mge->q);
  mge->ntotal=ntot;
  free(flags);
  return;
}


void free_mge(struct multigaussexp *mge){
  free(mge->area);
  free(mge->sigma);
  free(mge->q);


}

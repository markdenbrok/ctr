/************************************************************************
 * This version of tge code was written by Mark den Brok (Utah 2014)    *
 * denbrok@physics.utah.edu                                             *
 ************************************************************************/



#include<stdio.h>
#include"tr.h"



void print_mge(struct multigaussexp *mge){
  int i;
  printf("\n");
  for (i=0; i<mge->ntotal; i++)
    printf("%d %lf %lf %lf \n",i,mge->area[i],mge->sigma[i], mge->q[i]);
  printf("\n");
}

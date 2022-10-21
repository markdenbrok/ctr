/************************************************************************
 * Original part of code was written by Mark den Brok (Heidelberg 2011) *
 * denbrok@astro.rug.nl                                                 *
 ************************************************************************/




#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>


#define elem_swap(a,b) { register double t=(a);(a)=(b);(b)=t;}
#define median(a,n) kth_smallest(a,n,(((n)&1)?((n)/2):(((n)/2)-1)))


double kth_smallest(double *a,int n, int k){
  int i,j,l,m;
  double x;
  l=0;
  m=n-1;
  while(l<m) {
    x=a[k];
    i=l;
    j=m;
    do {
      while (a[i]<x) i++;
      while (x<a[j]) j--;
      if (i<=j) {
	elem_swap(a[i],a[j]);
	i++;
	j--;
      }
    } while (i<=j);
    if (j<k) l=i;
    if (k<i) m=j;
  }
  return a[k];
}
    

double *where_array(double *input_array, int *goodvals,int totalgood){
  int i;
  double *retvals;
  retvals = (double *) malloc(totalgood * sizeof(double));
  if (retvals==NULL){
    printf("Error: failed in memory declaration\n");
    exit(1);
  }
  for (i=0; i<totalgood; i++){
    retvals[i]=input_array[goodvals[i]];
  }
  return retvals;
}


double max_in_array(double x[], int n){
  int i;
  double mx;
  mx=x[0];
  for (i=0; i<n; i++)
    {
    if (x[i]>mx)
      {
	mx=x[i];
      }
    }
  return mx;
}



int max_in_array_index(double x[], int n){
  int i, imax;
  double mx;
  mx=x[0];
  imax=0;
  for (i=0; i<n; i++)
    {
    if (x[i]>mx)
      {
	mx=x[i];
	imax=i;
      }
    }
  return imax;
}


double minimum(double *x, int n){
  int i;
  double mx=x[0];
  for (i=0; i<n; i++)
    {
    if (x[i]<mx)
      {
	mx=x[i];
      }
    }
  return mx;
}



double *range_arr(double x1, double x2, int n, int open)
{
  int i;
  double *returnval;
  returnval = (double *) malloc(n*sizeof(double));
  if (open) {
    for (i=0; i<n; i++){
      returnval[i]= x1 + (x2-x1)*(0.5+i)/n;
    }
  }
  else {
    if (n<0) {
      n=(int) (x2 - x1);
    }
    for (i=0; i<n; i++){
      returnval[i]= x1 + (x2-x1)*i/(n-1.0);
    }
  }
  return (double *) returnval;
}


double median2(double *x, int n) {
    double temp;
    int i, j;
    // the following two loops sort the array x in ascending order
    for(i=0; i<n-1; i++) {
        for(j=i+1; j<n; j++) {
            if(x[j] < x[i]) {
                // swap elements
                temp = x[i];
                x[i] = x[j];
                x[j] = temp;
            }
        }
    }
 
    if(n%2==0) {
        // if there is an even number of elements, return mean of the two elements in the middle
        return((x[n/2] + x[n/2 - 1]) / 2.0);
    } else {
        // else return the element in the middle
        return x[n/2];
    }
}



void printdarr(double *arr,int n,char *title){
  int i;
  printf("%s: ",title);
  for (i=0; i<n; i++){
    printf("%lf, ",arr[i]);
  }
  printf("\n");
  return;
}


char *stringconcat(char *a, char *b){
  int len1, len2;
  char *retstr;
  len1=strlen(a);
  len2=strlen(b);
  retstr=(char *)malloc((len1+len2+1)*sizeof(char));
  sprintf(retstr,"%s%s",a,b);
  return retstr;
}
  
float *d2farray(double *a,int n){
  float *res;
  int i;
  res=(float *)malloc(n*sizeof(float));
  for (i=0;i<n;i++){
    res[i]=(float)a[i];
  }
  return res;
}
  

//Returns a normally distributed deviate with zero mean and unit variance, using ran1(idum)
//as the source of uniform deviates.
// From NR
double gasdev()
{
  //double ran1(long *idum);
  static int iset=0;
  static double gset;
  double fac,rsq,v1,v2;
  if (iset == 0) {
    //We don’t have an extra deviate handy, so
    do {
      v1=2.0*drand48()-1.0;
      //v1=2.0*ran1(idum)-1.0;
    //pick two uniform numbers in the square ex-
      v2=2.0*drand48()-1.0;
      //v2=2.0*ran1(idum)-1.0;
      //tending from -1 to +1 in each direction,
      rsq=v1*v1+v2*v2;
    //see if they are in the unit circle,
  } while (rsq >= 1.0 || rsq == 0.0);

  fac=sqrt(-2.0*log(rsq)/rsq);
  gset=v1*fac;
  iset=1;
  return v2*fac;
 } 
 else {
  iset=0;
  return gset;
 }
}

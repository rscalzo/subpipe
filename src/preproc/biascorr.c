#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fitsio.h"

#include "biascorr.h"

float biascorr_side (fitsfile *infptr,fitsfile *outfptr, const int side,
		     const int *bias, const int *trim, const long *naxes, 
		     const long naxes1, const int fpixel0, const double *xx,
		     double *fitxx, double *fitarray, 
		     float *array, float *newarray, int *status) {
  long fpixel[2];
  double  nulval=0.;
  double devarray, medarray;
  int ii, jj, ind, ngood, nreject, ct, anynul;
  float meanbias;
  
  meanbias=0.;
  for (ii = 1; ii <= naxes[1]; ii++) {
    fpixel[0] = 1;
    fpixel[1] = ii;
    fits_read_pix(infptr, TFLOAT, fpixel, naxes[0], &nulval, array,
		  &anynul, status);
    ind=0;
    
    if (side==0) {
      for (jj = bias[0]; jj <=bias[1]; jj++) {
	fitarray[ind] = array[jj];
	fitxx[ind] = jj;
	ind++;
      }
    } else {
      for (jj = bias[2]; jj <= bias[3]; jj++) {
	fitarray[ind] = array[jj];
	fitxx[ind] = jj;
	ind++;
      }
    }
    /*fit */
    devarray=0;
    medarray=0;
    ngood=ind;
    nreject=ind;
    while (nreject >0) {
      medarray=kth_smallest_dbl(fitarray,ngood,ngood/2);
      for (jj = 0; jj < ngood ;jj++) {
	devarray +=(fitarray[jj]-medarray)*(fitarray[jj]-medarray);
      }
      devarray =sqrt(devarray/(ngood-1));
      ct=0;
      for (jj = 0; jj < ngood ;jj++) {
	if (abs(fitarray[jj]-medarray)/devarray <= 3.) {
	  fitxx[ct]=fitxx[jj];
	  fitarray[ct]=fitarray[jj];
	  ct++;
	}
      }
      nreject=ngood-ct;
      ngood =ct;
    }
    
    meanbias+=medarray;
    for (jj = 0; jj < naxes1; jj++) {
      newarray[jj] = array[jj+trim[0]]- medarray;
    }
    
    fpixel[0]=fpixel0;
    fits_write_pix(outfptr, TFLOAT, fpixel, naxes1, newarray,
		   status);
  }
  meanbias/=naxes[1];
  return meanbias;
}

float biascorr (fitsfile *infptr,fitsfile *outfptr,
	      const int *bias, const int *trim, const long *naxes, 
	      const long naxes1, const int fpixel0, const double *xx,
	      double *fitxx, double *fitarray, 
	      float *array, float *newarray, int *status) {
  long fpixel[2];
  double  nulval=0.;
  double c1;
  double devarray, medarray, midl, midr,midxl,midxr;
  int ii, jj, ind, ngood, nreject, ct, anynul;
  float meanbias;
  
  meanbias=0.;
  for (ii = 1; ii <= naxes[1]; ii++) {
    fpixel[0] = 1;
    fpixel[1] = ii;
    fits_read_pix(infptr, TFLOAT, fpixel, naxes[0], &nulval, array,
		  &anynul, status);
    ind=0;
    /*left median*/
    for (jj = bias[0]; jj <=bias[1]; jj++) {
      fitarray[ind] = array[jj];
      fitxx[ind] = jj;
      ind++;
    }
    midxl=(bias[0]+bias[1])/2;
    devarray=0;
    medarray=0;
    ngood=ind;
    nreject=ind;
    while (nreject >0) {
      medarray=kth_smallest_dbl(fitarray,ngood,ngood/2);
      for (jj = 0; jj < ngood ;jj++) {
	devarray +=(fitarray[jj]-medarray)*(fitarray[jj]-medarray);
      }
      devarray =sqrt(devarray/(ngood-1));
      ct=0;
      for (jj = 0; jj < ngood ;jj++) {
	if (abs(fitarray[jj]-medarray)/devarray <= 3.) {
	  fitxx[ct]=fitxx[jj];
	  fitarray[ct]=fitarray[jj];
	  ct++;
	}
      }
      nreject=ngood-ct;
      ngood =ct;
    }
    midl=medarray;
    
    ind=0;
    /*right*/
    for (jj = bias[2]; jj <= bias[3]; jj++) {
      fitarray[ind] = array[jj];
      fitxx[ind] = jj;
      ind++;
    }    
    midxr=(bias[2]+bias[3])/2;
    devarray=0;
    medarray=0;
    ngood=ind;
    nreject=ind;
    while (nreject >0) {
      medarray=kth_smallest_dbl(fitarray,ngood,ngood/2);
      for (jj = 0; jj < ngood ;jj++) {
	devarray +=(fitarray[jj]-medarray)*(fitarray[jj]-medarray);
      }
      devarray =sqrt(devarray/(ngood-1));
      ct=0;
      for (jj = 0; jj < ngood ;jj++) {
	if (abs(fitarray[jj]-medarray)/devarray <= 3.) {
	  fitxx[ct]=fitxx[jj];
	  fitarray[ct]=fitarray[jj];
	  ct++;
	}
      }
      nreject=ngood-ct;
      ngood =ct;
    }
    midr=medarray;
    
    c1=(midr-midl)/(midxr-midxl);
    
    meanbias+=midl +c1* (xx[naxes1/2]-midxl);
    for (jj = 0; jj < naxes1; jj++) {
      newarray[jj] = array[jj+trim[0]]- midl - c1* (xx[jj]-midxl);
    }
    
    
    fpixel[0]=fpixel0;
    fits_write_pix(outfptr, TFLOAT, fpixel, naxes1, newarray,
		   status);
  }
  meanbias/=naxes[1];
  return meanbias;
}


/* Use the STL algorithm if using C++ compiler  */
 
#if defined(_cplusplus) || defined(__cplusplus)

#include <algorithm>

double kth_smallest_dbl(double *a, int n, int k) 
{
  double *pk = a + k;
  std::nth_element(a, pk, (a  + n));
  return *pk;
}

float kth_smallest_flt(float *a, int n, int k) 
{
  float *pk = a + k;
  std::nth_element(a, pk, (a  + n));
  return *pk;
}

#else

double kth_smallest_dbl(double *a, int n, int k) 
{
  int i,j,l,m ;
  double x,temp;

  l=0 ; m=n-1 ; 
  while (l<m) { 
    x=a[k] ; 
    i=l ; 
    j=m ; 
    do { 
      while (a[i]<x) i++ ; 
      while (x<a[j]) j-- ; 
      if (i<=j) {
	temp=a[i];
	a[i]=a[j];
	a[j]=temp;
	i++ ; j-- ; 
      } 
    } while (i<=j) ; 
    if (j<k) l=i ; 
    if (k<i) m=j ; 
  } 
  return a[k] ; 
} 

float kth_smallest_flt(float *a, int n, int k) 
{ 
  /* Code commented by RS 2012/04/18.
   * Start off with array a of length n and an integer k < n.
   */
  int i,j,l,m ;
  float x,temp;

  l=0 ; m=n-1 ; 
  while (l<m) { 
    /* start i and j at the ends of the interval pick x as a "pivot" element */
    x=a[k] ; 
    i=l ; 
    j=m ; 
    do { 
      /* look for two elements which belong on opposite sides of x */
      while (a[i]<x) i++ ; 
      while (x<a[j]) j-- ; 
      if (i<=j) {
        /* swap i and j */
        temp=a[i];
        a[i]=a[j];
        a[j]=temp;
        i++ ; j-- ; 
      } 
    } while (i<=j) ;
    /* k now separates the array into elements greater than a[k], and elements
     * smaller than a[k].  If all the elements live on one side of k, restrict
     * all future iterations to this smaller space.
     * This way of doing it may be more efficient than sorting the whole array
     * because you only care where the nth largest element is.  But I bet that
     * it isn't that much more efficient than an N log N sort.
     */
    if (j<k) l=i ; 
    if (k<i) m=j ; 
  } 
  return a[k] ; 
}

#endif

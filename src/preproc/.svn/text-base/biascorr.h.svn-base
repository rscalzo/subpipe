#ifndef _BIASCORR_H_GUARD
#define _BIASCORR_H_GUARD

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fitsio.h"

double kth_smallest_dbl(double *a, int n, int k);
float kth_smallest_flt(float *a, int n, int k);

float biascorr_side (fitsfile *infptr,fitsfile *outfptr, const int side,
		     const int *bias, const int *trim, const long *naxes, 
		     const long naxes1, const int fpixel0, const double *xx,
		     double *fitxx, double *fitarray, 
		     float *array, float *newarray, int *status);

float biascorr (fitsfile *infptr,fitsfile *outfptr,
		const int *bias, const int *trim, const long *naxes, 
		const long naxes1, const int fpixel0, const double *xx,
		double *fitxx, double *fitarray, 
		float *array, float *newarray, int *status);

#endif

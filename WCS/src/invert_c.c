#include <stdio.h>
#include <math.h>

#define MAXSIZE 1000

#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define ABS(a) (((a) > 0) ? (a) : -(a))

invert(int n, double *a, double *det)
{

/*
      SUBROUTINE INVERT(N,A,RST,DET)

     Subroutine to invert a matrix, and compute the determinant.
     John Tonry, 9/2/80; transcribed to C 020904.

     INVERT's arguments:

     N - The dimension of the matrix
     A - The NxN matrix to be inverted. Upon successful inversion, A
         contains the inverse. A must be real*8
     DET - The determinant of the matrix. DET is set to 0 for a
         singular matrix, and in that case, A contains garbage.

     return values are 0 -- normal return
                      -2 -- matrix too large
                      -1 -- singular matrix
*/
   double save, pivot, onrow, cprev, cnow, decr;
   short int rst[MAXSIZE][2];
   int i, j, k, l, mrank, isign, nrow, ncol;

   if(n > MAXSIZE) {
      fprintf(stderr, "Matrix size too big: %d > %d\n", n, MAXSIZE);
      return(-2);
   }
   mrank = 0;
   isign = 1;
   *det = 0.0;
   for(j=0; j<n; j++) rst[j][0] = rst[j][1] = -1;
/* Loop over columns, reducing each */
   for(i=0; i<n; i++) {

/* Find the pivot element */
      pivot = -1.0;
      for(j=0; j<n; j++) {
	 if(rst[j][0] != -1) continue;
	 for(k=0; k<n; k++) {
	    if(rst[k][0] != -1) continue;
	    if(pivot >= ABS(a[j+k*n])) continue;
	    pivot = ABS(a[j+k*n]);
	    nrow = j;
	    ncol = k;
	 }
      }
      pivot = a[nrow+ncol*n];
      if(pivot == 0.0) {
	 *det = 0;
	 return(-1);
      }
      rst[ncol][0] = nrow;
      rst[ncol][1] = i;
/* Swap pivot element onto the diagonal */
      for(k=0; k<n; k++) {
	 save = a[nrow+k*n];
	 a[nrow+k*n] = a[ncol+k*n];
	 a[ncol+k*n] = save;
      }
/*   Reduce pivot column */
      for(j=0; j<n; j++) a[j+ncol*n] = -a[j+ncol*n]/pivot;
      a[ncol+ncol*n] = 1/pivot;

/*   Reduce other columns */
      for(k=0; k<n; k++) {
	 if(k == ncol) continue;
/*     Find maximum of column to check for singularity */
	 cprev = 0;
	 for(j=0; j<n; j++) cprev = MAX(cprev,ABS(a[j+k*n]));
/*     Reduce the column */
	 onrow = a[ncol+k*n];
	 a[ncol+k*n] = 0;
	 for(j=0; j<n; j++) a[j+k*n] = a[j+k*n] + onrow*a[j+ncol*n];
/*     Find the new maximum of the column */
	 cnow = 0;
	 for(j=0; j<n; j++) cnow = MAX(cnow,ABS(a[j+k*n]));

/*     Quit if too many figures accuracy were lost (singular) */
	 decr = cprev / cnow;
	 if(cnow == 0.0 || decr > 1e8) {
	    *det = 0;
	    return(-1);
	 }
      }
      *det = *det + log(ABS(pivot));
      if(pivot < 0) isign *= -1;
      mrank++;
   }

/*     Now untangle the mess */
   for(j=0; j<n; j++) {
      for(k=0; k<n; k++) {
	 if(rst[k][1] != (n-1-j)) continue;
	 ncol = rst[k][0];
	 if(ncol == k) break;
	 for(l=0; l<n; l++) {
	    save = a[l+ncol*n];
	    a[l+ncol*n] = a[l+k*n];
	    a[l+k*n] = save;
	 }
	 break;
      }
   }
	 
   if(ABS(*det) < 88) *det = isign * exp(*det);
   return(0);
}

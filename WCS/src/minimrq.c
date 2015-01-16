#include <stdio.h>
#include <math.h>

#define MAXPAR  16		/* Max number of parameters */
#define MAXITER 20		/* Max iterations */
#define QFRAC   1e-8		/* Fractional change for quitting */

#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define ABS(a) (((a) > 0) ? (a) : -(a))


mini(int npar, double *par, int *usepar,
     int func(int npar, double *par, int *usepar, double *value,
	      int doderiv, double *deriv, double *curve),
     double *cov, int nderiv, int *lambda, int verbose)
{
/*
*     This is John Tonry's black box minimization (and fitting) program,
*     implementing a Marquard algorithm.
*     Transcribed to C, 020904.
*     Revision 2.0, 11/17/82.
*
*     MINI's arguments are as follows:
*     NPAR - The number of parameters to be varied in searching for a
*        minimum.
*     PAR - A double vector, it serves three purposes:
*        (1) It passes MINI the initial estimate for the parameters
*        (2) It passes arguments to the function to be minimized
*        (3) It returns the minimum parameter values
*     USEPAR - set to 0/1 to turn on that parameter for variation
*     FUNC - The function to be minimized. FUNC takes NPAR, PAR and 
*        USEPAR as arguments and return the value of the function as
*        VALUE.  If DODERIV = 0, this is all that is wanted.  If
*        DODERIV = 1, the first partial derivatives wrt the active
*        parameters should be returned in DERIV (close packed).  If
*        DODERIV = 2, the second partial derivatives wrt the active
*        parameters should be returned in CURVE (close packed).
*     COV - A NxN real*8 matrix in which the covariance matrix of the
*        fit is returned.
*     NDERIV - Variable to govern how MINI gets its derivatives:
*        NDERIV = 0 for a function of arbitrary A and numerical derivatives.
*        NDERIV = 1 for a function of arbitrary A, analytic first derivatives
*                  provided by DFUNK, and numerical second derivatives.
*        NDERIV = 2 for a function of arbitrary A, analytic first and second
*                  derivatives provided by DFUNK and by D2FUNK.
*        NDERIV = 3 for a function that is quadratic in A, and MINI
*                  will iterate once, computing the minimum exactly.
*     VERBOSE - governs whether MINI will print each iteration
*
*     LAMBDA - starting value for lambda (suggest -2), returns final value
*
*     Return value - 0 successful completion
*                   -1 error from invert: singular matrix
*                   -2 error from invert: matrix too large
*                   -3 too many parameters
*                    N maximum iteration reached (N)
*
*     Descriptions of some of the variables:
*     PAR - Argument for the function
*     A0 - Current guess for the minimum
*     AI - increments for a0 in computing derivatives
*     DA - Vector from A0 to new guess for the minimum
*     DF - First derivatives of the function
*     D2F - Second derivatives of the function
*     LAMBDA - Governs mix of gradient and analytic searches
*     ITER - Maximum number of iterations
*     QFRAC - Maximum fractional change for successful exit
*
*     The calling program should be as follows (eg):
*********************************************************
*     REAL*8 A(4)
*     REAL*8 COV(4,4)
*     EXTERNAL CHI, DCHI, D2CHI
*     ... (Initialize A to the guess for the minimum)
*     CALL MINI(4,A,CHI,DCHI,D2CHI,COV,NDERIV)
*     ...
*     FUNCTION CHI(A)
*     REAL*8 A(4), CHI
*     ... (define the function)
*     SUBROUTINE DCHI(A,DF)
*     REAL*8 A(4), DF(4)
*     ... (dummy or compute the derivatives)
*     SUBROUTINE D2CHI(A,COV)
*     REAL*8 A(4), COV(4,4)
*     ... (dummy or compute the derivatives)
*************************************************************
*/
   double df[MAXPAR], a0[MAXPAR], da[MAXPAR], ai[MAXPAR];
   double d2f[MAXPAR*MAXPAR], fnow, fthen, fminus, curve, det, err;
   double dfrac=0.02, places=1e-7, vary=0.1, base=10.0;
   int on[MAXPAR];
   int i, j, k, lam, nuse, error, iter, lamit;
      
   if(npar > MAXPAR) {
      fprintf(stderr, "Too many parameters, max = %d\n", MAXPAR);
      return(-3);
   }

/*     Define a few parameters */
   lam = *lambda;
   for(i=0, nuse=0; i<npar; i++) {
      if(usepar[i]) {
	 on[nuse] = i;
	 nuse++;
      }
   }

/*     If NDERIV = 3, compute the minimum directly and exit. */
   if(nderiv == 3) {
      for(i=0; i<nuse; i++) par[on[i]] = 0.0;
      error = func(npar, par, usepar, &fnow, 0, NULL, NULL);
      for(i=0; i<nuse; i++) {
	 par[on[i]] = 1.0;
	 error = func(npar, par, usepar, &df[i], 0, NULL, NULL);
	 for(j=0; j<=i; j++) {
	    par[on[j]] = par[on[j]] + 1;
	    error = func(npar, par, usepar, &d2f[i+nuse*j], 0, NULL, NULL);
	    d2f[i+nuse*j] += fnow - df[i] - df[j];
	    cov[i+nuse*j] = cov[j+nuse*i] = d2f[j+nuse*i] = d2f[i+nuse*j];
	    par[on[j]] = 0.0;
	 }
	 par[on[i]] = 0.0;
      }
      for(i=0; i<nuse; i++) df[i] -= 0.5*d2f[i+i*nuse] + fnow;

      if( (error=invert(nuse, cov, &det)) != 0) return(error);

      for(i=0; i<nuse; i++) {
	 for(j=0, par[on[i]]=0.0; j<nuse; j++) 
	    par[on[i]] -= cov[i+nuse*j] * df[j];
      }
      error = func(npar, par, usepar, &fnow, 0, NULL, NULL);

      goto DOCOVAR;
      fprintf(stderr, "minimrq: Should never see this...\n");
   }

/*     Initialize A0 */
   for(i=0; i<npar; i++) a0[i] = par[i];
   error = func(npar, par, usepar, &fnow, 0, NULL, NULL);

   if(verbose) vprint(npar, a0, fnow, 0, lam);

/*     Initialize AI */
   for(i=0; i<nuse; i++) {
      ai[on[i]] = ABS(vary*a0[on[i]]);
      if(ai[on[i]] == 0.0) ai[on[i]] = vary;
   }

/*     Begin iteration to find minimum */
   for(iter=0; iter<MAXITER; iter++) {

      fthen = fnow;

/* Compute the function derivatives. */
      for(j=0; j<nuse; j++) par[on[j]] = a0[on[j]];

      switch(nderiv) {

/* First case: NDERIV = 0 so entirely numerical derivatives are required */
	 case 0:
/* First the 1st derivatives */
	    for(j=0; j<nuse; j++) {
	       par[on[j]] = a0[on[j]] + ai[on[j]];
	       error = func(npar, par, usepar, &df[j], 0, NULL, NULL);
	       par[on[j]] = a0[on[j]];
	    }
/* Next, get the off diagonal 2nd derivatives. */
	    for(j=1; j<nuse; j++) {
	       for(k=0; k<j; k++) {
		  par[on[k]] = a0[on[k]] + ai[on[k]];
		  par[on[j]] = a0[on[j]] + ai[on[j]];
		  error = func(npar, par, usepar, &d2f[j+k*nuse], 0, NULL, NULL);
		  d2f[j+k*nuse] = d2f[k+j*nuse] = 
		     (d2f[j+k*nuse]-df[k]-df[j]+fnow) / (ai[on[j]]*ai[on[k]]);
		  
		  par[on[k]] = a0[on[k]];
		  par[on[j]] = a0[on[j]];
	       }
	    }
/* Finally do the on diagonal 2nd derivatives, and fix the 1st ones. */
	    for(j=0; j<nuse; j++) {
	       par[on[j]] = a0[on[j]] - ai[on[j]];
	       error = func(npar, par, usepar, &fminus, 0, NULL, NULL);
	       d2f[j+j*nuse] = (fminus+df[j]-2*fnow) / (ai[on[j]]*ai[on[j]]);
	       df[j] = (df[j] - fminus) / (2*ai[on[j]]);
	       par[on[j]] = a0[on[j]];
	    }
	    break;
	 
/* Second case: NDERIV = 1 so analytic first derivatives are available */
	 case 1:
	    error = func(npar, par, usepar, &fminus, 1, df, NULL);
	    for(j=0; j<nuse; j++) {
	       par[on[j]] = a0[on[j]] + ai[on[j]];
	       error = func(npar, par, usepar, &fminus, 1, da, NULL);
	       par[on[j]] = a0[on[j]];
	       for(i=0; i<=j; i++) {
		  d2f[i+j*nuse] = d2f[j+i*nuse] = (da[i]-df[i]) / ai[on[j]];
	       }
	    }
	    break;

/* Third case: NDERIV = 2 so analytic derivatives are available */
	 case 2:
	    error = func(npar, par, usepar, &fminus, 2, df, d2f);
	    break;

	 default:
	    fprintf(stderr, "unknown nderiv = %d\n", nderiv);
	    return(-3);
      }

/* Compute better estimates for the increments. */
      for(j=0; j<nuse; j++) {
	 curve = d2f[j+j*nuse];
	 if(curve == 0.0) curve = vary;
	 ai[on[j]] = sqrt(pow(df[j]*dfrac/curve,2.0)+ABS(fnow*places/curve));
      }
/*     Begin loop to find a direction along which function decreases */
      for(lamit=0; lamit<MAXITER; lamit++) {

/*     Get weight matrix */
	 for(j=0; j<nuse; j++) {
	    for(i=0; i<j; i++) cov[i+j*nuse] = cov[j+i*nuse] = d2f[i+j*nuse];
	    cov[j+j*nuse] = ABS(d2f[j+j*nuse]*(1+pow(base,(double)lam)));
	 }
	 if( (error=invert(nuse, cov, &det)) != 0) return(error);

/*     Multiply to get dA */
	 for(j=0; j<nuse; j++) {
	    for(i=0,da[j]=0.0; i<nuse; i++) da[j] -= cov[j+i*nuse]*df[i];
	 }
/*     Now get new function value */
	 for(j=0; j<nuse; j++) par[on[j]] = a0[on[j]] + da[j];
	 error = func(npar, par, usepar, &fnow, 0, NULL, NULL);

/*
 *     Test for whether the function has decreased
 *     If so, adopt the new point and decrement LAMBDA
 *     Else, increment LAMBDA, and get a new weight matrix
 */
	 if(fnow < fthen) break;
	 lam++;
      }
/*     Normal exit, the function at A0 + DA is less than at A0 */
      for(j=0; j<nuse; j++) a0[on[j]] = par[on[j]];
      lam--;
/*
 *     Print the current status and test to see whether the function
 *     has varied fractionally less than QFRAC.
 */
      if(verbose) vprint(npar, a0, fnow, iter, lam);
      if(ABS(fthen-fnow)/fnow < QFRAC) break;
   }
/*     This is the final computation of the covariance matrix */
/*     Quit if no minimum was found in the allowed number of iterations */

DOCOVAR:

   if(iter >= MAXITER) {
      if(verbose) printf("Maximum iteration exceeded\n");
      error = iter;
   } else {
      error = 0;
   }
/*     Finally, compute the covariance matrix */
   for(j=0; j<nuse; j++) {
      for(k=0; k<nuse; k++) cov[k+j*nuse] = d2f[k+j*nuse] / 2;
   }

   if( (error=invert(nuse, cov, &det)) != 0) return(error);
   for(j=0; j<nuse; j++) {
      err = sqrt(ABS(cov[j+j*nuse]));
      if(cov[j+j*nuse] < 0) err = -err;
      cov[j+j*nuse] = err;
   }

   for(j=1; j<nuse; j++) {
      for(k=0; k<j; k++) {
	 if(!usepar[k]) continue;
	 cov[k+j*nuse] /= cov[k+k*nuse]*cov[j+j*nuse];
	 cov[j+k*nuse] = cov[k+j*nuse];
      }
   }
   *lambda = lam;
   return(error);
}

vprint(int n, double *a, double f, int niter, int lambda)
{
   int i;
   for(i=0; i<n; i++) {
      if((i%7) == 0) printf(" A(I) =");
      printf("%11.4g", a[i]);
      if((i%7) == 6) printf("\n");
   }
   if((i%7) != 6) printf("\n");
   printf(" F = %14.7g     ITER = %3d    LAMBDA = %3d\n", f, niter, lambda);
   return(0);
}

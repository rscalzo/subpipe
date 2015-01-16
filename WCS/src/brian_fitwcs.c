#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <cirdr.h>
#include <cir_wcssubs.h>

#define MAXPAR 20
#define RAD2DEG 180.0/M_PI
#define TWOPI 2*M_PI
#define INITALLOC 64
#define MAXLIST 1000
#define MAXGRID 10000

static double bps_ra[10000],bps_dec[10000],bps_x[10000],bps_y[10000];
static double bps_eta[10000],bps_xi[10000];
static double residra[10000],residdec[10000];
static int bps_usept[10000],bps_npts,calcResXflag=1;
struct wcsprm *wcs = NULL;
static int verbosity;
static char cvsid[] = "";
static char errmsg[120];

/* RS 2011/11/08:  Modified slightly to ensure that we return a zero value
   upon successful termination (this is important for pipeline operations). */

main (argc,argv)
     int argc;
     char *argv[];
{
  double par[MAXPAR], cov[MAXPAR*MAXPAR];
  double parx[MAXPAR],pary[MAXPAR];
  double rpar[MAXPAR], rcov[MAXPAR*MAXPAR];
  double xgridcoord[MAXGRID],ygridcoord[MAXGRID],ragridcoord[MAXGRID],decgridcoord[MAXGRID];
  double xierr,etaerr;
  double dummyderiv,dummycurve;
  int npar=MAXPAR, usepar[MAXPAR], nderiv=0, lambda=-2, verbose=0;
  int nrpar=MAXPAR, rusepar[MAXPAR];
  int calcerr(int npar, double *par, int *usepar, double *value,
	   int doderiv, double *deriv, double *curve);
  int calcresid(int npar, double *par, int *usepar, double *value,
	   int doderiv, double *deriv, double *curve);
  int calcerrlin(int npar, double *par, int *usepar, double *value,
	   int doderiv, double *deriv, double *curve);
  int applyresid(int npar, double *par);
  int xy2wcs(int npar, double *parx, double *pary, int nlist, double *xlist, double *ylist, double *ra, double *dec);
  double fnow,v1,v2,v3,v4,v5,v6;
  double a,b,c,d,e,f,xi,eta;
  double x,y,xgrid,ygrid,dxgrid,dygrid;
  double ra[MAXLIST], dec[MAXLIST],xlist[MAXLIST],ylist[MAXLIST];
  double crpix1,crpix2,cd1_1,cd1_2,cd2_1,cd2_2;
  double pc1_1,pc1_2,pc2_1,pc2_2,cdelt1,cdelt2;
  int i,retval,npts,err,fitpv,fitpv5,itersig,fitcrval,fitresid,nlist,outputgrid;
  FILE *fd;
  char msg[BUFSIZ],image[100],posfile[100],line[150];
  char outfile[100],outresidfile[100],outgridfile[100];
  float sigmatch,floormatch;
  int nexclude,iter,maxit,usemedsig;
  void syntax(),residcalc(),medianerr(),printobj(),transform(),gausserr();
  int excludeobjs();
 
  
  if (argc < 3) {
    syntax(argv[0]);
    exit(-1);
  }
  
  sprintf(image,argv[1]);
  sprintf(posfile,argv[2]);

  /* default values for everything*/

  strcpy(outfile,posfile);
  strcat(outfile,".fit");

  strcpy(outresidfile,posfile);
  strcat(outresidfile,".resid");

  strcpy(outgridfile,posfile);
  strcat(outgridfile,".grid");

  sigmatch=3.;
  floormatch=0.06; 
  maxit = 3;
  verbosity=0;
  fitpv=0;
  fitpv5=0;
  itersig=0;
  fitcrval=0;
  usemedsig=1;
  fitresid=0;
  outputgrid=0;

  /* parse optional command line arguments */
  if (argc > 3) { 
    for (i=3;i<argc;i++) {


      /* Sigclip value */

      if ( strcasecmp(argv[i],"-sigclip") == 0 ) {
	if ( ++i >= argc ) {
	  fprintf(stderr,"-sigclip option requires argument: Sigma value for iteratively removing stars from fit\n");
	  syntax(argv[0]);
	  exit(-1);
	}
	sigmatch = atof(argv[i]);
	if (sigmatch==1) {
	  itersig=1;
	}
	if (sigmatch < 1.0000000) {
	  fprintf(stderr,"-sigmatch  must be more than 1 or =1 (where it iteratively raises the stddev [%f]\n",sigmatch);
	  syntax(argv[0]);
	  exit(-1);
	}
	continue;
      }
      
      /* maxnumber of iterations value */
      
      if ( strcasecmp(argv[i],"-maxit") == 0 ) {
	if ( ++i >= argc ) {
	  fprintf(stderr,"-maxit option requires argument: max number of iterations for iteratively removing stars from fit\n");
	  syntax(argv[0]);
	  exit(-1);
	}
	maxit = atoi(argv[i]);
	if (maxit < 1 || maxit > 100) {
	  fprintf(stderr,"maxit  must be between 1 an 100 [%d]\n",maxit);
	  syntax(argv[0]);
	  exit(-1);
	}
	continue;
      }

      /* verbosity */
      
      if ( strcasecmp(argv[i],"-verbose") == 0 ) {
	if ( ++i >= argc ) {
	  fprintf(stderr,"-verbose option requires argument: level of verbosity 0,1,2\n");
	  syntax(argv[0]);
	  exit(-1);
	}
	verbosity = atoi(argv[i]);
	if (verbosity>0) {verbose=1;}
	if (verbosity < 0 && verbosity > 2) {
	  fprintf(stderr,"verbosity  must be between 0 and 2 [%d]\n",verbosity);
	  syntax(argv[0]);
	  exit(-1);
	}
	continue;
      }

      if ( strcasecmp(argv[i],"-usegausssig") == 0 ) {
	usemedsig=0;
 	continue;
      }

      /*do a residual calc at the end*/
      if ( strcasecmp(argv[i],"-fitresid") == 0 ) {
	fitresid=1;
 	continue;
      }

      /* fit non-linear terms? */
      
      if ( strcasecmp(argv[i],"-fitpv") == 0 ) {
	fitpv=1;
 	continue;
      }

      /* fit non-linear terms 5th order? */
      
      if ( strcasecmp(argv[i],"-fitpv5") == 0 ) {
	fitpv5=1;
 	continue;
      }

      if ( strcasecmp(argv[i],"-fitcrval") == 0 ) {
	fitcrval=1;
 	continue;
      }
     

      /* floormatch value */
      
      if ( strcasecmp(argv[i],"-floormatch") == 0 ) {
	if ( ++i >= argc ) {
	  fprintf(stderr,"-floormatch option requires argument: Minimum value for removing objects\n");
	  syntax(argv[0]);
	  exit(-1);
	}
	floormatch = atof(argv[i]);
	if (floormatch < 0) {
	  fprintf(stderr,"-floormatch  must be more than 0 [%f]\n",floormatch);
	  syntax(argv[0]);
	  exit(-1);
	}
	continue;
      }
      


      /* outputgrid */
      
      if ( strcasecmp(argv[i],"-outputgrid") == 0 ) {
	outputgrid=1;
	if ( ++i >= argc ) {
	  fprintf(stderr,"-outputgrid  option requires 4 arguments: endpixvalx endpixvaly dx dy\n");
	  syntax(argv[0]);
	  exit(-1);
	}
	xgrid = atof(argv[i++]);
	ygrid = atof(argv[i++]);
	dxgrid = atof(argv[i++]);
	dygrid = atof(argv[i++]);

	if (dxgrid*dygrid > MAXGRID) {
	  fprintf(stderr,"-outputgrid has too many points [%f] MAX allowed [%d]\n",dxgrid*dygrid,MAXGRID);
	  syntax(argv[0]);
	  exit(-1);
	}
	continue;
      }
      
      
      /*Outputfile */
      if ( strcasecmp(argv[i],"-outfile") == 0 ) {
	if ( ++i >= argc ) {
	  fprintf(stderr,"-outfile option requires argument: filename to output matched list\n");
	  syntax(argv[0]);
	  exit(-1);
	}
	strcpy(outfile,argv[i]);
	continue;
      }

      /* unrecognized command line option */
      fprintf(stderr,"unrecognized command line option: %s\n", argv[i]);
      syntax(argv[0]);
      exit(-1);
    }     

  }



  /* Open the fits image and read the current WCS info into WCS structure */
  
  retval = cir_wcsopen(image,&wcs,msg);
  if (retval != CIR_OK) {
    printf("PLATESOL: Couldn't read file input WCS %s -- %s\n", image,msg);
    exit(-1);
  }

  /* read in file that contains x,y on image and then RA,DEC in degrees*/

  if ((fd = fopen(posfile,"r")) == NULL) {
    printf("Brian_fitwcs: No input position file: ");
    exit(-1);
  }

  i=npts = 0;
  while(feof(fd)==0) {
    fgets(line,150, fd);
    if(sscanf(line,"%lf %lf %lf %lf",&v1,&v2,&v3,&v4)==4) {
      bps_x[i] = v1;
      bps_y[i] = v2;
      bps_ra[i]= v3;
      bps_dec[i]=v4;
      bps_usept[i]=1;
      cir_radectoxieta(wcs,v3,v4,&xi,&eta);
      bps_xi[i]=xi*180/3.1415926535897*3600;
      bps_eta[i]=eta*180/3.1415926535897*3600;
      i++;
    }
    else {
      printf("skipped: %70s\n",line);
    }
  }
  bps_npts=i-1;
  close(fd);

 
  /*first lets find a linear solution*/
  transform(bps_npts,bps_xi,bps_eta,bps_x,bps_y,&c,&a,&b,&f,&d,&e);

  cd1_1 = a/3600;
  cd1_2 = b/3600;
  cd2_1 = d/3600;
  cd2_2 = e/3600;
  crpix1 = (e*c - b*f)/(d*b - e*a);
  crpix2 = (a*f - d*c)/(d*b - e*a);

  
  if (verbosity > 0) {
    fprintf (stderr,"LINEAR FIT: %9.5e %9.5e %9.5e %9.5e %9.5e %9.5e\n",cd1_1,cd1_2,cd2_1,cd2_2,crpix1,crpix2);
  }
  
  /*put this fit into the wcs structure*/

  wcs->cd[0]=cd1_1  ; /*cd1_1*/
  wcs->cd[1]=cd1_2  ; /*cd1_2*/
  wcs->cd[2]=cd2_1  ; /*cd2_1*/
  wcs->cd[3]=cd2_2  ; /*cd2_2*/
  wcs->crpix[0]=crpix1  ; /*crpix1*/
  wcs->crpix[1]=crpix2  ;/*crpix2*/
  wcsset(wcs);
  
  /* now do the minimisation if it is a non-linear fit*/
 
  /*we can find up to 7 parameter fits to crpix, and the cd matrix, and the pv3 term*/
  

  npar=10;
  for(i=0; i<npar; i++) {
    usepar[i] = 1;
  }
  nderiv=0.;
    
  /*use linear fit as first guess plus whatever is in the header*/
    

  par[0] = wcs->cd[0];/*cd_1_1*/
  par[1] = wcs->cd[1];/*cd1_2*/
  par[2] = wcs->cd[2];/*cd2_1*/
  par[3] = wcs->cd[3];/*cd2_2*/
  par[4] = wcs->crpix[0]; /*crpix1*/
  par[5] = wcs->crpix[1]; /*crpix2*/
  if (strcmp(wcs->ctype[0],"RA---ZPN")==1) {
    par[6]= wcs->pv[3].value; /*pv_3*/   
    par[7]= wcs->pv[5].value; /*pv_5*/   
  }
  else {
    par[6]=par[7]=0.;
  }
  par[8]= wcs->crval[0]; /**/   
  par[9]= wcs->crval[1]; /**/   


  if (par[6] < 0.01) { par[6]=0.0;}
  if (par[7] < 0.01) { par[7]=0.0;}



  if (verbosity>1) {
    fprintf (stderr,"WCS READ IN: %9.5e %9.5e %9.5e %9.5e %9.5e %9.5e %9.5e %9.5e\n",par[0],par[1],par[2],par[3],par[4],par[5],par[6],par[7]);
  }

  npar=11;
  for(i=0; i<npar; i++) {
    usepar[i] = 0;
  }
  nderiv=0.;
  
  if (fitpv==1) {
    usepar[6]=1;
    usepar[7]=0;  
  }
    
  if (fitpv5==1) {
    fitpv=1;
    usepar[7]=1;
    usepar[6]=1;
  }
  

  if (fitcrval==1) {
    if (verbosity > 1) {
      fprintf(stderr,"Doing fit for crval\n");
    }
    usepar[8]=1;
    usepar[9]=1;
  }
  

  /* iterate until we converge or too many times*/
  nexclude=1;iter=0;
  while((nexclude>0 && iter < maxit) || (itersig==1) && sigmatch<3.5) {
    iter++;
  
    if (fitpv==1 || fitcrval==1) { /*do non-linear fit if we are fitting pv terms*/
      if (verbosity > 1) {
	fprintf(stderr,"Doing Non-linear fit\n");
      }
      /*      err = mini(npar, par, usepar, calcerr, cov, nderiv, &lambda, verbose);*/
      lambda=-1;
      err = mini(npar, par, usepar, calcerrlin, cov, nderiv, &lambda, verbose);
      lambda=-1;
    }
    else { /*do linear fit*/
        transform(bps_npts,bps_xi,bps_eta,bps_x,bps_y,&c,&a,&b,&f,&d,&e);
	
	cd1_1 = a/3600;
	cd1_2 = b/3600;
	cd2_1 = d/3600;
	cd2_2 = e/3600;
	crpix1 = (e*c - b*f)/(d*b - e*a);
	crpix2 = (a*f - d*c)/(d*b - e*a);
	

	cdelt1=sqrt(cd1_1*cd1_1+cd2_1*cd2_1);
	cdelt2=sqrt(cd1_2*cd1_2+cd2_2*cd2_2);
	pc1_1=cd1_1/cdelt1;
	pc1_2=cd1_2/cdelt1;
	pc2_1=cd2_1/cdelt2;
	pc2_2=cd2_2/cdelt2;
	
	
	if (verbosity > 0) {
	  fprintf (stderr,"\nCDMATRIX: %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f\n",cd1_1,cd1_2,cd2_1,cd2_2,crpix1,crpix2);
	  fprintf (stderr,"\nPCMATRIX: %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f\n",pc1_1,pc1_2,pc2_1,pc2_2,cdelt1,cdelt2,crpix1,crpix2);
	  fprintf (stderr,"LINEAR FIT: %9.5e %9.5e %9.5e %9.5e %9.5e %9.5e\n",cd1_1,cd1_2,cd2_1,cd2_2,crpix1,crpix2);
	}
    
	/*put this fit into the wcs structure*/
	
	wcs->cd[0]=cd1_1  ; /*cd1_1*/
	wcs->cd[1]=cd1_2  ; /*cd1_2*/
	wcs->cd[2]=cd2_1  ; /*cd2_1*/
	wcs->cd[3]=cd2_2  ; /*cd2_2*/
	wcs->crpix[0]=crpix1  ; /*crpix1*/
	wcs->crpix[1]=crpix2  ;/*crpix2*/

	wcsset(wcs);

    }

    /*find offset of stars from fit*/
    residcalc(residra,residdec);

    /*find the 1-sig error of the fit*/
    /*using a robust method or not*/
    if (usemedsig==1) {
      medianerr(residra,residdec,&xierr,&etaerr);
    }
    else {
      gausserr(residra,residdec,&xierr,&etaerr);
    }
    
    if (verbosity > 0) {
      fprintf(stderr,"ERROR in XI: %7.3f ETA: %7.3f\n",xierr,etaerr);
    }
    /*exclude those objects outside of acceptable range*/
    /*return how many objects are ecluded in n*/
    (nexclude)=excludeobjs(residra,residdec,xierr,etaerr,sigmatch,floormatch);
    if (verbosity > 0) {
      fprintf(stderr,"Excluded %d objects\n",nexclude);
    }
    
    if (verbosity>1) {
      printobj(residra,residdec,"stderr");
    }
    if (itersig==1) {
      sigmatch+=0.5;
    }
    
  }
  
  printobj(residra,residdec,outfile);
  
  /*print out the parameters*/
  for(i=0; i<npar; i++) {
    printf ("%12.9f ",par[i]);
  }
  printf(" RMS: (xi)%7.3f (eta)%7.3f arcsec %d stars\n",xierr,etaerr,bps_npts-nexclude);


  /*do we want to calculate a residual fit*/
    nrpar=6;
    for(i=0; i<nrpar; i++) {
      rusepar[i] = 1;
      parx[i]=0.0;
      pary[i]=0.0;
    }

  if (fitresid==1) {
    
    /*find residuals from last fit*/
    residcalc(residra,residdec);
    
    /*run linear fitter*/
    lambda=-2;
    /*do it first in xi space*/

    /*set up the parameters*/
   

    calcResXflag=1;
    /*gives us the function value*/
    err = mini(nrpar, parx, rusepar, calcresid, rcov, 3, &lambda, verbose);

    c=parx[0];
    a=parx[1];
    b=parx[2];


    printf ("\nXI RESID-CORR: ");
    for(i=0; i<nrpar; i++) {
      printf ("%12.9f ",parx[i]);
    }
    calcresid(nrpar,parx,rusepar,&xierr,0,&dummyderiv,&dummycurve);
    xierr=sqrt(xierr);

    applyresid(nrpar,parx);

    /* next do it in eta space*/

    calcResXflag=0;

    err = mini(nrpar, pary, rusepar, calcresid, rcov, 3, &lambda, verbose);
    f=pary[0];
    d=pary[1];
    e=pary[2];
  
    cd1_1=a/3600.;
    cd1_2=b/3600.;
    cd2_1=d/3600.;
    cd2_2=e/3600.;
    crpix1=(e*c - b*f)/(d*b - e*a);
    crpix2=(a*f - d*c)/(d*b - e*a);
  
    cdelt1=sqrt(cd1_1*cd1_1+cd2_1*cd2_1);
    cdelt2=sqrt(cd2_1*cd2_1+cd1_2*cd1_2);
    pc1_1=cd1_1/cdelt1;
    pc1_2=cd1_2/cdelt1;
    pc2_1=cd2_1/cdelt2;
    pc2_2=cd2_2/cdelt2;
  
    printf ("\nCDMATRIX: %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f\n",cd1_1,cd1_2,cd2_1,cd2_2,crpix1,crpix2);
    printf ("\nPCMATRIX: %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f\n",pc1_1,pc1_2,pc2_1,pc2_2,cdelt1,cdelt2,crpix1,crpix2);

    printf ("\nETA RESID-CORR: ");
    for(i=0; i<nrpar; i++) {
      printf ("%12.9f ",pary[i]);
    }
    calcresid(nrpar,pary,rusepar,&etaerr,0,&dummyderiv,&dummycurve);

    etaerr=sqrt(etaerr);

    applyresid(nrpar,pary);

    printobj(residra,residdec,outresidfile);
    printf("\nAFTER RESID-CORR RMS: (xi)%7.3f (eta)%7.3f arcsec\n",xierr,etaerr);
  }
  if (outputgrid==1) {
    i=0;
    for (x=1;x<=xgrid;x+=(xgrid-1)/dxgrid) {
      for (y=1;y<=ygrid;y+=(ygrid-1)/dygrid) {
	xgridcoord[i]=x;
	ygridcoord[i]=y;
	i++;
      }
    }
    
    nlist=i;
    xy2wcs(nrpar, parx, pary, nlist, xgridcoord, ygridcoord, ragridcoord, decgridcoord);
    if ((fd = fopen(outgridfile,"w")) == NULL) {
      fprintf(stderr,"OUTGRIDFILE: Cannot open %s for writing\n",outgridfile);
      exit(-1);
    }
    for (i=0;i<nlist;i++) {
      fprintf(fd,"%9.3f %9.3f %14.8f %14.8f\n",xgridcoord[i], ygridcoord[i], ragridcoord[i], decgridcoord[i]);
    }
    close (fd);
    
  }
  /* RS:  we got here, so presumably we haven't crashed; return zero */
  exit(0);
 }

/*calculates the average residual of points to a WCS fit*/
/*here the whole shooting match is thrown at the minimiser */
/*used by minimrq program*/
/*doesn't work well so is not used anymore*/
int calcerr(int npar, double *par, int *usepar, double *value,
	int doderiv, double *deriv, double *curve)
{
  int i, k;
  double used=1,c=0,ra,dec,ra1,dec1;

  /*update the wcs structure with values*/
  wcs->cd[0]=par[0]  ; /*cd1_1*/
  wcs->cd[1]=par[1]  ; /*cd1_2*/
  wcs->cd[2]=par[2]  ; /*cd2_1*/
  wcs->cd[3]=par[3]  ; /*cd2_2*/
  wcs->crpix[0]=par[4]  ; /*crpix1*/
  wcs->crpix[1]=par[5]  ;/*crpix2*/
  wcs->pv[3].value=par[6]  ;/*pv_3*/   
  wcs->pv[5].value=par[7]  ;/*pv_5*/   
  wcs->crval[0]=par[8];
  wcs->crval[1]=par[9];


  wcsset(wcs);
  
  /*find residuals of all stars*/
  for(i=0;i<bps_npts;i++) {
    if (bps_usept[i]) {
      cir_xytoradec(wcs,bps_x[i],bps_y[i],&ra,&dec);
      /*this is residual in arcseconds (squared)*/
      c+=3600*3600*((ra-bps_ra[i])*(ra-bps_ra[i])+(dec-bps_dec[i])*(dec-bps_dec[i]));
      used+=1;
    }
  }
  *value = c/used;
  return(0);
}

/*calculates the average residual of points to a WCS fit*/
/*here the non-linear terms are thrown at the minimizer */
/*with the linear transformation done as the first step*/
/*used by minimrq program*/

int calcerrlin(int npar, double *par, int *usepar, double *value,
	int doderiv, double *deriv, double *curve)
{
  int i, k;
  double used=1,chi=0,ra,dec,ra1,dec1;
  double a,b,c,d,e,f;
  void transform();

  /*update the wcs structure with values*/
  wcs->cd[0]=par[0]  ; /*cd1_1*/
  wcs->cd[1]=par[1]  ; /*cd1_2*/
  wcs->cd[2]=par[2]  ; /*cd2_1*/
  wcs->cd[3]=par[3]  ; /*cd2_2*/
  wcs->crpix[0]=par[4]  ; /*crpix1*/
  wcs->crpix[1]=par[5]  ;/*crpix2*/  
  wcs->pv[3].value=par[6]  ;/*pv_3*/   
  wcs->pv[5].value=par[7]  ;/*pv_5*/   
  wcs->crval[0]=par[8];
  wcs->crval[1]=par[9];

  wcsset(wcs);

  /*build up info for linear fit from all point*/
  for(i=0;i<bps_npts;i++) {
    cir_radectoxieta(wcs,bps_ra[i],bps_dec[i],&bps_xi[i],&bps_eta[i]);
    bps_xi[i]=bps_xi[i]*180/3.1415926535897*3600;
    bps_eta[i]=bps_eta[i]*180/3.1415926535897*3600;
  }
  /*fnd a linear solution*/
  transform(bps_npts,bps_xi,bps_eta,bps_x,bps_y,&c,&a,&b,&f,&d,&e);

  par[0] = a/3600;
  par[1] = b/3600;
  par[2] = d/3600;
  par[3]= e/3600;
  par[4] = (e*c - b*f)/(d*b - e*a);
  par[5] = (a*f - d*c)/(d*b - e*a);

  /*put into wcs structure*/
  wcs->cd[0]=par[0]  ; /*cd1_1*/
  wcs->cd[1]=par[1]  ; /*cd1_2*/
  wcs->cd[2]=par[2]  ; /*cd2_1*/
  wcs->cd[3]=par[3]  ; /*cd2_2*/
  wcs->crpix[0]=par[4]  ; /*crpix1*/
  wcs->crpix[1]=par[5]  ;/*crpix2*/  

  wcsset(wcs);
  
  
  /*find residuals of all stars*/
  for(i=0;i<bps_npts;i++) {
    if (bps_usept[i]) {
      cir_xytoradec(wcs,bps_x[i],bps_y[i],&ra,&dec);
      /*this is residual in arcseconds (squared)*/
      chi+=3600*3600*((ra-bps_ra[i])*(ra-bps_ra[i])+(dec-bps_dec[i])*(dec-bps_dec[i]));
      used+=1;
    }
  }
  *value = chi/used;
  return(0);
}

/* given a fit calculate the residual as a 6 term polynomial*/
/* used by minimrq in linear mode to find 6 terms*/

int calcresid(int npar, double *par, int *usepar, double *value,
	int doderiv, double *deriv, double *curve)
{
  int i, k;
  double used=0,chi=0,delta;
  double k00,k10,k20,k11,k01,k02;

  k00=par[0];
  k10=par[1];
  k20=par[2];
  k11=par[3];
  k01=par[4];
  k02=par[5];
  
  for(i=0;i<bps_npts;i++) {
    if (bps_usept[i]==1) {
      delta=k00+k10*bps_x[i]+k20*bps_x[i]*bps_x[i]+k11*bps_x[i]*bps_y[i]+k01*bps_y[i]+k02*bps_y[i]*bps_y[i];
      used++;
      if (calcResXflag==1) {
	chi+=(delta-residra[i])*(delta-residra[i]);
      }
      else {
	chi+=(delta-residdec[i])*(delta-residdec[i]);
      }
    }
  }
  *value = chi/used;
  return(0);
}

/*apply residual fit to the list of residuals a global yuk! (I am sorry) variable*/
int applyresid(int npar, double *par)
{
  int i, k;
  double delta;
  double k00,k10,k20,k11,k01,k02;
  
  k00=par[0];
  k10=par[1];
  k20=par[2];
  k11=par[3];
  k01=par[4];
  k02=par[5]; 

  for(i=0;i<bps_npts;i++) {
    delta=k00+k10*bps_x[i]+k20*bps_x[i]*bps_x[i]+k11*bps_x[i]*bps_y[i]+k01*bps_y[i]+k02*bps_y[i]*bps_y[i];
    if (calcResXflag==1) {
      residra[i]-=(delta);
    }
    else {
      residdec[i]-=(delta);
    }
  }
  return(0);
}

/*given a list of x and y positions, apply the residual fit and the WCS fit to the positions*/
/*and fill up ra and dec with their positions*/
int xy2wcs(int npar, double *parx, double *pary, int nlist, double *x, double *y, double *ra, double *dec)
{
  int i, k;
  double delta;
  double k00,k10,k20,k11,k01,k02;

  for (i=0;i<nlist;i++) {
    k00=parx[0];
    k10=parx[1];
    k20=parx[2];
    k11=parx[3];
    k01=parx[4];
    k02=parx[5]; 
    
    cir_xytoradec(wcs,x[i],y[i],&ra[i],&dec[i]);

      
    
    delta=k00+k10*x[i]+k20*x[i]*x[i]+k11*x[i]*y[i]+k01*y[i]+k02*y[i]*y[i];



    ra[i]-=delta/3600.;
    
    k00=pary[0];
    k10=pary[1];
    k20=pary[2];
    k11=pary[3];
    k01=pary[4];
    k02=pary[5]; 
    
    delta=k00+k10*x[i]+k20*x[i]*x[i]+k11*x[i]*y[i]+k01*y[i]+k02*y[i]*y[i];
    dec[i]-=delta/3600. ;

    if (verbosity > 2) {
      fprintf(stderr,"XY2WCS: %9.3f %9.3f %15.8f %15.8f %15.8f\n",x[i],y[i],ra[i],dec[i],delta);
    }

  }
}

void syntax(char *name) 
{
  printf("\nbrian_fitwcs: image posfile\n");
  printf("\twhere image is a FITS image with a ROUGH WCS system of (ZPN) or (TAN)\n");
  printf("\tposfile is a file with x y RA(deg) DEC(deg) on each line\n");
  printf("\n");
  printf("Optional arguments are (default value)\n");
  printf("\t-sigclip (3)\n");
  printf("\t-outfile (posfile.wcs)\n");
  printf("\t-floormatch (0.1) in arcsec is the lowest level a sig clip with prune objects\n");
  printf("\t-maxit (3) maximum sigclip iterations\n");
  printf("\t-usegaussig use gaussian estimate of error for sig-clipping instead of median\n");
  printf("\t-fitresid do a 6-term polynomial residual fit after doing ZPN Fit output to posfile.resid file\n");
  printf("\t-fitpv Fit for the pv3 term in the ZPN - this is a non-linear fit\n");
  printf("\t-fitpv5 Fit for the pv3 and pv5 terms in the ZPN - this is a non-linear fit\n");
  printf("\t-fitcrval Fit for the crval terms - Probably doesn't work\n");
  printf("\t-outputgrid endpixx endpixy dx dy - 4 arguments required output to posfile.grid file");
  printf("\t-outfile (posfile.fit) for wcs fit output");
  printf("\t-verbose (0)level of verbosity 0-3\n");
  return;
}


/*calculates the residuals of points to a WCS fit and output this to a list.*/
void residcalc(double *resxi, double *reseta)
{
  int i, k;
  double used=1,c=0,ra,dec,ra1,dec1;

  /*find residuals of all stars currently in list*/
  for(i=0;i<bps_npts;i++) {
    cir_xytoradec(wcs,bps_x[i],bps_y[i],&ra,&dec);
    if (verbosity > 2) {
      fprintf(stderr,"RESIDCALC: %f %f %f %f\n",bps_x[i],bps_y[i],ra,dec);
    }      
    resxi[i]=3600*(ra-bps_ra[i]);
    reseta[i]=3600*(dec-bps_dec[i]);
  }
}

/*returns the error determined from the median of a list*/
void medianerr(double *resxi, double *reseta,double *xierr, double *etaerr)
{
  int i,j,k=0;
  float x[10000];
  
  /*start with resxi*/
  /*find ones we are are using*/
  for(i=0;i<bps_npts;i++) {
    if (bps_usept[i]) { 
      x[k++]=fabs(resxi[i]);
    }
  }
  
    /*now select 68%*/
  /*j is the element we select out of the k that remain*/
  j=(k*.683);
  select(j,k,x);
  *xierr=(double)x[j];


  /*reseta*/
  /*find ones we are are using*/
  k=0;
  for(i=0;i<bps_npts;i++) {
    if (bps_usept[i]) { 
      x[k++]=fabs(reseta[i]);
    }
  }
  
  /*now select 68%*/
  /*j is the element we select out of the k that remain*/
  j=(k*.683);
  select(j,k,x);
  *etaerr=(double)x[j];
}
/*returns the RMS error*/ 
void gausserr(double *resxi, double *reseta,double *xierr, double *etaerr)
{
  int i,j,k=0;
  double sumxx=0.,sumx=0.,nx=0.,meanx;
  double sumyy=0.,sumy=0.,ny=0.,meany;
  
  for(i=0;i<bps_npts;i++) {
    if (bps_usept[i]) { 
      sumx+=resxi[i];
      nx++;
      sumy+=reseta[i];
      ny++;
    }
  }
  
  meanx=sumx/nx;
  meany=sumy/ny;

  for(i=0;i<bps_npts;i++) {
    if (bps_usept[i]) { 
      sumxx+=(resxi[i]-meanx)*(resxi[i]-meanx);
      sumyy+=(reseta[i]-meany)*(reseta[i]-meany);
    }
  }
  *xierr=sqrt(sumxx/nx);
  *etaerr=sqrt(sumyy/ny);
}

/*exclude objects in a list based on a sigma clip parameter and floor and error*/

int excludeobjs (double *resxi, double *reseta,double xierr, double etaerr,double sig, double floor)
{
  int i,j,k=0;
  
  for(i=0;i<bps_npts;i++) {
    if ((fabs(resxi[i]) > sig*xierr && fabs(resxi[i]) > floor) || (fabs(reseta[i]) > sig*etaerr && fabs(reseta[i]) > floor)) {
      bps_usept[i]=0;
      k++;
    }
    else {
      bps_usept[i]=1;
    }
    if (verbosity>1 || (verbosity > 0 && bps_usept[i]==0)) {
      fprintf (stderr,"EXLUDEOBJ: %d %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n",bps_usept[i],resxi[i],reseta[i],resxi[i]/xierr,reseta[i]/etaerr,xierr,etaerr);
    }
  }
  return(k);
}

/*output a WCS fit to a file*/
void printobj(double *resxi, double *reseta,char *OUTPUT)
{
  int i,j,k=0;
  FILE *fd;

  if (strcmp(OUTPUT,"stderr")!=0) {
    
    if (verbosity > 0) {
      fprintf(stderr,"PRINTOBJ: opened file %s for output\n",OUTPUT);
    }
    
    if ((fd = fopen(OUTPUT,"w")) == NULL) {
      fprintf(stderr,"BRIAN_FITWCS.PRINTOBJ Cannot open %s for writing\n",OUTPUT);
      exit(-1);
    }
    
    
    for(i=0;i<bps_npts;i++) {
      fprintf(fd,"%12.4f %12.4f %12.7f %12.7f %7.3f %7.3f %d\n",bps_x[i],bps_y[i],bps_ra[i],bps_dec[i],resxi[i],reseta[i],bps_usept[i]);
    }
    close(fd);
  }
  
  else {
    for(i=0;i<bps_npts;i++) {
      fprintf(stderr,"ITERATION: %12.4f %12.4f %12.7f %12.7f %7.3f %7.3f %d\n",bps_x[i],bps_y[i],bps_ra[i],bps_dec[i],resxi[i],reseta[i],bps_usept[i]);
    }
  } 
}


/* calculate the*/
/* a linear transform of the form */
/* x= ax +bx*X + cx*Y */
/* y= ay +by*X + cy*Y */

void transform(nstar,x2,y2,x1,y1,ax,bx,cx,ay,by,cy)
     int nstar;
     double x2[],y2[],x1[],y1[],*ax,*bx,*cx,*ay,*by,*cy;
{

  int i;
  double sum=0.,sumx=0.,sumy=0.,sumxy=0.,sumu=0.,sumv=0.,sumx2=0.,sumy2=0.;
  double sumux=0.,sumuy=0.,sumvx=0.,sumvy=0.;
  double delta,sumdx,sumdy,sumabsdx,sumabsdy;
  double x1p,y1p,dx,dy;
  for (i=1; i<=nstar;i++) {
    if (bps_usept[i]==1) {
      sum ++;
      sumx += x1[i];
      sumy += y1[i];
      sumu += x2[i];
      sumv += y2[i];
      sumx2 += x1[i]*x1[i];
      sumxy += x1[i]*y1[i];
      sumy2 += y1[i]*y1[i];
      sumux += x1[i]*x2[i];
      sumvx += x1[i]*y2[i];
      sumuy += y1[i]*x2[i];
      sumvy += y1[i]*y2[i];
    }
  }
  delta = sum*(sumx2*sumy2 - sumxy*sumxy) + sumx*(sumxy*sumy - sumx*sumy2) + sumy*(sumx*sumxy-sumx2*sumy);
  *ax = (sumu*(sumx2*sumy2 - sumxy*sumxy) + sumx*(sumxy*sumuy - sumux*sumy2) + sumy*(sumux*sumxy - sumx2*sumuy))/delta;
  *bx = (sum*(sumux*sumy2 - sumxy*sumuy) + sumu*(sumxy*sumy - sumx*sumy2) + sumy*(sumuy*sumx - sumy*sumux))/delta;
  *cx = (sum*(sumx2*sumuy - sumxy*sumux) + sumx*(sumux*sumy - sumx*sumuy) + sumu*(sumx*sumxy - sumx2*sumy))/delta;
  *ay = (sumv*(sumx2*sumy2 - sumxy*sumxy) + sumx*(sumxy*sumvy - sumvx*sumy2) + sumy*(sumvx*sumxy - sumx2*sumvy))/delta;
  *by = (sum*(sumvx*sumy2 - sumxy*sumvy) + sumv*(sumxy*sumy - sumx*sumy2) + sumy*(sumvy*sumx - sumy*sumvx))/delta;
  *cy = (sum*(sumx2*sumvy - sumxy*sumvx) + sumx*(sumvx*sumy - sumx*sumvy) + sumv*(sumx*sumxy - sumx2*sumy))/delta;
  
}


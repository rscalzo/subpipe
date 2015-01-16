#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <cirdr.h>
#include <cir_wcssubs.h>

#define RAD2DEG 180.0/M_PI
#define TWOPI 2*M_PI
#define INITALLOC 64

static double ra,dec;
static double crval1,crval2,crpix1,crpix2,cd1_1,cd1_2,cd2_1,cd2_2,projp3;

static struct wcsprm *wcs = NULL;
static char cvsid[] = "$Id: cir_platesol.c,v 1.9 2004/09/07 14:18:54 jim Exp $";
static char errmsg[120];

main (argc,argv)
     int argc;
     char *argv[];
{
  
  int i,nalloc,npts,nrej,status,niter,ngood,nc2,retval;
  FILE *fd;
  float averr;
  double dxi,deta,daverr,xi,eta,r1,r2,d1,d2;
  double a,b,c,d,e,f,xifit,etafit,crpix1,crpix2,cd[4];
  double newcrval1,newcrval2;
  char msg[BUFSIZ],v1[16],v2[16],v3[16],v4[16],v5[16],v6[16],image[100],posfile[100],line[150];
  
  sprintf(image,argv[1]);
  sprintf(posfile,argv[2]);


  /* Open the fits image and read the current WCS info into WCS structure */
  
  cir_stamp(image,cvsid);
  retval = cir_wcsopen(image,&wcs,msg);
  if (retval != CIR_OK) {
    printf("PLATESOL: Couldn't read file input WCS %s -- %s\n", image,msg);
    exit(-1);
  }

  wcs->crpix[0]=1.0;

    if (retval != CIR_OK) {
        (void)sprintf(errmsg,"IMSTACK: Unable to create new WCS == %s\n",msg);
        return(CIR_FATAL);
    }

    /* Open the position file and read the relevant data. */

  if ((fd = fopen(posfile,"r")) == NULL) {
    printf("Brian_xi_eta: No input position file");
    exit(-1);
  }
  npts = 0;

  while(feof(fd)==0) {
    fgets(line,150, fd);
    sscanf(line,"%s %s %s",v1,v2,v3);
    
    ra = atof(v2);
    dec= atof(v3);
    cir_radectoxieta(wcs,ra,dec,&xi,&eta);
 
    xi *=180/3.1415926535897*3600;
    eta *=180/3.1415926535897*3600;
    
    printf("%15s\t%9.4f\t%9.4f\t%15s\t%12.9f\t%12.9f\n",v1,xi,eta,v1,ra,dec);
    npts++;
  }
}


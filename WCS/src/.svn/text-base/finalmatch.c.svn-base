#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#define MOBJ 4000
#define MTRI 25005
#define NPAR 7
#define MERRMAX 0.1
#define ERRMAX 1.0
#define ROTERR 0.02


main (argc,argv)
     int argc;
     char *argv[];
{
  
  FILE *stream1,*stream2,*stream3;

  float *x1, *y1, *x2,*y2;

  float *vector();

  int nobj,nstart1,nstart2,i,j,k,l,m,ii,jj,kk,ll,mm,itype,nstar1,nstar2,q,r,rmax,firsttime;
  int jbest,nm,nf,irepeat,anf,igotavg,iend,id,imtype,nsum1,nsum2,jlow,jhigh;
  float errmax,sepmin,tolerance,rtolerance,sigcut,rsummin; 
  float err,r12,r13,r23,xhunt,tol,delta,temp;
  float sum,sumrs,xtest,ytest,sumr2,sumr1,avgdiff,gamma;
  float ax,bx,cx,ay,by,cy,sumdx,sumdy,sumdx2,sumdy2;
  float sdx,sdy,xcut,ycut,XC,YC,xc,yc,TILT,AMINOR,AMAJOR,AREA,FMAG;
  float xx,yy,bmag,scaleratio,srtolerance,ratio,x,y,dummy;
  float ratiomap[10],starpar[NPAR],staravg[NPAR];
  float AMINOR1,AMAJOR1,TILT1;
  char ians,line[150],file1[80],file2[80],file3[80],file4[80],dummyc[1000];
  char dummy1[MOBJ][30],dummy2[MOBJ][30],dummy3[MOBJ][30];
  int ns1,ns2,nstar,nbr,nbl,nbb,nbt,gotit=0,skiprow,good=0;
  float mag2,emag2,FERR,xa1,ya1,xa2,ya2,xc1,xc2,yc1,yc2;
  int bestj,wcsout;
  float bestdist,limit,dist;


  void elliptint(),transform();

  if ( argc != 14 ) {
     printf( "Usage: %s templist objlist ax bx cx ay by cy closerthan totalmatchlimit skiprow wcsout\n", argv[0] );
     exit( 1 );
  }

  
  x1=vector(0,MOBJ);
  y1=vector(0,MOBJ);
  x2=vector(0,MOBJ);
  y2=vector(0,MOBJ);

  ns2=10000;
  

  sprintf(file1,argv[1]);
  sprintf(file2,argv[2]);
  sprintf(file3,argv[3]);
  bx = atof( argv[4] );
  cx = atof( argv[5] );
  ax = atof( argv[6] );
  by = atof( argv[7] );
  cy = atof( argv[8] );
  ay = atof( argv[9] );
  limit = atof( argv[10]);
  nobj = atof( argv[11] );
  skiprow = atof( argv[12] );
  wcsout = atof( argv[13] );

  if ((stream2 = fopen(file2,"r"))== NULL) {
    printf("\nFile %s Not found...exiting before I dump the core\n", file2);
    exit(0);
  }
  
  if ((stream1 = fopen(file1,"r"))== NULL) {
    printf("File %s Not found - exiting...\n",file1);
    exit(0);
  } 

  if ((stream3 = fopen(file3,"w"))== NULL) {
    printf("File %s cannot be openned...exiting\n",file1);
    exit(0);
  } 
  
  
  /* load up first file */ 
  for (i=1; feof(stream1) == 0 && i < nobj && i < MOBJ; i++) {
    fgets(line,150, stream1); 
    if (skiprow ==0) {
      if(strncmp(line, "#", 1) != 0) {
	if (sscanf(line,"%f %f %s %s %s",&x1[i],&y1[i],dummy1[i],dummy2[i],dummy3[i]) <1) i--;
      }
      else i--;
    }
    else {
      if(strncmp(line, "#", 1) != 0) {
	if (sscanf(line,"%s %f %f %s %s %s",dummyc,&x1[i],&y1[i],dummy1[i],dummy2[i],dummy3[i]) <1) i--;
      }
      else i--;
    }
  }
  nstar1 = i-1;
  if (nstar1 < 3) {
     exit( 1 );
  }
  
  /* load up 2nd file */
  for (i=1; feof(stream2) == 0 && i < nobj && i < MOBJ; i++) {
     fgets(line,150, stream2); 
         if (skiprow ==0) {
	   if(strncmp(line, "#", 1) != 0) {
	     if (sscanf(line,"%f %f",&x2[i],&y2[i]) < 1 ) i--;
	   }
	   else i--;
	 }
	 else {
	   if(strncmp(line, "#", 1) != 0) {
	     if (sscanf(line,"%s %f %f",dummyc,&x2[i],&y2[i]) < 1 ) i--;
	   }
	   else i--;
	 }
  }
  nstar2 = i-1;
  
  fclose(stream1);
  fclose(stream2);

  for ( i = 1; i <= nstar1; i++ ) {
    bestdist = 999.99;
    for ( j = 1; j <= nstar2; j++ ) {
      x = x2[j] * bx + y2[j] * cx + ax;
      y = x2[j] * by + y2[j] * cy + ay;
      dist = sqrt( ( x1[i] - x ) * ( x1[i] - x ) + ( y1[i] - y ) * ( y1[i] - y ) );
      if ( dist < bestdist ) {
	bestdist = dist;
	bestj = j;
      }
    }
    if ( bestdist < limit ) {
      if (wcsout==0) {
	fprintf(stream3,"%8.2f %8.2f %8.2f %8.2f %8.2f %s %s %s\n",x1[i],y1[i],x2[bestj],y2[bestj],bestdist,dummy1[i],dummy2[i],dummy3[i]);
      }
      else {
	fprintf(stream3,"%8.2f %8.2f %20s %20s %6.2f\n",x2[bestj],y2[bestj],dummy2[i],dummy3[i],bestdist);
      }
    }
  }
  fclose(stream3);
}

 


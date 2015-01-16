#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define MOBJ 5001
#define MTRI 2000000

/*Brian Schmidt Apr 2003 re-write of code from 1990-2002*/

main (argc,argv)
     int argc;
     char *argv[];
{
  
  FILE *stream1,*stream2,*stream3,*stream4;

  float *x1, *y1, *amag1, *rsum1, *x2,*y2, *amag2,*rsum2,*xi,*eta;
  float *xm1, *ym1,*xm2,*ym2,*xf1,*yf1,*xf2,*yf2,*xh1,*yh1,*xh2,*yh2,*amagdiff;
  float *ratiolist,*ratsky,*sky1,*sky2,*diffx,*diffy,*xpos1,*ypos1,*mag1,*emag1;
  float *alist,*blist,*clist,*dlist,*elist,*flist;
  float **r2,**r1,**rabs1,**rabs2;
  int *index1, *indexr1,  *index2,*indexr2,*index, **ijk1,**ijk2;
  float *vector(),**matrix();
  int   *ivector(),**imatrix();
  void syntax();


  int nobj=-1,nstart1=1,nstart2=1,i,j,k,l,m,ii,jj,kk,ll,mm,itype,nstar1,nstar2,q,r=0,rmax,firsttime;
  int iii,jjj,kkk,rnew;
  int jbest,nm,nf,irepeat,anf,igotavg,iend,id,imtype,nsum1,nsum2,jlow,jhigh;
  int STARMATCHES,skipcol=1; /*added BPS 29mar06 unitialised*/
  float abest=0.,bbest=0.,cbest=0.,dbest=0.,ebest=0.,fbest=0.,amiss,cbmiss,totalmiss=0.1;
  float atest,btest,ctest,dtest,etest,ftest,sumtest,scaletmp;
  float errmax,sepmin,tolerance,rtolerance,sigcut,rsummin; 
  float err,xhunt,tol,delta,temp;
  float sum,sumrs,xtest,ytest,sumr2,sumr1,avgdiff,gamma;
  float ax,bx,cx,ay,by,cy,sumdx,sumdy,sumdx2,sumdy2;
  float xcut,ycut,XC,YC,xc,yc,TILT,AMINOR,AMAJOR,AREA,FMAG;
  float xx,yy,bmag,scaleratio,srtolerance,ratio,x,y,dummy;
  float AMINOR1,AMAJOR1,TILT1,ROTERR;
  char ians,line[150],file1[80],file2[80],file3[80],file4[80],dummyc[1000],outfile[100];
  int ns1,ns2,nstar,nbr,nbl,nbb,nbt,gotit=0,diagnose,maxiterations=5,iterations;
  float mag2,emag2,FERR,xa1,ya1,xa2,ya2,xc1,xc2,yc1,yc2;

  double r12,r13,r23,R12,R13,sdx,sdy;
  double a,b,c,d,e,f,side1,side2,side3,longside1,longside2,scale,a1,a2,a3,squaretol=3,rotation=0;
  double OKx, OKy,floormatch=0.001,sigmatch=2.5;
  float **xfm1,**yfm1,**xfm2,**yfm2;

  int fail,solution,*matches,bestmatch=0,needmatches,matchid,nstars,neededmatches=5,outids,iwaussian;
  int nstarsold,xieta;

  int js=0,jumpstartref1[100],jumpstartref2[100],jstmp=0;
  char jumpstartid1[100][100],jumpstartid2[100][100];
  
  void transform();

  
  x1=vector(0,MOBJ);
  y1=vector(0,MOBJ);
  xi=vector(0,MOBJ);
  eta=vector(0,MOBJ);
  amag1=vector(0,MOBJ);
  rsum1=vector(0,MTRI);
  x2=vector(0,MOBJ);
  y2=vector(0,MOBJ);
  amag2=vector(0,MOBJ);
  rsum2=vector(0,MTRI);
  xm1=vector(0,MTRI);
  ym1=vector(0,MTRI);
  xm2=vector(0,MTRI);
  ym2=vector(0,MTRI);
  xf1=vector(0,MTRI);
  yf1=vector(0,MTRI);
  xf2=vector(0,MTRI);
  yf2=vector(0,MTRI);
  xh1=vector(0,MTRI);
  yh1=vector(0,MTRI);
  xh2=vector(0,MTRI);
  yh2=vector(0,MTRI);
  xpos1=vector(0,MTRI);
  ypos1=vector(0,MTRI);
  mag1=vector(0,MTRI);
  emag1=vector(0,MTRI);
  amagdiff=vector(0,MTRI);
  ratiolist=vector(0,MTRI);
  alist=vector(0,MTRI);
  blist=vector(0,MTRI);
  clist=vector(0,MTRI);
  dlist=vector(0,MTRI);
  elist=vector(0,MTRI);
  flist=vector(0,MTRI);
  ratsky=vector(0,MOBJ);
  sky1=vector(0,MOBJ);
  sky2=vector(0,MOBJ);
  diffx=vector(0,MOBJ);
  diffy=vector(0,MOBJ);
  index1=ivector(0,MOBJ);
  indexr1=ivector(0,MTRI);
  index2=ivector(0,MOBJ);
  indexr2=ivector(0,MTRI);
  index=ivector(0,MOBJ);
  r2=matrix(0,3,0,MTRI);
  r1=matrix(0,3,0,MTRI);
  rabs1=matrix(0,3,0,MTRI);
  rabs2=matrix(0,3,0,MTRI);
  ijk1=imatrix(0,3,0,MTRI);
  ijk2=imatrix(0,3,0,MTRI);
  xfm1=matrix(0,3,0,MTRI);
  xfm2=matrix(0,3,0,MTRI);
  yfm1=matrix(0,3,0,MTRI);
  yfm2=matrix(0,3,0,MTRI);
  matches=ivector(0,MTRI);

  if (argc <3) {
    syntax(argv[0]);
    exit(0);
  }

  sprintf(file1,argv[1]);
  sprintf(file2,argv[2]);

  /* default values for everything*/

  diagnose=0;
  scaleratio=0.0;
  srtolerance=0.01; 
  STARMATCHES=5;
  /*  nobj=50;*/
  iwaussian=0;       /* switch to use waussian output as input */
  outids=0;          /*append matched ids in the output as cols5 & 6*/
  amiss=0.05;
  xieta=0; /* 1 means first file is a .xieta file with RA and DEC in cols 5 and 6*/ 
  cbmiss=0.05;
  totalmiss=0.1;
  ROTERR=3.;
  rotation=0.;
  squaretol=3;
  sepmin = 8.0;      /* Min Separation in pixels for Triangle */
  tolerance = 0.001;  /* Tolerance in Similarity of Triangle */
  rtolerance = 0.001;   /* other Tolerance in Similarity of Triangle */
  rsummin = 1.10;     /* Min obliqueness of triangle */
  strcpy(outfile,"matchstar.out");

  if (argc > 3) { /* parse optional command line arguments */
    for (i=3;i<argc;i++) {

      /*diagnose option*/
      if ( strcasecmp(argv[i],"-diagnose") == 0 ) {
	if ( ++i >= argc ) {
	  fprintf(stderr,"-diagnose option requires argument: 0 1 2 3 (0 low level 3 everything and the kitchen sink)\n");
	  syntax(argv[0]);
	  exit(0);
	}
	diagnose = atoi(argv[i]);
	if (diagnose < 0 || diagnose > 3 ) {
	  fprintf(stderr,"invalid diagnose argument (%d)\n",diagnose);
	  syntax(argv[0]);
	  exit(0);
	}
	continue;
      }
      
      /*scale option*/
      if ( strcasecmp(argv[i],"-scale") == 0 ) {
	if ( ++i >= argc ) {
	  fprintf(stderr,"-scale option requires argument: scale ratio of image / template (0 means figure it out)\n");
	  syntax(argv[0]);
	  exit(0);
	}
	scaleratio = atof(argv[i]);
	if (scaleratio < 0 ) {
	  fprintf(stderr,"scale cannot be less than 0 [%f]\n",scaleratio);
	  syntax(argv[0]);
	  exit(0);
	}
	continue;
      }


            /*scale tolerance option*/
      if ( strcasecmp(argv[i],"-tol") == 0 ) {
	if ( ++i >= argc ) {
	  fprintf(stderr,"-tol option requires argument: tolerance in specified scale ratio of image / template\n");
	  syntax(argv[0]);
	  exit(0);
	}
	srtolerance = atof(argv[i]);
	if (srtolerance < 0 ) {
	  fprintf(stderr,"scale tolerance (-tol) cannot be less than 0 [%f]\n",srtolerance);
	  syntax(argv[0]);
	  exit(0);
	}
	continue;
      }

      /*no id in column 1 in star file*/
      if ( strcasecmp(argv[i],"-noid") == 0 ) {
	skipcol=0;
      	continue;
      }

      /*Would like id columns in the output*/
      if ( strcasecmp(argv[i],"-outids") == 0 ) {
	outids=1;
      	continue;
      }

      /*First file is a xieta file, and we then output RA DEC */
      if ( strcasecmp(argv[i],"-xieta") == 0 ) {
	xieta=1;
      	continue;
      }

      /*Use waussian input*/
      if ( strcasecmp(argv[i],"-wauss") == 0 ) {
	iwaussian=1;
      	continue;
      }


      /*how many triangles is a success*/
      if ( strcasecmp(argv[i],"-ntri") == 0 ) {
	if ( ++i >= argc ) {
	  fprintf(stderr,"-ntri option requires argument: total number of triangles match to declare victory\n");
	  syntax(argv[0]);
	  exit(0);
	}
	STARMATCHES = atoi(argv[i]);
	if (STARMATCHES < 1 ) {
	  fprintf(stderr,"-ntri  cannot be less than 1 [%d]\n",STARMATCHES);
	  syntax(argv[0]);
	  exit(0);
	}
	continue;
      }

      /*how many stars do we match at a time*/
      if ( strcasecmp(argv[i],"-nmatch") == 0 ) {
	if ( ++i >= argc ) {
	  fprintf(stderr,"-nmatch option requires argument: total number of stars to match at a time\n");
	  syntax(argv[0]);
	  exit(0);
	}
	nobj = atoi(argv[i]);
	if (nobj < 5 ) {
	  fprintf(stderr,"-nmatch  cannot be less than 5 [%d]\n",nobj);
	  syntax(argv[0]);
	  exit(0);
	}
	continue;
      }

      /*how close do the offsets have to be (as percentage)*/
      if ( strcasecmp(argv[i],"-amiss") == 0 ) {
	if ( ++i >= argc ) {
	  fprintf(stderr,"-amiss option requires argument: fractional miss allowed in constants\n");
	  syntax(argv[0]);
	  exit(0);
	}
	amiss = atof(argv[i]);
	if (amiss < 0.000001 ) {
	  fprintf(stderr,"-amiss  must be more than 0 [%f]\n",amiss);
	  syntax(argv[0]);
	  exit(0);
	}
	continue;
      }


      /*how close do the rotation matrix have to be (as percentage)*/
      if ( strcasecmp(argv[i],"-cbmiss") == 0 ) {
	if ( ++i >= argc ) {
	  fprintf(stderr,"-cbmiss option requires argument: fractional miss allowed in rot matrix\n");
	  syntax(argv[0]);
	  exit(0);
	}
	cbmiss = atof(argv[i]);
	if (cbmiss < 0.000001 ) {
	  fprintf(stderr,"-cbmiss  must be more than 0 [%f]\n",cbmiss);
	  syntax(argv[0]);
	  exit(0);
	}
	continue;
      }

      /*Allow Rotation ?*/
      if ( strcasecmp(argv[i],"-roterr") == 0 ) {
	if ( ++i >= argc ) {
	  fprintf(stderr,"-roterr option requires argument: angular rotation allowed from 0,90,180,270\n");
	  syntax(argv[0]);
	  exit(0);
	}
	ROTERR = atof(argv[i]);
	if (ROTERR < 0.000001) {
	  fprintf(stderr,"-roterr  must be more than 0 [%f]\n",ROTERR);
	  syntax(argv[0]);
	  exit(0);
	}
	continue;
      }

      /*total agreement in matrix allowed*/
      if ( strcasecmp(argv[i],"-totmiss") == 0 ) {
	if ( ++i >= argc ) {
	  fprintf(stderr,"-totmiss option requires argument: total sum of rotational matrix errors allowed\n");
	  syntax(argv[0]);
	  exit(0);
	}
	totalmiss = atof(argv[i]);
	if (totalmiss < 0.000001) {
	  fprintf(stderr,"-totmiss  must be more than 0 [%f]\n",totalmiss);
	  syntax(argv[0]);
	  exit(0);
	}
	continue;
      }


      /*how square does matrix have to be*/
      if ( strcasecmp(argv[i],"-squaretol") == 0 ) {
	if ( ++i >= argc ) {
	  fprintf(stderr,"-squaretol option requires argument: in degrees how far away from square allowed\n");
	  syntax(argv[0]);
	  exit(0);
	}
	squaretol = atof(argv[i]);
	if (squaretol < 0.000001) {
	  fprintf(stderr,"-squaretol  must be more than 0 [%f]\n",squaretol);
	  syntax(argv[0]);
	  exit(0);
	}
	continue;
      }

      /*expected angle of the rotation*/
      if ( strcasecmp(argv[i],"-rotation") == 0 ) {
	if ( ++i >= argc ) {
	  fprintf(stderr,"-rotation option requires argument: in degrees relative rotation of template and image\n");
	  syntax(argv[0]);
	  exit(0);
	}
	rotation = atof(argv[i]);
	continue;
      }

      /* Min Separation in pixels for Triangle */

      if ( strcasecmp(argv[i],"-sepmin") == 0 ) {
	if ( ++i >= argc ) {
	  fprintf(stderr,"-sepmin option requires argument: in pixels what is the minimum size of a triangle side\n");
	  syntax(argv[0]);
	  exit(0);
	}
	sepmin = atof(argv[i]);
	if (squaretol < 0.000001) {
	  fprintf(stderr,"-sepmin  must be more than 0 [%f]\n",sepmin);
	  syntax(argv[0]);
	  exit(0);
	}
	continue;
      }


      /* Sigclip value */

      if ( strcasecmp(argv[i],"-sigclip") == 0 ) {
	if ( ++i >= argc ) {
	  fprintf(stderr,"-sigclip option requires argument: Sigma value for iteratively removing stars from fit\n");
	  syntax(argv[0]);
	  exit(0);
	}
	sigmatch = atof(argv[i]);
	if (sigmatch < 1.00000001) {
	  fprintf(stderr,"-sigmatch  must be more than 1 [%f]\n",sigmatch);
	  syntax(argv[0]);
	  exit(0);
	}
	continue;
      }

      /* maxnumber of iterations value */

      if ( strcasecmp(argv[i],"-maxit") == 0 ) {
	if ( ++i >= argc ) {
	  fprintf(stderr,"-maxit option requires argument: max number of iterations for iteratively removing stars from fit\n");
	  syntax(argv[0]);
	  exit(0);
	}
	maxiterations = atoi(argv[i]);
	if (maxiterations < 1 || maxiterations > 100) {
	  fprintf(stderr,"maxit  must be between 1 an 100 [%d]\n",maxiterations);
	  syntax(argv[0]);
	  exit(0);
	}
	continue;
      }

      /* floormatch value */

      if ( strcasecmp(argv[i],"-floormatch") == 0 ) {
	if ( ++i >= argc ) {
	  fprintf(stderr,"-floormatch option requires argument: Minimum value for removing objects\n");
	  syntax(argv[0]);
	  exit(0);
	}
	floormatch = atof(argv[i]);
	if (floormatch < 0) {
	  fprintf(stderr,"-floormatch  must be more than 0 [%f]\n",sigmatch);
	  syntax(argv[0]);
	  exit(0);
	}
	continue;
      }


      /* Tolerance in Similarity of Triangle */
      if ( strcasecmp(argv[i],"-tolerance") == 0 ) {
	if ( ++i >= argc ) {
	  fprintf(stderr,"-tolerance option requires argument: fraction that triangles can deviate from each other\n");
	  syntax(argv[0]);
	  exit(0);
	}
	tolerance = atof(argv[i]);
	if (tolerance < 0.000001) {
	  fprintf(stderr,"-tolerance  must be more than 0 [%f]\n",tolerance);
	  syntax(argv[0]);
	  exit(0);
	}
	continue;
      }


     
      /* Min obliqueness of triangle */
      if ( strcasecmp(argv[i],"-rsummin") == 0 ) {
	if ( ++i >= argc ) {
	  fprintf(stderr,"-rsummin option requires argument: min obliqueness of triangle allowed (1-3)\n");
	  syntax(argv[0]);
	  exit(0);
	}
	rsummin = atof(argv[i]);
	if (rsummin < 1) {
	  fprintf(stderr,"-rsummin  must be more than 1 [%f]\n",rsummin);
	  syntax(argv[0]);
	  exit(0);
	}
	continue;
      }



           
      /* jumpstart the process*/
      if ( strcasecmp(argv[i],"-jumpstart") == 0 ) {
	if ( ++i >= argc ) {
	  fprintf(stderr,"-jumpstart requires two arguments: idfile1 idfile2)\n");
	  syntax(argv[0]);
	  exit(0);
	}
	strcpy(jumpstartid1[js],argv[i]);
	if ( ++i >= argc ) {
	  fprintf(stderr,"-jumpstart requires two arguments: idfile1 idfile2)\n");
	  syntax(argv[0]);
	  exit(0);
	}
	strcpy(jumpstartid2[js++],argv[i]);
	continue;
      }

      
      /*Outputfile */
      if ( strcasecmp(argv[i],"-outfile") == 0 ) {
	if ( ++i >= argc ) {
	  fprintf(stderr,"-outfile option requires argument: filename to output matched list\n");
	  syntax(argv[0]);
	  exit(0);
	}
	strcpy(outfile,argv[i]);
	continue;
      }

      /* unrecognized command line option */
      fprintf(stderr,"unrecognized command line option: %s\n", argv[i]);
      syntax(argv[0]);
      exit(0);
    }     

  }

  /* test if user wants output ids but doesn't have them */
  if (skipcol ==0 && outids ==1) {
    fprintf(stderr,"output ids requested with no input ids\n");
    exit(0);
  }


 /* test if user wants output ids but doesn't have them */
  if (skipcol ==0 && xieta ==1) {
    fprintf(stderr,"skipcol and xieta not compatable\n");
    exit(0);
  }

 /* test if user wants output ids but doesn't have them */
  if (iwaussian ==1 && xieta ==1) {
    fprintf(stderr,"-iwaussian and -xieta not compatable\n");
    exit(0);
  }


  if (nobj<5 && js <1) { /*not given, set the number of triangles to 50 if we are jumpstarting*/
    nobj=50;             /*roughly 10000 triangles*/
  }

  if (nobj<5 && js==1) { /*not given, using 1 star, set the number of triangles to 1000*/
    nobj=500;           /*we have N*N-1/2 triangles, so this is a big but tenable number*/
  }                     /*increase to as large as 2000 if you are having problems*/ 
  
  if (nobj<5 && js==2) { /*not given, using 1 star, set the number of triangles to 100000*/
    nobj=100000;         /*e.g., use all the objects since there are only N triangles*/ 
  }
  
  
  
 
  neededmatches=STARMATCHES;

  if ((stream2 = fopen(file2,"r"))== NULL) {
    printf("\nTemplate File Not found...exiting before I dump the core\n");
    exit(0);
  }
  
  if ((stream1 = fopen(file1,"r"))== NULL) {
    printf("File %s Not found - exiting...\n",file1);
    exit(0);
  } 
  
  for (j=0;j<js;j++) { /* initialise jumpstarts*/
    jumpstartref1[j]=-1;
  }


  /* load up first file */ 
  for (i=1; feof(stream1) == 0 && i < 5000; i++) {
    if (iwaussian) {
          fgets(line,200, stream1); 
    }
    else {
      fgets(line,150, stream1); 
    }
    if(strncmp(line, "#", 1) != 0) {
      if (skipcol == 0)  {
	sscanf(line,"%f %f %f",&x1[i],&y1[i],&amag1[i]);
      }
      else {
	if (xieta==1) {
	  sscanf(line,"%f %f %f %s %f %f",&x1[i],&y1[i],&amag1[i],dummyc,&xi[i],&eta[i]);
	}
	if (iwaussian) {
	  sscanf(line,"%s %f %f %s %s %s %s %s %s %s %s %f",dummyc,&x1[i],&y1[i],dummyc,dummyc,dummyc,dummyc,dummyc,dummyc,dummyc,dummyc,&amag1[i]);
	  if (amag1[i] > 0) {
	    amag1[i] = -2.5*log(amag1[i])/log(10);
	  }
	  else {
	    amag1[i]=99.0;	    
	  }
	}
	else {
	  sscanf(line,"%s %f %f %f",dummyc,&x1[i],&y1[i],&amag1[i]);
	}
	  
      }
      
      if (js > 0) { /*we are jumpstarting*/
	for (j=0;j<js;j++) {
	  if (strcasecmp(jumpstartid1[j],dummyc) == 0 ) {
	    jumpstartref1[j]=i;
	    if (diagnose) {
	      printf("jumpstart star %s given reference %d for %s\n",jumpstartid1[j],i,file1);
	    }
	  }
	}
      }
      if (amag1[i] ==0) amag1[i]=99.99;
      if (x1[i]==0 && y1[i]==0) i--;
    }
    else i--;
  }
  nstar1 = i-1;
  if (nstar1 < 3)  {
    printf("\n%s only has %d Star/s...Cannot proceed (need at least 3)\n",file1,nstar1);
    exit(-1);
  }

 
  indexx(nstar1,amag1,index1); /* sort by magnitudes  - give sort into index1 */

  fail=0;
  /*test that we got our jump stars in list and update with sorted value*/
  for (j=0;j<js;j++) {
    if (jumpstartref1[j] < 0) {
      printf("Didn't find jumpstart star %s in %s\n",jumpstartid1[j],file1);
      fail=1;
    }
    for (k=1;k<=nstar1;k++) {
      if (index1[k] == jumpstartref1[j]) {
	jumpstartref1[j]=k;
	k=nstar1+1;
      }
    }
  }
  if (fail) {
    exit(-1);
  }



  if (diagnose) {
    printf("read in %d stars from %s\n",nstar1,file1);
  }  
  


  for (j=0;j<js;j++) { /* initialise jumpstarts*/
    jumpstartref2[j]=-1;
  }
   
  /* load up 2nd file */
  for (i=1; feof(stream2) == 0 && i < 5000; i++) {
    if (iwaussian) {
          fgets(line,200, stream2); 
    }
    else {
      fgets(line,150, stream2); 
    }
    if(strncmp(line, "#", 1) != 0) {
      if (skipcol == 0)   {
	sscanf(line,"%f %f %f",&x2[i],&y2[i],&amag2[i]);
      }
      else {
	if (iwaussian) {	 
	  sscanf(line,"%s %f %f %s %s %s %s %s %s %s %s %f",dummyc,&x2[i],&y2[i],dummyc,dummyc,dummyc,dummyc,dummyc,dummyc,dummyc,dummyc,&amag2[i]);
	  if (amag2[i] > 0) {
	    amag2[i] = -2.5*log(amag2[i])/log(10);
	  }
	  else {
	    amag2[i]=99.0;	    
	  }
	}
	else {
	  sscanf(line,"%s %f %f %f",dummyc,&x2[i],&y2[i],&amag2[i]);
	}
      }
      if (js > 0) { /*we are jumpstarting*/
	for (j=0;j<js;j++) {
	  if (strcasecmp(jumpstartid2[j],dummyc) == 0 ) {
	    jumpstartref2[j]=i;
	    if (diagnose) {
	      printf("jumpstart star %s given reference %d for %s\n",jumpstartid2[j],i,file2);
	    }
	  }
	}
      }


      if (amag2[i] ==0) amag2[i]=99.99;
      if (x2[i]==0 && y2[i]==0) i--;
    }
    else i--;
  }
  nstar2 = i-1;
  if (nstar2 < 3)  {
    printf("\n%s only has %d Star/s...skipping it\n",file2,nstar2);
    exit(-1);
  }


  indexx(nstar2,amag2,index2); /* sort by magnitudes  - give sort into index1 */

  fail=0;
  /*test that we got our jump stars in list*/
  for (j=0;j<js;j++) {
    if (jumpstartref2[j] < 0) {
      printf("Didn't find jumpstart star %s in %s\n",jumpstartid2[j],file2);
      fail=1;
    }
    for (k=1;k<=nstar2;k++) {
      if (index2[k] == jumpstartref2[j]) {
	jumpstartref2[j]=k;
	k=nstar2+1;
      }
    }
  }
  if (fail) {
    exit(-1);
  }


   
  fclose(stream1);
  fclose(stream2);
    

  if (diagnose) {
    printf("read in %d stars from %s\n",nstar2,file2);
  }  



  /* now we have four cases mainly because I am simple minded in my
     programming*/

  /*These cases are 0,1,2, and 3+ jumpstart stars given*/
  /* Case 0 - go through and calculate all triangles, starting
     with the brightest first

     Case 1 - Make all triangles that include star 1 - this means nstar*(nstar
     -1)/2 triangles
  
     Case 2- Make all triangles that include stars 1 and 2 - only N triangles
       
     Case 3+ - just use these stars as though we had successfully got our 
     solution via cases 0-2 */



  /* Calculate Triangle Ratios */


  solution=0;
  fail=0;

  if (js >3) jstmp=3;
  else jstmp=js;
  
  for (nstart2=1;nstart2 < nstar2 && (! solution) ;nstart2+=nobj) {          
    for (nstart1=1;nstart1 < nstar1 && (! solution);nstart1 +=nobj) {
      
      switch(jstmp) {
      case 0:
	
	/* 1st file */ 
	for (kk=1,i=nstart1; i<nstart1+nobj-2-1 && i <= nstar1 && (! solution); i++) {
	  for (j=i+1; j < nstart1+nobj-1-1 && j <= nstar1; j++) {
	    r12=R12 = sqrt(((x1[index1[i]]-x1[index1[j]])*(x1[index1[i]]-x1[index1[j]])) +
		       ((y1[index1[i]]-y1[index1[j]])*(y1[index1[i]]-y1[index1[j]])));
	    if (r12 > sepmin && i != j) {
	      for(k=j+1;k<nstart1+nobj-1 && k <= nstar1;k++) {
		R13=r13 = sqrt(((x1[index1[i]]-x1[index1[k]])*(x1[index1[i]]-x1[index1[k]])) +
			   ((y1[index1[i]]-y1[index1[k]])*(y1[index1[i]]-y1[index1[k]])));
		if (r13 > sepmin && i != k ) {
		  r23 = sqrt(((x1[index1[j]]-x1[index1[k]])*(x1[index1[j]]-x1[index1[k]])) 
			     + ((y1[index1[j]]-y1[index1[k]])*(y1[index1[j]]-y1[index1[k]])));
		  if (r23 > sepmin && j != k ) {
		    r12=R12;
		    r13=R13;
		    if (r12 > r13) { /*sort from to small large*/
		      temp=r12;
		      r12=r13;
		      r13=temp;
		    }
		    
		    if (r12 > r23) {		  
		      temp=r12 ;
		      r12=r23;
		      r23=temp;
		    }
		    
		    if (r13 > r23) {		  
		      temp=r13;
		      r13=r23;
		      r23=temp;
		    }
		    
		    
		    r1[1][kk] = r12/r13;
		    r1[2][kk] = r12/r23;
		    r1[3][kk] = r13/r23;
		    
		    rsum1[kk] = r1[1][kk]+r1[2][kk]+r1[3][kk];
		    
		    if (rsum1[kk] > rsummin) {
		      rabs1[1][kk] = r12;
		      rabs1[2][kk] = r13;
		      rabs1[3][kk] = r23;
		      ijk1[1][kk] = i;
		      ijk1[2][kk] = j;
		      ijk1[3][kk] = k;
		      kk ++;
		    }
		  }
		}
	      }
	    }
	  }
	}
	nsum1 = kk - 1;
	if (nsum1 <2) nsum1 =2; 
	indexx(nsum1,rsum1,indexr1);
	if (diagnose > 1) {
	  printf("%d triangles from starid %d %d in file %s\n",nsum1,nstart1,nstart1+nobj-1,file1);
	}  
	

      
	/* 2nd file */
      
	for (kk=1,i=nstart2; i<nstart2+nobj-2-1 && i <= nstar2 && (! solution); i++) {
	  for (j=i+1; j < nstart2+nobj-1-1 && j <=nstar2 && (! solution); j++) {
	    R12=r12 = sqrt(((x2[index2[i]]-x2[index2[j]])*(x2[index2[i]]-x2[index2[j]])) +
		       ((y2[index2[i]]-y2[index2[j]])*(y2[index2[i]]-y2[index2[j]]))); 
	    if (r12 > sepmin && i != j) {  
	      for(k=j+1;k<nstart2+nobj-1 && k <= nstar2 && (! solution);k++) {
		R13=r13 = sqrt(((x2[index2[i]]-x2[index2[k]])*(x2[index2[i]]-x2[index2[k]])) +
			   ((y2[index2[i]]-y2[index2[k]])*(y2[index2[i]]-y2[index2[k]])));
		if (r13 > sepmin && i !=k) {
		  r23 = sqrt(((x2[index2[j]]-x2[index2[k]])*(x2[index2[j]]-x2[index2[k]])) +
			     ((y2[index2[j]]-y2[index2[k]])*(y2[index2[j]]-y2[index2[k]])));
		  if (r23 > sepmin && j !=k) {
		    r12=R12;
		    r13=R13;
		    if (r12 > r13) { /*sort from to small large*/
		      temp=r12;
		      r12=r13;
		      r13=temp;
		    }
		    
		    if (r12 > r23) {		  
		      temp=r12 ;
		      r12=r23;
		      r23=temp;
		    }
		    
		    if (r13 > r23) {		  
		      temp=r13;
		      r13=r23;
		      r23=temp;
		    }
		    
		  
		    r2[1][kk] = r12/r13;
		    r2[2][kk] = r12/r23;
		    r2[3][kk] = r13/r23;
		    
		    rsum2[kk] = r2[1][kk]+r2[2][kk]+r2[3][kk];
		    if (rsum2[kk] > rsummin) {
		      rabs2[1][kk] = r12;
		      rabs2[2][kk] = r13;
		      rabs2[3][kk] = r23;
		      ijk2[1][kk] = i;
		      ijk2[2][kk] = j;
		      ijk2[3][kk] = k;
		      kk ++;
		    }
		  }
		}
	      }
	    }
	  }
	}
	nsum2 = kk -1;
	
	if (nsum2 < 2) nsum2 =2; 
	/*	indexx(nsum2,rsum2,indexr2);*/
	
	if (diagnose > 1) {
	  printf("%d triangles from starid %d %d in file %s\n",nsum2,nstart2,nstart2+nobj-1,file2);
	}  
	break;
      
      case 1:
	/* 1st file */ 
	for (kk=1,i=nstart1; i<nstart1+nobj-2-1 && i <= nstar1 && (! solution); i++) {
	  for (j=i+1; j < nstart1+nobj-1-1 && j <= nstar1; j++) {
	    r12 = sqrt(((x1[index1[i]]-x1[index1[j]])*(x1[index1[i]]-x1[index1[j]])) +
		       ((y1[index1[i]]-y1[index1[j]])*(y1[index1[i]]-y1[index1[j]])));
	    if (r12 > sepmin) {
	      k=jumpstartref1[0]; /*this is always one of our points ERROR HERE*/
	      r13 = sqrt(((x1[index1[i]]-x1[index1[k]])*(x1[index1[i]]-x1[index1[k]])) +
			 ((y1[index1[i]]-y1[index1[k]])*(y1[index1[i]]-y1[index1[k]])));
	      if (r13 > sepmin && i != k ) {
		r23 = sqrt(((x1[index1[j]]-x1[index1[k]])*(x1[index1[j]]-x1[index1[k]])) 
			   + ((y1[index1[j]]-y1[index1[k]])*(y1[index1[j]]-y1[index1[k]])));
		if (r23 > sepmin && j != k ) {
		  if (r12 > r13) { /*sort from to small large*/
		    temp=r12;
		    r12=r13;
		    r13=temp;
		  }
		  
		  if (r12 > r23) {		  
		    temp=r12 ;
		    r12=r23;
		    r23=temp;
		  }
		  
		  if (r13 > r23) {		  
		    temp=r13;
		    r13=r23;
		    r23=temp;
		  }
		  
		  
		  r1[1][kk] = r12/r13;
		  r1[2][kk] = r12/r23;
		  r1[3][kk] = r13/r23;
		  
		  rsum1[kk] = r1[1][kk]+r1[2][kk]+r1[3][kk];
		  
		  if (rsum1[kk] > rsummin) {
		    rabs1[1][kk] = r12;
		    rabs1[2][kk] = r13;
		    rabs1[3][kk] = r23;
		    ijk1[1][kk] = i;
		    ijk1[2][kk] = j;
		    ijk1[3][kk] = k;
		    kk ++;
		  }
		}
	      }
	    }
	  }
	}
	nsum1 = kk - 1;
	if (nsum1 <2) nsum1 =2; 
	indexx(nsum1,rsum1,indexr1);
	if (diagnose > 1) {
	  printf("%d triangles from starid %d %d in file %s\n",nsum1,nstart1,nstart1+nobj-1,file1);
	}  
	
      
	/* 2nd file */
	for (kk=1,i=nstart2; i<nstart2+nobj-2-1 && i <= nstar2 && (! solution); i++) {      
	  for (j=i+1; j < nstart2+nobj-1-1 && j <=nstar2 && (! solution); j++) {
	    r12 = sqrt(((x2[index2[i]]-x2[index2[j]])*(x2[index2[i]]-x2[index2[j]])) +
		       ((y2[index2[i]]-y2[index2[j]])*(y2[index2[i]]-y2[index2[j]])));
	    if (r12 > sepmin) {
	      k=jumpstartref2[0]; /*this is always one of our points*/
	      r13 = sqrt(((x2[index2[i]]-x2[index2[k]])*(x2[index2[i]]-x2[index2[k]])) +
			 ((y2[index2[i]]-y2[index2[k]])*(y2[index2[i]]-y2[index2[k]])));
	      if (r13 > sepmin && i !=k) {
		r23 = sqrt(((x2[index2[j]]-x2[index2[k]])*(x2[index2[j]]-x2[index2[k]])) +
			   ((y2[index2[j]]-y2[index2[k]])*(y2[index2[j]]-y2[index2[k]])));
		if (r23 > sepmin && j !=k) {
		  if (r12 > r13) { /*sort from to small large*/
		    temp=r12;
		    r12=r13;
		    r13=temp;
		  }
		  
		  if (r12 > r23) {		  
		    temp=r12 ;
		    r12=r23;
		    r23=temp;
		  }
		  
		  if (r13 > r23) {		  
		    temp=r13;
		    r13=r23;
		    r23=temp;
		  }
		  
		  
		  r2[1][kk] = r12/r13;
		  r2[2][kk] = r12/r23;
		  r2[3][kk] = r13/r23;
		  
		  rsum2[kk] = r2[1][kk]+r2[2][kk]+r2[3][kk];
		  if (rsum2[kk] > rsummin) {
		    rabs2[1][kk] = r12;
		    rabs2[2][kk] = r13;
		    rabs2[3][kk] = r23;
		    ijk2[1][kk] = i;
		    ijk2[2][kk] = j;
		    ijk2[3][kk] = k;
		    kk ++;
		  }
		}
	      }
	    }
	  }
	}
      
      
	nsum2 = kk -1;
      
	if (nsum2 < 2) nsum2 =2; 
	/*	indexx(nsum2,rsum2,indexr2);*/
	
	if (diagnose > 1) {
	  printf("%d triangles from starid %d %d in file %s\n",nsum2,nstart2,nstart2+nobj-1,file2);
	}  
	break;
    
    
      case 2:
	/* 1st file */ 
	for (kk=1,i=nstart1; i<nstart1+nobj-2-1 && i <= nstar1 && (! solution); i++) {
	  j=jumpstartref1[1];
	  r12 = sqrt(((x1[index1[i]]-x1[index1[j]])*(x1[index1[i]]-x1[index1[j]])) +
		     ((y1[index1[i]]-y1[index1[j]])*(y1[index1[i]]-y1[index1[j]])));
	  if (r12 > sepmin) {
	    k=jumpstartref1[0]; /*this is always one of our points*/
	    r13 = sqrt(((x1[index1[i]]-x1[index1[k]])*(x1[index1[i]]-x1[index1[k]])) +
		       ((y1[index1[i]]-y1[index1[k]])*(y1[index1[i]]-y1[index1[k]])));
	    if (r13 > sepmin && i != k ) {
	      r23 = sqrt(((x1[index1[j]]-x1[index1[k]])*(x1[index1[j]]-x1[index1[k]])) 
			 + ((y1[index1[j]]-y1[index1[k]])*(y1[index1[j]]-y1[index1[k]])));
	      if (r23 > sepmin && j != k ) {
		if (r12 > r13) { /*sort from to small large*/
		  temp=r12;
		  r12=r13;
		  r13=temp;
		}
		
		if (r12 > r23) {		  
		  temp=r12 ;
		  r12=r23;
		  r23=temp;
		}
		
		if (r13 > r23) {		  
		  temp=r13;
		  r13=r23;
		  r23=temp;
		}
		
		
		r1[1][kk] = r12/r13;
		r1[2][kk] = r12/r23;
		r1[3][kk] = r13/r23;
		
		rsum1[kk] = r1[1][kk]+r1[2][kk]+r1[3][kk];
		
		if (rsum1[kk] > rsummin) {
		  rabs1[1][kk] = r12;
		  rabs1[2][kk] = r13;
		  rabs1[3][kk] = r23;
		  ijk1[1][kk] = i;
		  ijk1[2][kk] = j;
		  ijk1[3][kk] = k;
		  kk ++;
		}
	      }
	    }
	  }
	}
	nsum1 = kk - 1;
	if (nsum1 <2) nsum1 =2; 
	indexx(nsum1,rsum1,indexr1);
	if (diagnose > 1) {
	  printf("%d triangles from starid %d %d in file %s\n",nsum1,nstart1,nstart1+nobj-1,file1);
	}  
	
      	/* 2nd file */
	for (kk=1,i=nstart2; i<nstart2+nobj-2-1 && i <= nstar2 && (! solution); i++) {      
	  j=jumpstartref2[1];
	  r12 = sqrt(((x2[index2[i]]-x2[j])*(x2[index2[i]]-x2[j])) +
		     ((y2[index2[i]]-y2[j])*(y2[index2[i]]-y2[j])));
	  if (r12 > sepmin) {
	    k=jumpstartref2[0]; /*this is always one of our points*/
	    r13 = sqrt(((x2[index2[i]]-x2[index2[k]])*(x2[index2[i]]-x2[index2[k]])) +
		       ((y2[index2[i]]-y2[index2[k]])*(y2[index2[i]]-y2[index2[k]])));
	    if (r13 > sepmin && i !=k) {
	      r23 = sqrt(((x2[index2[j]]-x2[index2[k]])*(x2[index2[j]]-x2[index2[k]])) +
			 ((y2[index2[j]]-y2[index2[k]])*(y2[index2[j]]-y2[index2[k]])));
	      if (r23 > sepmin && j !=k) {
		if (r12 > r13) { /*sort from to small large*/
		  temp=r12;
		  r12=r13;
		  r13=temp;
		}
		
		if (r12 > r23) {		  
		  temp=r12 ;
		  r12=r23;
		  r23=temp;
		}
		
		if (r13 > r23) {		  
		  temp=r13;
		  r13=r23;
		  r23=temp;
		}
		
		  
		r2[1][kk] = r12/r13;
		r2[2][kk] = r12/r23;
		r2[3][kk] = r13/r23;
		
		rsum2[kk] = r2[1][kk]+r2[2][kk]+r2[3][kk];
		if (rsum2[kk] > rsummin) {
		  rabs2[1][kk] = r12;
		  rabs2[2][kk] = r13;
		  rabs2[3][kk] = r23;
		  ijk2[1][kk] = i;
		  ijk2[2][kk] = j;
		  ijk2[3][kk] = k;
		  kk ++;
		}
	      }
	    }
	  }
	}
	nsum2 = kk -1;
	
	if (nsum2 < 2) nsum2 =2; 
	/*indexx(nsum2,rsum2,indexr2); */ /*why do this? we really do not want to search by triange shape do we?*/
	
	
	if (diagnose > 1) {
	  printf("%d triangles from starid %d %d in file %s\n",nsum2,nstart2,nstart2+nobj-1,file2);
	}  
	break;
      
      case 3:
	solution=2;/*two means successful, but no triangles*/
	
	for (j=1;j<=js;j++) {
	  xf1[j]= x1[jumpstartref1[j-1]]; /*we are mixing 0,1 indices here*/
	  yf1[j]= y1[jumpstartref1[j-1]]; /* yuk!*/
	  
	  xf2[j]= x2[jumpstartref2[j-1]];
	  yf2[j]= y2[jumpstartref2[j-1]];
	}	  
	nf=js;
	break;
      }
      
      /* I know have two lists of triangles for the two images*/
      

      /*now start looking for similar triangles*/
      /*go over each triangle in 2nd image*/
      /*assuming we do not already have a solution*/
      
      for (i=1;i <= nsum2 && (! solution); i++) {
	
	xhunt = rsum2[i] - tolerance; /*Search list in template file for closest */
	hunt(rsum1,nsum1,xhunt,&jlow); /* Match to the sums in active file */
	
	if (jlow > 1) jlow--;
	xhunt = rsum2[i] + tolerance;
	hunt(rsum1,nsum1,xhunt,&jhigh);
	if (jhigh < nsum1) jhigh++;
	
	/*so relevant triangles are between jlo and jhi*/
	
	tol = tolerance;
	jbest = 0;
	
	for (jbest=jlow; jbest <= jhigh+1 && jbest > 0 && (! solution) ; jbest++) {  /* find best in the range give by Hunt */
	  
	  
	  xf1[1]= a= x1[index1[ijk1[1][indexr1[jbest]]]];
	  yf1[1]= b= y1[index1[ijk1[1][indexr1[jbest]]]];
	  
	  xf1[2]= c= x1[index1[ijk1[2][indexr1[jbest]]]];
	  yf1[2]= d= y1[index1[ijk1[2][indexr1[jbest]]]];
	  
	  xf1[3]= e= x1[index1[ijk1[3][indexr1[jbest]]]];
	  yf1[3]= f= y1[index1[ijk1[3][indexr1[jbest]]]];
	  
	  
	  side1=sqrt((xf1[2]-xf1[3])*(xf1[2]-xf1[3])+(yf1[2]-yf1[3])*(yf1[2]-yf1[3]));
	  side2=sqrt((xf1[3]-xf1[1])*(xf1[3]-xf1[1])+(yf1[3]-yf1[1])*(yf1[3]-yf1[1]));
	  side3=sqrt((xf1[2]-xf1[1])*(xf1[2]-xf1[1])+(yf1[2]-yf1[1])*(yf1[2]-yf1[1]));
	  
	  
	  /*dumb sort to figure out to order the triangle points*/
	  if (side1 > side2 && side1 > side3 && side2 > side3) { /*3,2,1*/
	    xf1[1]=e; yf1[1]=f; xf1[2]=c;yf1[2]=d;xf1[3]=a;yf1[3]=b;
	    longside1=side1;
	  }
	  if (side1 > side2 && side1 > side3 && side3 > side2) { /*2,3,1*/
	    xf1[1]=c; yf1[1]=d; xf1[2]=e;yf1[2]=f;xf1[3]=a;yf1[3]=b;
	    longside1=side1;
	  }
	  if (side2 > side1 && side2 > side3 && side3 > side1) { /*1,3,2*/
	    xf1[1]=a; yf1[1]=b; xf1[2]=e;yf1[2]=f;xf1[3]=c;yf1[3]=d;
	    longside1=side2;
	  }
	  if (side2 > side1 && side2 > side3 && side1 > side3) { /*3,1,2*/
	    xf1[1]=e; yf1[1]=f; xf1[2]=a;yf1[2]=b;xf1[3]=c;yf1[3]=d;
	    longside1=side2;
	  }
	  if (side3 > side1 && side3 > side2 && side1 > side2) { /*2,1,3*/
	    xf1[1]=c; yf1[1]=d; xf1[2]=a;yf1[2]=b;xf1[3]=e;yf1[3]=f;
	    longside1=side3;
	  }
	  if (side3 > side1 && side3 > side2 && side2 > side1) { /*1,2,3*/
	    xf1[1]=a; yf1[1]=b; xf1[2]=c;yf1[2]=d;xf1[3]=e;yf1[3]=f;
	    longside1=side3;
	  }
	  
	  
	  /*	  xf2[1]=	  a = x2[index2[ijk2[1][indexr2[i]]]];
	  yf2[1]=	  b = y2[index2[ijk2[1][indexr2[i]]]];
	  
	  xf2[2]=	  c = x2[index2[ijk2[2][indexr2[i]]]];
	  yf2[2]=	  d = y2[index2[ijk2[2][indexr2[i]]]];
	  
	  xf2[3]=	  e = x2[index2[ijk2[3][indexr2[i]]]];
	  yf2[3]=	  f = y2[index2[ijk2[3][indexr2[i]]]];  */

	  xf2[1]=	  a = x2[index2[ijk2[1][i]]];
	  yf2[1]=	  b = y2[index2[ijk2[1][i]]];
	  					  
	  xf2[2]=	  c = x2[index2[ijk2[2][i]]];
	  yf2[2]=	  d = y2[index2[ijk2[2][i]]];
	  					  
	  xf2[3]=	  e = x2[index2[ijk2[3][i]]];
	  yf2[3]=	  f = y2[index2[ijk2[3][i]]];


	 
	  side1=sqrt((xf2[2]-xf2[3])*(xf2[2]-xf2[3])+(yf2[2]-yf2[3])*(yf2[2]-yf2[3]));
	  side2=sqrt((xf2[3]-xf2[1])*(xf2[3]-xf2[1])+(yf2[3]-yf2[1])*(yf2[3]-yf2[1]));
	  side3=sqrt((xf2[2]-xf2[1])*(xf2[2]-xf2[1])+(yf2[2]-yf2[1])*(yf2[2]-yf2[1]));
	  
	  
	  /*dumb sort to figure out to order the triangle points*/
	  if (side1 > side2 && side1 > side3 && side2 > side3) { /*3,2,1*/
	    xf2[1]=e; yf2[1]=f; xf2[2]=c;yf2[2]=d;xf2[3]=a;yf2[3]=b;
	    longside2=side1;
	  }
	  if (side1 > side2 && side1 > side3 && side3 > side2) { /*2,3,1*/
	    xf2[1]=c; yf2[1]=d; xf2[2]=e;yf2[2]=f;xf2[3]=a;yf2[3]=b;
	    longside2=side1;
	  }
	  if (side2 > side1 && side2 > side3 && side3 > side1) { /*1,3,2*/
	    xf2[1]=a; yf2[1]=b; xf2[2]=e;yf2[2]=f;xf2[3]=c;yf2[3]=d;
	    longside2=side2;
	  }
	  if (side2 > side1 && side2 > side3 && side1 > side3) { /*3,1,2*/
	    xf2[1]=e; yf2[1]=f; xf2[2]=a;yf2[2]=b;xf2[3]=c;yf2[3]=d;
	    longside2=side2;
	  }
	  if (side3 > side1 && side3 > side2 && side1 > side2) { /*2,1,3*/
	    xf2[1]=c; yf2[1]=d; xf2[2]=a;yf2[2]=b;xf2[3]=e;yf2[3]=f;
	    longside2=side3;
	  }
	  if (side3 > side1 && side3 > side2 && side2 > side1) { /*1,2,3*/
	    xf2[1]=a; yf2[1]=b; xf2[2]=c;yf2[2]=d;xf2[3]=e;yf2[3]=f;
	    longside2=side3;
	  }
	
	  /*matching triangles are now in vectors xf1,xf2; yf1,yf2*/
	  
	  if (scaleratio == 0 || (fabs(longside1/longside2/scaleratio-1.) < srtolerance)) {
	    /*scale is OK, now do the transform*/
	    transform(3,xf1,yf1,xf2,yf2,&ax,&bx,&cx,&ay,&by,&cy);
	    
	    /*is it still right scale*/
	    if ((bx >-1e10 && bx < 1e10) && (scaleratio == 0 || fabs((sqrt((bx*bx+cx*cx+by*by+cy*cy)/2)/scaleratio)-1) < srtolerance)) { /*occaisional problem in transform causes corruption of xf2*/
	    /* is transform square*/
	      
	      a1=atan2(bx,cx)*180/3.1415926535-90;
	      a2=atan2(by,cy)*180/3.1415926535-90;
	      a3=a1-a2-90;
	      while (a3 < -180) {a3+=360;} /*normalise our angle*/
	      while (a3 > 180) {a3-=360;}
	      

	      if (fabs((bx*bx+cx*cx)/(by*by+cy*cy)-1) < srtolerance) { /*is it square*/
		
		/* does transform have Ok angle e.g. is it E-W etc or anything*/
		/* within 90 degrees of this*/
		a3=a1-rotation;
		
		while (a3 < -45) {a3+=90;} /*normalise our angle*/
		while (a3 > 45) {a3-=90;} /*normalise our angle*/
		
		if (fabs(a3) < ROTERR) {
		  
		  if (diagnose > 1) {
		    printf("OKtriangle:%5d %10.4f %10.4f %10.2f %10.4f %10.4f %10.3f\n",r,bx,cx,ax,by,cy,ay);
		    transform(3,xf1,yf1,xf2,yf2,&ax,&bx,&cx,&ay,&by,&cy);
		  }
		  
		  /* load up our 6 param fit and starx,ys*/
		  alist[r]=ax;
		  blist[r]=bx;
		  clist[r]=cx;
		  dlist[r]=ay;
		  elist[r]=by;
		  flist[r]=cy;
		  xfm1[1][r]=xf1[1];
		  xfm1[2][r]=xf1[2];
		  xfm1[3][r]=xf1[3];
		  yfm1[1][r]=yf1[1];
		  yfm1[2][r]=yf1[2];
		  yfm1[3][r]=yf1[3];
		  
		  xfm2[1][r]=xf2[1];
		  xfm2[2][r]=xf2[2];
		  xfm2[3][r]=xf2[3];
		  yfm2[1][r]=yf2[1];
		  yfm2[2][r]=yf2[2];
		  yfm2[3][r]=yf2[3];
		  
		  matches[r]=0;
		  
		  r++;
		  
		  
		  /*check for matching ones*/

		  scaletmp=sqrt(bx*bx+cx*cx);
		  for (j=0;j<r-1;j++) {
	 
		    atest=fabs(alist[j]-ax)/(fabs(ax)+10*(scaletmp/0.5));
		    dtest=fabs(dlist[j]-ay)/(fabs(ay)+10*(scaletmp/0.5));
		    btest=fabs(blist[j]-bx)/scaletmp;
		    ctest=fabs(clist[j]-cx)/scaletmp;
		    etest=fabs(elist[j]-by)/scaletmp;
		    ftest=fabs(flist[j]-cy)/scaletmp;
		    sumtest=atest+dtest+btest+ctest+etest+ftest;

		    if (atest < amiss && dtest < amiss && btest < cbmiss && ctest < cbmiss && etest < cbmiss &&
			ftest < cbmiss && sumtest < totalmiss) { 
		      matches[j]++;
		      if (matches[j] > bestmatch ) {
			bestmatch=matches[j];
			matchid=j;
			if (diagnose > 0) {
			  printf("Now have %d matches\n",matches[j]);
			  printf("Caught a:%7.3f d:%7.3f b:%7.3f c:%7.3f e:%7.3f f:%7.3f\n",atest,dtest,btest,ctest,dtest,etest,ftest);
			  printf("This a:%7.3f d:%7.3f b:%7.3f c:%7.3f e:%7.3f f:%7.3f\n",alist[j],dlist[j],blist[j],clist[j],elist[j],flist[j]);
			  printf("With a:%7.3f d:%7.3f b:%7.3f c:%7.3f e:%7.3f f:%7.3f\n",ax,ay,bx,cx,by,cy);
			}
			if (matches[j] > neededmatches) { /*we are done here*/
			  solution=1;
			}
		      }
		    }
		    else {
		      if (diagnose > 2 ) {
			printf("missed a:%7.2f d:%7.2f b:%7.2f c:%7.2f e:%7.2f f:%7.2f\n",atest,dtest,btest,ctest,dtest,etest,ftest);
		      }
		    }
		  }
		}
	      }
	    }
	  } 
	} 
      }
   }
  }
  if (! solution) {
    printf("\nMATCH STAR FAILED\n"); 
    exit(0);      
  }
  
  
  /* I have enough matching triangles. */
  /*get list of stars*/
  if (solution ==1) { /* solution=2 skip this*/
    jj=1;
    
    scaletmp=sqrt(blist[matchid]*blist[matchid]+clist[matchid]*clist[matchid]);
    for (j=0;j<=r;j++) {
      if (fabs(alist[j]-alist[matchid])/((fabs(alist[matchid])+10*(scaletmp/0.5))) < amiss &&  
	  fabs(dlist[j]-dlist[matchid])/((fabs(dlist[matchid])+10*(scaletmp/0.5))) < amiss &&
	  fabs(blist[j]-blist[matchid])/scaletmp < cbmiss &&
	  fabs(clist[j]-clist[matchid])/scaletmp < cbmiss &&
	  fabs(elist[j]-elist[matchid])/scaletmp < cbmiss &&
	  fabs(flist[j]-flist[matchid])/scaletmp < cbmiss) { /* this is one of the matching*/

	
	xf1[jj]=xfm1[1][j];
	xf1[jj+1]=xfm1[2][j];
	xf1[jj+2]=xfm1[3][j];
	
	yf1[jj]=yfm1[1][j];
	yf1[jj+1]=yfm1[2][j];
	yf1[jj+2]=yfm1[3][j];
	
	xf2[jj]=xfm2[1][j];
	xf2[jj+1]=xfm2[2][j];
	xf2[jj+2]=xfm2[3][j];
	
	yf2[jj]=yfm2[1][j];
	yf2[jj+1]=yfm2[2][j];
	yf2[jj+2]=yfm2[3][j];
	
	jj+=3;
      }
    }
    nf=jj-1;
  }
  /* Do Multvariate regression to find tranformation coefficients */
  
  /*use regression to find all stars that are probably acceptable*/
  transform(nf,xf1,yf1,xf2,yf2,&ax,&bx,&cx,&ay,&by,&cy);
  scaleratio=(float) sqrt(bx*bx+cx*cx);
  
  sumdx = 0.0;
  sumdy = 0.0;
  sumdx2 = 0.0;
  sumdy2 = 0.0;
    
  for (i=1; i<=nf;i++) {       /* Calc sigma */
    x = ax +bx*xf2[i] + cx*yf2[i];
    y = ay +by*xf2[i] + cy*yf2[i];
    
    diffx[i] = x - xf1[i];
    diffy[i] = y - yf1[i];
    
    sumdx += diffx[i];
    sumdx2 += diffx[i]*diffx[i];
    sumdy += diffy[i];
    sumdy2 += diffy[i]*diffy[i];
    if (diagnose >1) {
      printf("\tStar in fit %5d%9.2f%9.2f%9.2f%9.2f%9.3f%9.3f\n",i,xf2[i],yf2[i],x,y,diffx[i],diffy[i]);
    }
  }  
  anf = nf;
  sdx = sqrt((1.0/(anf-1.0))*(sumdx2-(sumdx*sumdx/anf)))+scaleratio*0.01;
  sdy = sqrt((1.0/(anf-1.0))*(sumdy2-(sumdy*sumdy/anf)))+scaleratio*0.01;
  if (sdx < scaleratio*0.1) { /*don't let get smaller than 0.1 of a pixel*/
    sdx=scaleratio*0.1;
  }
  if (sdy < scaleratio*0.1) { /*don't let get smaller than 0.1 of a pixel*/
    sdy=scaleratio*0.1;
  }
  
  /*transform completed first iteration from triangles*/
  /*now output all stars that seem reasonable doing a iterative sig clip*/
  nstarsold=9999999;
  nstars=nf;
  iterations=0;

  while(nstarsold != nstars && iterations < maxiterations) { /*have we converged yet?*/
    if (diagnose>0) {
      printf("%.4f %.4f %.3f\n%.4f %.4f %.3f\nsigx=%.2f sigy=%.2f tri=%d stars=%d\n",bx,cx,ax,by,cy,ay,sdx,sdy,bestmatch,nstars); 
    }
    /* match stars in the two files using our transform*/
    nstarsold=nstars;
    nstars=1;

    OKx=sdx*sigmatch+floormatch;
    OKy=sdy*sigmatch+floormatch;
    if (diagnose > 1) {
      printf("OKx:%7.3f sigmax: %7.3f OKy: %7.3f sigmay: %7.3f\n",OKx,sdx,OKy,sdy);
    }
    for (i=0;i<=nstar2;i++) {
      x = x2[i];
      y = y2[i];
     
      XC = ax + bx*x +cx*y;  
      YC = ay + by*x + cy*y; /* determine this star's position in other file */

      for (j=0; j < nstar1; j++) { /* does it match a star in file 1?, if yes, keep it*/
	if ( fabs(XC-x1[j]) < OKx && fabs(YC-y1[j]) < OKy) {
	  if (diagnose>1) {
	    printf("\tStar in fit %5d%9.2f%9.2f%9.2f%9.2f%9.3f%9.3f\n",nstars,x,y,x1[j],y1[j],XC-x1[j],YC-y1[j]);
	  }
	  xf1[nstars]=x1[j];
	  xf2[nstars]=x;
	  yf1[nstars]=y1[j];
	  yf2[nstars]=y;
	  nstars++;
	  j = nstar1+1;
	}
      } 
    }
    nstars-=1;
    
    /*recalculate with whole enchilada of stars*/
    if( diagnose >2) {
      for (i=1;i<=nstars;i++) {
	printf("\tstar in transform %5d%9.2f%9.2f%9.2f%9.2f\n",i,xf2[i],yf2[i],xf1[i],yf1[i]);
      }
    }
    transform(nstars,xf1,yf1,xf2,yf2,&ax,&bx,&cx,&ay,&by,&cy);
    
    /*recalcuate sigma*/
    
    nf=nstars;
    sumdx = 0.0;
    sumdy = 0.0;
    sumdx2 = 0.0;
    sumdy2 = 0.0;
    
    for (i=1; i<=nf;i++) {       /* Calc sigma */
      x = ax +bx*xf2[i] + cx*yf2[i];
      y = ay +by*xf2[i] + cy*yf2[i];
      
      diffx[i] = x - xf1[i];
      diffy[i] = y - yf1[i];
      
      sumdx += diffx[i];
      sumdx2 += diffx[i]*diffx[i];
      sumdy += diffy[i];
      sumdy2 += diffy[i]*diffy[i];
      
    }  
    anf = nf;
    sdx = sqrt((1.0/(anf-1.0))*(sumdx2-(sumdx*sumdx/anf)));
    sdy = sqrt((1.0/(anf-1.0))*(sumdy2-(sumdy*sumdy/anf)));
    iterations++;
  }
  
  printf("%.4f %.4f %.3f %.4f %.4f %.3f sigx=%.2f sigy=%.2f tri=%d stars=%d\n",bx,cx,ax,by,cy,ay,sdx,sdy,bestmatch,nstars); 
  stream3 = fopen(outfile,"w");          
  for (i=0;i<=nstar2;i++) {
    x = x2[i];
    y = y2[i];
    
    XC = ax + bx*x +cx*y;  
    YC = ay + by*x + cy*y; /* determine this star's position in other file */
     
    OKx=sdx*sigmatch+floormatch;
    OKy=sdy*sigmatch+floormatch;
     
    for (j=0; j < nstar1; j++) { /* does it match a star in file 1?, if yes, write it out*/
      if ( fabs(XC-x1[j]) < OKx && fabs(YC-y1[j]) < OKy) {
	if (outids == 1) { 
	  fprintf(stream3,"\n%9.2f%9.2f%9.2f%9.2f %8i %8i",x,y,x1[j],y1[j],i,j);
	}
	else {
	  if (xieta==1) {
	    fprintf(stream3,"%9.2f\t%9.2f\t%11.7f\t%11.7f\t%9.3f\t%9.3f\n",x,y,xi[j],eta[j],XC-x1[j],YC-y1[j]);}
	  else {
	    fprintf(stream3,"\n%9.2f%9.2f%9.2f%9.2f%9.3f%9.3f",x,y,x1[j],y1[j],XC-x1[j],YC-y1[j]);
	  }
	}  
      }
    }
  }
  fclose(stream3);
}


void transform(nstar,x2,y2,x1,y1,ax,bx,cx,ay,by,cy)
     int nstar;
     float x2[],y2[],x1[],y1[],*ax,*bx,*cx,*ay,*by,*cy;
{

  int i;
  double sum=0.,sumx=0.,sumy=0.,sumxy=0.,sumu=0.,sumv=0.,sumx2=0.,sumy2=0.;
  double sumux=0.,sumuy=0.,sumvx=0.,sumvy=0.;
  double delta,sumdx,sumdy,sumabsdx,sumabsdy;
  double x1p,y1p,dx,dy;
  for (i=1; i<=nstar;i++) {
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
  delta = sum*(sumx2*sumy2 - sumxy*sumxy) + sumx*(sumxy*sumy - sumx*sumy2) + sumy*(sumx*sumxy-sumx2*sumy);
  *ax = (sumu*(sumx2*sumy2 - sumxy*sumxy) + sumx*(sumxy*sumuy - sumux*sumy2) + sumy*(sumux*sumxy - sumx2*sumuy))/delta;
  *bx = (sum*(sumux*sumy2 - sumxy*sumuy) + sumu*(sumxy*sumy - sumx*sumy2) + sumy*(sumuy*sumx - sumy*sumux))/delta;
  *cx = (sum*(sumx2*sumuy - sumxy*sumux) + sumx*(sumux*sumy - sumx*sumuy) + sumu*(sumx*sumxy - sumx2*sumy))/delta;
  *ay = (sumv*(sumx2*sumy2 - sumxy*sumxy) + sumx*(sumxy*sumvy - sumvx*sumy2) + sumy*(sumvx*sumxy - sumx2*sumvy))/delta;
  *by = (sum*(sumvx*sumy2 - sumxy*sumvy) + sumv*(sumxy*sumy - sumx*sumy2) + sumy*(sumvy*sumx - sumy*sumvx))/delta;
  *cy = (sum*(sumx2*sumvy - sumxy*sumvx) + sumx*(sumvx*sumy - sumx*sumvy) + sumv*(sumx*sumxy - sumx2*sumy))/delta;
  
}



void syntax(char *name) {
  
  fprintf(stderr,"Syntax: %s templatestarfile imagestarfile\n\n",name);
  fprintf(stderr,"-scale (scaleratio):\n\t ratio of image/template. 0 means figure it out [0]\n");
  fprintf(stderr,"-tol    (scaletolerance):\n\t tolerance in specified scale ratio [0.01]\n");
  fprintf(stderr,"-ntri (matchingtriangles):\n\t  total number of triangles match to declare victory [5]\n");
  fprintf(stderr,"-nmatch (starstomatchatonce):\n\t total number of stars to match at a time [50]\n");
  fprintf(stderr,"-amiss (percentage):\n\t fractional miss allowed in constants [0.05]\n");
  fprintf(stderr,"-cbmiss (percentage):\n\t fractional miss allowed in rotation constants [0.05]\n");
  fprintf(stderr,"-totmiss (percentage):\n\t total sum of errors allowed in rotation matrix [0.1]\n");
  fprintf(stderr,"-roterr (degrees):\n\t angular rotation allowed away from -rotation +/- 0,90,180,270 [3]\n");
  fprintf(stderr,"-rotation (degrees):\n\t expected relative rotation of template and image [0]\n");
  fprintf(stderr,"-squaretol (degrees):\n\t how far away from square allowed [3]\n");
  fprintf(stderr,"-rsummin  (mintrisum):\n\t min allowed obliqueness of triangle [1.1]\n");
  fprintf(stderr,"-tolerance:\n\t fractional tolerance in triangle matching [0.001]\n");
  fprintf(stderr,"-sepmin (pixels):\n\t min length of allowable triangle side [8]\n");
  fprintf(stderr,"-outfile (filename):\n\t file name with output matched list [matchstar.out]\n");
  fprintf(stderr,"-jumpstart3 id1 id2  :\n\t provide program matching stars... can invoke 1-N times\n\tid is 1st col string, or   position (starting with 0) if -noid invoked\n");
  fprintf(stderr,"-noid :\n\t x,y in column 1 and 2, no id column\n");
  fprintf(stderr,"-floormatch :\n\t distance which sigclipping will not go below[0.1]\n");
  fprintf(stderr,"-sigclip :\n\t sigclipping value[2.5]\n");
  fprintf(stderr,"-maxit :\n\t max number of sigclip iterations[5]\n");
  fprintf(stderr,"-outids :\n\t id1,id2 in column 5 and 6 of output\n");
  fprintf(stderr,"-diagnose (0 1 2 3):\n\t how much diagnosis stuff to spew out... [0]\n\n");
  return;
}


#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fitsio.h"
#include <sys/stat.h>
#include <errno.h>
#include "skymap.h"
#include "skymap_cross.h"

#include "biascorr.h"

/* Local function definitions */

int update_header (fitsfile *infptr, fitsfile *outfptr, const int hdupos, const long *outnaxes, int *status);
int update_primary_header (fitsfile *infptr, fitsfile *outfptr, int *status);
int head2range (const char *head, long *range);
int range2head (const long *range, char *head);
int getbiassec(int hdupos, int *bias);
int getcoefflist(int hdupos, double *coeff);
int create_directory(const char *path);
int flag_combine(unsigned char *flag1, unsigned char flag2, char *type);

/*******************************************************************
* FILENAME :   fits64to32_cal.c
*
* DESCRIPTION : 
*       Remove overscan by linear fit or median from one side, combine the amplifiers.
*       Option to align the two amps.
*       Option to flat field.
*       Option to mask saturated pixels and cross-talk.
* Notes :
*
* CHANGES :
* 2012-08-15: add option to flag cross-talk of satuated pixels.
*
**/


int main(int argc, char *argv[])
{
  fitsfile *infptr, *outfptr;  /* FITS file pointers */
  char infile[PATH_MAX], outbase[PATH_MAX], outdir[PATH_MAX], outfile[PATH_MAX];
  char *temp;
  int bias[4], nbias;
  long naxes[2], naxes1;
  long outnaxes[2];
  int hdupos, bitpix, naxis, fpixel0, hdupos1;
  int status = 0;
  int ii, jj;
  float array[1124], newarray[1024], meanbias;
  double xx[1024],fitxx[100],fitarray[100];
  
  long fpixel[2];
  double nulval=0.;
  int anynul;
  int doflat=0;
  int doscale=0;
  int isempty=0;

  fitsfile *flatptr;
  char flatbase[PATH_MAX],flatfile[PATH_MAX],flatpath[PATH_MAX];
  long flatnaxes[2];
  float flatarray[2048],outarray[2048];

  int domask=0;
  fitsfile *maskptr, *inmaskptr;
  char inmaskbase[PATH_MAX],inmaskpath[PATH_MAX],inmaskfile[PATH_MAX],outmaskfile[PATH_MAX];
  int mbitpix;
  long masknaxes[2];
  unsigned char maskarray[2048];
  int cind;
  float satlimit=55000.;

  /* correct cross talk */
  double coeff[3];
  fitsfile *fsatfptr;
  char fsatmap[PATH_MAX];
  int hduother;
  int correctcross=0;
  float satarray[1124];
  int tempjj;
  long nsatpix,nsatarray[32];

  /* RS 2012/04/18: for calculating imagewide median */
  long nmed;
  float medarray[2048*4096];
 
  float medl, medr, meansky, scalel, scaler;
  float checkarray[40960];
  long inc[2]={1,1};
  long lpixel[2],ncheck1,nmid1,ncheck2,nmid2;

  unsigned char flag_good=1, flag_satur=0, flag_xtalk=0;
  char flag_type[]="AND";
  
  /* RS 2012/04/20:  updated the usage message to be more informative */
  if (argc < 3) { 
    printf("Usage: fits64to32_cal infile outdir [-flatpath flatpath -flatbase flatbase] -align -domask [-maskpath maskpath -maskbase maskbase]\n");
    printf("   infile     the input *mosaic* image to split\n");
    printf("   outdir     the output directory into which to place the split images\n");
    printf("   -flatpath  the desired flats are stored in subdirectories of this directory\n");
    printf("              corresponding to each CCD (<flatpath>/01 through <flatpath>/32)\n");
    printf("   -flatbase  the form of the flat filename; e.g. for 'flat_latest_g_27.fits',\n");
    printf("              flatbase = 'flat_latest_g'\n");
    printf("   -align     (boolean) offset the background levels of L + R amps to match?\n");
    printf("   -domask    (boolean) flag satuated pixels and cross-talk affected pixels\n");
    printf("   -maskpath  optional, path of input bad pixel mask\n");
    printf("   -maskbase  optional, name of input bad pixel mask\n");
    printf("   -fsatmap    optional, name of input map for fully saturated pixels\n");
    printf("   -good      flag value for good pixel (default=1)\n");
    printf("   -bad_satur flag value for saturated pixel (default=0)\n");
    printf("   -bad_xtalk flag value for cross-talk affected pixel (default=0)\n");
    printf("   -flag_type how the flags are combinationed, can be AND or OR (default=AND)\n");
    printf("\n");
    return(2);
  }  
  /*get output file name*/
  strcpy(infile,argv[1]);
  strcpy(outbase,argv[2]);  
  temp=strchr(infile,'/');
  if (temp!=NULL){
    temp=strtok(infile,"/");
    do {
      strcpy(infile,temp);
      temp=strtok(NULL,"/");
    } while (temp!=NULL);
  }
  temp=strstr(infile,".fit");
  if (temp==NULL) {
    printf("input file does not have a proper fits extension name\n");
    return(1);
  }
  *temp = '\0'; /* effectively strip off the .fits extension */
  
  /* RS 2012/04/20:  initialize some of these strings to zero */
  strcpy(flatpath,"");
  strcpy(flatbase,"");
  /* for mask */
  strcpy(inmaskpath,"");
  strcpy(inmaskbase,"");

  
  ii = 3; /* start processing options from here */
  while (ii < argc) {
    /* RS 2012/04/20:  added -flatdir to account for how flats are stored
     * (i.e. split into separate directories for each CCD)
     */
    if (strcmp(argv[ii],"-flatpath")==0) {
      if (++ii < argc) {
        strcpy(flatpath,argv[ii]);
        doflat=1;
      }
      else {
	printf("-flatpath requires an argument\n");
	return 1;
      }
    }
    else if (strcmp(argv[ii],"-flatbase")==0) {
      if (++ii < argc) {
        strcpy(flatbase,argv[ii]);
        doflat=1;
      }
      else {
	printf("-flatbase requires an argument\n");
	return 1;
      }
    }
    else if (strcmp(argv[ii],"-maskpath")==0) {
      if (++ii < argc) {
        strcpy(inmaskpath,argv[ii]);
        domask=1;
      } 
      else {
	printf("-maskpath requires an argument\n");
	return 1;
      }
    }
    else if (strcmp(argv[ii],"-maskbase")==0) {
      if (++ii < argc) {
        strcpy(inmaskbase,argv[ii]);
        domask=1;
      }
      else {
	printf("-maskbase requires an argument\n");
	return 1;
      }
    }
    else if (strcmp(argv[ii],"-align")==0) {
      doscale=1;
    }
    else if (strcmp(argv[ii],"-domask")==0) {
      domask=1;
    }
    else if (strcmp(argv[ii],"-fsatmap")==0) {
      if (++ii < argc) {
        strcpy(fsatmap,argv[ii]);
        correctcross=1;
      }
    }
    else if (strcmp(argv[ii],"-good")==0) {
      if (++ii < argc) {
	flag_good=(unsigned char)atoi(argv[ii]);
      }
      else {
	printf("-good requires an argument\n");
	return 1;
      }
    }
    else if (strcmp(argv[ii],"-bad_xtalk")==0) {
      if (++ii < argc) {
	flag_xtalk=(unsigned char)atoi(argv[ii]);
      }
      else {
	printf("-bad_xtalk requires an argument\n");
	return 1;
      }
    }
    else if (strcmp(argv[ii],"-bad_satur")==0) {
      if (++ii < argc) {
	flag_satur=(unsigned char)atoi(argv[ii]);
      }
      else {
	printf("-bad_satur requires an argument\n");
	return 1;
      }
    }
    else if (strcmp(argv[ii],"-flag_type")==0) {
      if (++ii < argc) {
        strcpy(flag_type,argv[ii]);
	if (strcmp(flag_type,"AND")!=0 && strcmp(flag_type,"OR")!=0) {
	  printf("-flag_type must be either AND or OR\n");
	  return 1;
	}
      }
      else {
        printf("-flag_type requires an argument\n");
        return 1;
      }
    }
    ++ii;
  }

  /* RS:  make sure both, or neither, of flatpath and flatbase are defined */
  if (doflat && (strlen(flatpath) == 0 || strlen(flatbase) == 0)) {
    printf("both -flatpath and -flatbase must be defined to do flats\n");
    return(1);
  }
  /* FY: same for maskpath and maskbase (but can have domask without input mask) */
  if ((strlen(inmaskpath) == 0) != (strlen(inmaskbase) == 0)) {
    printf("none or both of maskpath and maskbase have to be defined\n");
    return(1);
  }
  
  naxes1=trimsec[1]-trimsec[0]+1;
  outnaxes[0]=naxes1*2;
  
  if (create_directory(outbase) != 0) {
    fprintf(stderr, "Failed to make directory %s\n", outbase);
    return 3;
  }

  if (domask !=0) {
    /* Open the input file */
    fits_open_file(&infptr, argv[1], READONLY, &status);
    if (correctcross !=0) {
      /* Fully saturated pixels for corrections of ringing */
      fits_open_file(&fsatfptr, fsatmap, READONLY, &status); 
    }
    if (status != 0) {
      fits_report_error(stderr, status);
      return(status);
    }
    
    for (hdupos=1; hdupos<=32; hdupos++) {
      sprintf(outdir,"%s/%02d",outbase,hdupos);
      if (create_directory(outdir) != 0) {
	fprintf(stderr, "Failed to make directory %s\n", outbase);
	return 3;
      }
      
      if (hdupos<=16) {hdupos1=hdupos*2+1;} else {hdupos1=hdupos*2;}
      fits_movabs_hdu(infptr, hdupos1, NULL, &status);
      fits_get_hdu_num(infptr,&hdupos1);
      if (status) break;
      /* check input image size */
      fits_get_img_param(infptr, 2, &bitpix, &naxis, naxes, &status);
      outnaxes[1]=naxes[1];

      /* open input mask */
      if (strlen(inmaskpath) >0) {
	sprintf(inmaskfile,"%s/%02d/%s_%d.fits",inmaskpath,hdupos,inmaskbase,hdupos);
	printf("using input mask %s\n",inmaskfile);
	fits_open_file(&inmaskptr, inmaskfile, READONLY, &status);
	if (status) {
	  printf("Error opening input mask file %s\n", inmaskfile);
	  return(1);
	}
	fits_get_img_param(inmaskptr, 2, &mbitpix, &naxis, masknaxes, &status);
	if (masknaxes[0] != outnaxes[0] || masknaxes[1] !=outnaxes[1]) {
	  printf("Incompatible size of input mask and output image\n");
	  fits_close_file(inmaskptr, &status);
	  return(1);
	}
      }
      /* create output mask */
      sprintf(outmaskfile,"%s/%s_%d_mask.fits",outdir,infile,hdupos);
      printf("writing pixel mask to %s\n",outmaskfile);
      remove(outmaskfile);
      fits_create_file(&maskptr, outmaskfile, &status);
      fits_create_img(maskptr, BYTE_IMG, naxis, outnaxes, &status);
      /* update header */
      fits_write_key_lng(maskptr,"FLAGGD",flag_good,
			 "Flag for good pixel",&status);
      fits_write_key_lng(maskptr,"FLAGST",flag_satur,
			 "Flag for saturated pixel",&status);
      fits_write_key_lng(maskptr,"FLAGXT",flag_xtalk,
			 "Flag for cross-talk affected pixel",&status);
      fits_write_key_flt(maskptr,"SATMASK",satlimit,-6,
			 "Raw count limit for masking saturated and cross-talk",&status);    
      if (strlen(inmaskpath) >0) {
	fits_write_key_str(maskptr,"INMASK",inmaskbase,
			   "Input (bad pixel) mask used",&status);
      }
      /* prefill output mask */
      fpixel[0] = 1;
      for (ii = 1; ii <= outnaxes[1]; ii++) {
	fpixel[1] = ii;
	if (strlen(inmaskpath) >0) {
	  fits_read_pix(inmaskptr, TBYTE, fpixel, outnaxes[0], 
			&nulval, maskarray, &anynul, &status);
	} else {
	  for (jj = 0; jj < outnaxes[0]; jj++){ 
	    maskarray[jj] = flag_good;
	  }
	}
	fits_write_pix(maskptr,TBYTE,fpixel,outnaxes[0],maskarray,&status);
      }
      /* close files */      
      if (strlen(inmaskpath) >0) {
	fits_close_file(inmaskptr, &status);
      }
      fits_close_file(maskptr, &status);
      
      nsatpix = 0;
      fits_open_file(&maskptr, outmaskfile, READWRITE, &status);
      if (correctcross !=0) {
	fits_movabs_hdu(fsatfptr, hdupos1, NULL, &status);
      }
      for (ii = 1; ii <= naxes[1]; ii++) {
	fpixel[0] = 1;
	fpixel[1] = ii;
	fits_read_pix(infptr, TFLOAT, fpixel, naxes[0], &nulval, array,
		      &anynul, &status);
	fits_read_pix(maskptr, TBYTE, fpixel, outnaxes[0], &nulval, maskarray,
		      &anynul, &status);
	if (correctcross !=0) {
	  fits_read_pix(fsatfptr, TFLOAT, fpixel, naxes[0], &nulval, satarray,
		      &anynul, &status);
	}
	for (jj = 0; jj <naxes1; jj++) {
	  if (array[jj+trimsec[0]] >= satlimit) {
	    if (flag_combine(&maskarray[jj],flag_satur,flag_type)!=0){
	      return 1;
	    }
	    /* mask ringing next to saturated pixels*/
	    if (correctcross !=0) {
	      if (satarray[jj+trimsec[0]]>0) {
		flag_combine(&maskarray[jj+1],flag_satur,flag_type);
		flag_combine(&maskarray[jj+2],flag_satur,flag_type);
		flag_combine(&maskarray[jj+3],flag_satur,flag_type);
	      }}
	    else{
	      flag_combine(&maskarray[jj+1],flag_satur,flag_type);
	      flag_combine(&maskarray[jj+2],flag_satur,flag_type);
	      flag_combine(&maskarray[jj+3],flag_satur,flag_type);
	    }	      
	    cind=2047-jj;
	    flag_combine(&maskarray[cind],flag_xtalk,flag_type);
	    nsatpix ++;
	  }
	}
	fits_write_pix(maskptr,TBYTE,fpixel,outnaxes[0],maskarray,&status);
      }
      fits_close_file(maskptr, &status);
      
      if (hdupos<=16) {hdupos1=hdupos*2;} else {hdupos1=hdupos*2+1;}
      fits_movabs_hdu(infptr, hdupos1, NULL, &status);
      fits_get_hdu_num(infptr,&hdupos1);
      if (status) break;
      /* check input image size */
      fits_get_img_param(infptr, 2, &bitpix, &naxis, naxes, &status);
      if (naxes[1] !=outnaxes[1]) {
	printf("Image size is not consistent with the previous amp\n");
	return(1);
      }

      if (correctcross !=0) {
	fits_movabs_hdu(fsatfptr, hdupos1, NULL, &status);
      }
      fits_open_file(&maskptr, outmaskfile, READWRITE, &status);
      for (ii = 1; ii <= naxes[1]; ii++) {
        fpixel[0] = 1;
        fpixel[1] = ii;
        fits_read_pix(infptr, TFLOAT, fpixel, naxes[0], &nulval, array,
                      &anynul, &status);
        fits_read_pix(maskptr, TBYTE, fpixel, outnaxes[0], &nulval, maskarray,
                      &anynul, &status);
	if (correctcross !=0) {
	  fits_read_pix(fsatfptr, TFLOAT, fpixel, naxes[0], &nulval, satarray,
			&anynul, &status);
	}
        for (jj = 0;jj <naxes1;jj++) {
          if (array[jj+trimsec[0]]>=satlimit) {
            flag_combine(&maskarray[jj+naxes1],flag_satur,flag_type);
	    if (correctcross !=0) {
	      if (satarray[jj+trimsec[0]]>0) {
		flag_combine(&maskarray[jj+naxes1-1],flag_satur,flag_type);
		flag_combine(&maskarray[jj+naxes1-2],flag_satur,flag_type);
		flag_combine(&maskarray[jj+naxes1-3],flag_satur,flag_type);
	      }}
	    else{
		flag_combine(&maskarray[jj+naxes1-1],flag_satur,flag_type);
		flag_combine(&maskarray[jj+naxes1-2],flag_satur,flag_type);
		flag_combine(&maskarray[jj+naxes1-3],flag_satur,flag_type);
	    }	      
            cind=2047-(jj+naxes1);
            flag_combine(&maskarray[cind],flag_xtalk,flag_type);
	    nsatpix ++;
	  }
	}
        fits_write_pix(maskptr,TBYTE,fpixel,outnaxes[0],maskarray,&status);
      }
      fits_write_key_lng(maskptr,"NSATPIX",nsatpix,
			 "Number of saturated pixels",&status);
      nsatarray[hdupos-1]=nsatpix;
      fits_close_file(maskptr, &status);
    }
    fits_close_file(infptr, &status);
    if (correctcross !=0) {
      fits_close_file(fsatfptr, &status); 
    }
    
  }
  
  /* Overscan correction and optional cross talk correction */
  fits_open_file(&infptr, argv[1], READONLY, &status);
  if (correctcross !=0) {
    fits_open_file(&fsatfptr, fsatmap, READONLY, &status); 
  }
  if (status != 0) {    
    fits_report_error(stderr, status);
    return(status);
  }
  for (ii = 0; ii <naxes1; ii++){
    xx[ii] = ii+trimsec[0];
  }
  for (hdupos=1; hdupos<=32; hdupos++) {
    sprintf(outdir,"%s/%02d",outbase,hdupos);
    if (create_directory(outdir) != 0) {
      fprintf(stderr, "Failed to make directory %s\n", outbase);
      return 3;
    }
    sprintf(outfile,"%s/%s_%d.fits",outdir,infile,hdupos);
    printf("writing %s\n",outfile);
    
    if (hdupos<=16) {hdupos1=hdupos*2+1;} else {hdupos1=hdupos*2;}
    getbiassec(hdupos1,bias);
    nbias = (bias[3]-bias[2]+1) + (bias[1]-bias[0]+1);
    
    fits_movabs_hdu(infptr, hdupos1, NULL, &status);
    if (status) break;
    /* check input image size */
    fits_get_img_param(infptr, 2, &bitpix, &naxis, naxes, &status);
    if (naxes[0] < bias[3]) {
      printf("Image size is smaller than defined bias region \n");
      return(1);
    }        
    
    outnaxes[1]=naxes[1];
    remove(outfile);
    fits_create_file(&outfptr, outfile, &status);
    fits_create_img(outfptr, FLOAT_IMG, naxis, outnaxes, &status);
    fits_movabs_hdu(infptr, 1, NULL, &status);
    update_primary_header(infptr,outfptr,&status);
    fits_close_file(outfptr, &status);
 
    fits_movabs_hdu(infptr, hdupos1, NULL, &status);
    fits_get_hdu_num(infptr,&hdupos1);
    /* printf("input extension %d\n",hdupos1-1); */
    fits_open_file(&outfptr, outfile, READWRITE, &status);
    update_header(infptr,outfptr,hdupos,outnaxes,&status);
    fpixel0=1;
    if (biasmethod[hdupos1-2]==1) {
      meanbias=biascorr_side(infptr,outfptr,1,bias,trimsec,naxes,naxes1,fpixel0,xx,
                             fitxx,fitarray,array,newarray,&status);
    } else {
      meanbias=biascorr(infptr,outfptr,bias,trimsec,naxes,naxes1,fpixel0,xx,
                        fitxx,fitarray,array,newarray,&status);
    }
    fits_write_key_flt(outfptr,"BIASL",meanbias,-5,
                       "Mean bias removed at center column",&status);
    if (correctcross !=0) {
      /*correct cross talk here*/
      if (hdupos<=16) {hduother=hdupos*2;} else {hduother=hdupos*2+1;}
      fits_movabs_hdu(infptr, hduother, NULL, &status);
      fits_movabs_hdu(fsatfptr, hduother, NULL, &status);
      getcoefflist(hduother,coeff);
      for (ii = 1; ii <= naxes[1]; ii++) {
	fpixel[0] = 1;
	fpixel[1] = ii;
	fits_read_pix(infptr, TFLOAT, fpixel, naxes[0], &nulval, array,
		      &anynul, &status);
	fits_read_pix(fsatfptr, TFLOAT, fpixel, naxes[0], &nulval, satarray,
		      &anynul, &status);
	fits_read_pix(outfptr, TFLOAT, fpixel, naxes1, &nulval, newarray,
		      &anynul, &status);
	for (jj = 0;jj <naxes1;jj++) {
	  tempjj=1023-jj+trimsec[0];
	  newarray[jj]-=coeff[0]*pow(array[tempjj],2)+coeff[1]*array[tempjj];
	  /*correct cross talk for fully saturated pixels*/
	  if (satarray[tempjj]>0) {
	    newarray[jj]-=coeff[2];
	  }
	}
	fits_write_pix(outfptr, TFLOAT, fpixel, naxes1, newarray,&status);
      }
      fits_write_key_flt(outfptr,"CROSS_L1",coeff[1],3,"First order cross-talk correction",&status);
      fits_write_key_flt(outfptr,"CROSS_L2",coeff[0],3,"Second order cross-talk correction",&status);
      fits_write_key_flt(outfptr,"CROSS_L0",coeff[2],3,"Addtional correction for full-well saturation",&status);
    }
    fits_close_file(outfptr, &status);
    
    if (hdupos<=16) {hdupos1=hdupos*2;} else {hdupos1=hdupos*2+1;}
    getbiassec(hdupos1,bias);
    nbias = (bias[3]-bias[2]+1) + (bias[1]-bias[0]+1);
    
    fits_movabs_hdu(infptr, hdupos1, NULL, &status);
    fits_get_hdu_num(infptr,&hdupos1);
    /* printf("input extension: %d\n",hdupos1-1); */
    if (status) break;
    /* check input image size */
    fits_get_img_param(infptr, 2, &bitpix, &naxis, naxes, &status);
    if (naxes[0] < bias[3]) {
      printf("Image size is smaller than defined bias region\n");
      return(1);
    }        
    if (naxes[1] !=outnaxes[1]) {
      printf("Image size is not consistent with the previous amp\n");
      return(1);
    }
    
    fits_open_file(&outfptr, outfile, READWRITE, &status);
    fpixel0=naxes1+1;
    if (biasmethod[hdupos1-2]==1) {
      meanbias=biascorr_side(infptr,outfptr,0,bias,trimsec,naxes,naxes1,fpixel0,xx,
                             fitxx,fitarray,array,newarray,&status);
    } else {
      meanbias=biascorr(infptr,outfptr,bias,trimsec,naxes,naxes1,fpixel0,xx,
                        fitxx,fitarray,array,newarray,&status);
    }    
    fits_write_key_flt(outfptr,"BIASR",meanbias,-5,
                       "Mean bias removed at center column",&status);    
    if (correctcross !=0) {
      /*correct cross talk here*/
      if (hdupos<=16) {hduother=hdupos*2+1;} else {hduother=hdupos*2;}
      fits_movabs_hdu(infptr, hduother, NULL, &status);
      fits_movabs_hdu(fsatfptr, hduother, NULL, &status);
      getcoefflist(hduother,coeff);
      for (ii = 1; ii <= naxes[1]; ii++) {
	fpixel[0] = 1;
	fpixel[1] = ii;
	fits_read_pix(infptr, TFLOAT, fpixel, naxes[0], &nulval, array,
		      &anynul, &status);
	fits_read_pix(fsatfptr, TFLOAT, fpixel, naxes[0], &nulval, satarray,
		      &anynul, &status);
	fpixel[0]=fpixel0;
	fits_read_pix(outfptr, TFLOAT, fpixel, naxes1, &nulval, newarray,
		      &anynul, &status);
	for (jj = 0;jj <naxes1;jj++) {
	  tempjj=1023-jj+trimsec[0];
	  newarray[jj]-=coeff[0]*pow(array[tempjj],2)+coeff[1]*array[tempjj];
	  /*correct cross talk for fully saturated pixels*/
	  if (satarray[tempjj]>0) {
	    newarray[jj]-=coeff[2];
	  }
	}
	fits_write_pix(outfptr, TFLOAT, fpixel, naxes1, newarray,&status);
      }
      fits_write_key_flt(outfptr,"CROSS_R1",coeff[1],3,"First order cross-talk correction",&status);
      fits_write_key_flt(outfptr,"CROSS_R2",coeff[0],3,"Second order cross-talk correction",&status);
      fits_write_key_flt(outfptr,"CROSS_R0",coeff[2],3,"Addtional correction for full-well saturation",&status);
    }
    /* update header */
    fits_write_key_flt(outfptr,"SATMASK",satlimit,-6,
		       "Raw count limit for masking saturated and cross-talk",&status);    
    fits_write_key_lng(outfptr,"NSATPIX",nsatarray[hdupos-1],
		       "Number of saturated pixels",&status);
    fits_close_file(outfptr, &status);
    
    
    if (doflat !=0) {
      /* apply flat field */
      /* RS 2012/04/02:  Added some other sanity checks here. */
      sprintf(flatfile,"%s/%02d/%s_%d.fits",flatpath,hdupos,flatbase,hdupos);
      printf("applying flat %s\n",flatfile);
      /* RS:  first make sure the flat is really there and we can open it! */
      fits_open_file(&flatptr, flatfile, READONLY, &status);
      if (status)
      {
         printf("Error opening flat field file %s\n", flatfile);
         return(1);
      }
      /* RS:  now find its shape */
      fits_get_img_param(flatptr, 2, &bitpix, &naxis, flatnaxes, &status);
      if (flatnaxes[0] != outnaxes[0] || flatnaxes[1] !=outnaxes[1]) {
        printf("Incompatible size of flat and output image\n");
        fits_close_file(flatptr, &status);
        return(1);
      }   
      fits_open_file(&outfptr, outfile, READWRITE, &status);
      for (ii = 1; ii <= outnaxes[1]; ii++) {
        fpixel[0] = 1;
        fpixel[1] = ii;
        fits_read_pix(outfptr, TFLOAT, fpixel, outnaxes[0], &nulval, outarray,
                      &anynul, &status);
        fits_read_pix(flatptr, TFLOAT, fpixel, outnaxes[0], &nulval, flatarray,
                      &anynul, &status);
        for(jj=0; jj< outnaxes[0]; jj++) {
          if (flatarray[jj] !=0.)
            outarray[jj] /= flatarray[jj];
          else
            outarray[jj] = 0.;
        }
        fits_write_pix(outfptr, TFLOAT, fpixel, outnaxes[0], outarray,
                       &status);
      }
      fits_write_key_str(outfptr,"FLAT",flatbase,
                         "flat field applied",&status);
      fits_close_file(outfptr, &status);
      fits_close_file(flatptr, &status);
    }

    if (doscale !=0) {
      /*align the right part to the left */
      /*printf("aligning right and left parts to middle...\n"); */
      ncheck1=(checksec[1]-checksec[0]+1)*outnaxes[1];
      nmid1=ncheck1/2;
      ncheck2=(checksec[3]-checksec[2]+1)*outnaxes[1];
      nmid2=ncheck2/2;
      fits_open_file(&outfptr, outfile, READWRITE, &status);
      
      lpixel[1]=outnaxes[1];   
      fpixel[0]=checksec[0];
      lpixel[0]=checksec[1];
      fpixel[1]=1;
      fits_read_subset(outfptr, TFLOAT, fpixel, lpixel, inc,
                       &nulval, checkarray, &anynul, &status);
      medl=kth_smallest_flt(checkarray,ncheck1,nmid1);
      
      fpixel[0]=checksec[2];
      lpixel[0]=checksec[3];
      fits_read_subset(outfptr, TFLOAT, fpixel, lpixel, inc,
                       &nulval, checkarray, &anynul, &status);
      medr=kth_smallest_flt(checkarray,ncheck2,nmid2);
      
      /*update header */
      if ((medl <10) || (medr < 10)) {
        printf("mean sky count is too low for scaling to be accurate, so no aligning\n");
	if ((medl <0.001) || (medr < 0.001)) {
	  isempty=1;
	}
      }
      else {
        meansky=(medl+medr)/2;
        scalel=meansky/medl;
        scaler=meansky/medr;
        fits_write_key_flt(outfptr,"AMPSCL_L",scalel,3,
                           "Left amplifier adjusted by",&status);    
        fits_write_key_flt(outfptr,"AMPSCL_R",scaler,3,
                           "Right amplifier adjusted by",&status);    

        /*rescale*/
        fpixel[0]=1;
        for (ii=1;ii<=outnaxes[1];ii++) {
          fpixel[1]=ii;
          fits_read_pix(outfptr,TFLOAT,fpixel,naxes1,&nulval,newarray,&anynul,&status);
          for (jj=0;jj<naxes1;jj++){
            newarray[jj] *=scalel;
          }
          fits_write_pix(outfptr,TFLOAT,fpixel,naxes1,newarray,&status);
        }
        fpixel[0]=naxes1+1;
        for (ii=1;ii<=outnaxes[1];ii++) {
          fpixel[1]=ii;
          fits_read_pix(outfptr,TFLOAT,fpixel,naxes1,&nulval,newarray,&anynul,&status);
          for (jj=0;jj<naxes1;jj++){
            newarray[jj] *=scaler;
	  }
          fits_write_pix(outfptr,TFLOAT,fpixel,naxes1,newarray,&status);
        }
        fits_write_history(outfptr, "CCD combined with scaling", &status);
	fits_close_file(outfptr, &status);
      }
    }

    /* RS:  calculate sky median and write out to the header.
     * Note:  this implicitly assumes the image isn't a crowded field.
     * But since the SN survey is avoiding the galaxy and we'll only ever
     * use this keyword for things like sky flats, it's probably ok.
     * FY: probably consider taking the mean of left and right?
     */
    fits_open_file(&outfptr, outfile, READWRITE, &status);
    nmed=outnaxes[0]*outnaxes[1];
    lpixel[1]=outnaxes[1];   
    lpixel[0]=outnaxes[0];
    fpixel[1]=1;
    fpixel[0]=1;
    fits_read_subset(outfptr, TFLOAT, fpixel, lpixel, inc,
		     &nulval, medarray, &anynul, &status);
    meansky=kth_smallest_flt(medarray,nmed,nmed/2);
    fits_write_key_flt(outfptr,"MEDPIX",meansky,-8,
		       "Median bias-subtracted pixel value",&status);
    fits_write_key_flt(outfptr,"FLATNORM",meansky,5,
		       "SWarp flux scaling for flat field = 1/FLATNORM",&status);
    
    fits_close_file(outfptr, &status);
  }
  fits_close_file(infptr, &status);
  if (correctcross !=0) {
    fits_close_file(fsatfptr, &status); 
  }
  
  if (isempty==0) {
    /* if error occurred, print out error message */
    if (status)
      fits_report_error(stderr, status);
    return(status);
  } else {
    return(10);
  }
  
 }

int update_primary_header (fitsfile *infptr, fitsfile *outfptr, int *status) {
  int ii, nkeys;
  char card[FLEN_CARD];
  
  fits_get_hdrspace(infptr, &nkeys, NULL, status); 
  for (ii = 1; ii <= nkeys; ii++) {
    fits_read_record(infptr, ii, card, status);
    if (fits_get_keyclass(card) > TYP_SCAL_KEY)
      fits_write_record(outfptr, card, status);
  }
  fits_write_date(outfptr, status);
  fits_delete_key(outfptr,"NEXTEND", status);
  fits_update_key_log(outfptr,"EXTEND",0,NULL, status);

  return *status;
} 

int update_header (fitsfile *infptr, fitsfile *outfptr, const int hdupos, const long *outnaxes, int *status) {
  int ii, nkeys;
  char card[FLEN_CARD], keyvaluestr[FLEN_VALUE];
  long atm22;
  /*long dtv1, dtv2, detsec[4];
    int xccd, yccd;*/
  long ccdsec[4], temp;
  
  fits_get_hdrspace(infptr, &nkeys, NULL, status); 
  for (ii = 1; ii <= nkeys; ii++) {
    fits_read_record(infptr, ii, card, status);
    if (fits_get_keyclass(card) > TYP_SCAL_KEY)
      fits_write_record(outfptr, card, status);
  }
  
  /* update extension name and number */
  sprintf(keyvaluestr,"im%d",hdupos);
  fits_update_key_str(outfptr,"EXTNAME",keyvaluestr,NULL,status);
  fits_update_key_lng(outfptr,"EXTVER",hdupos,NULL,status);
  fits_update_key_lng(outfptr,"IMAGEID",hdupos,NULL,status);
  
  /* update CCDSIZE LTV1 TRIMSEC DATASEC */
  sprintf(keyvaluestr,"[1:%ld,1:%ld]",outnaxes[0],outnaxes[1]);
  fits_update_key_str(outfptr,"CCDSIZE",keyvaluestr,NULL,status);
  fits_update_key_flt(outfptr,"LTV1",0.,6,NULL,status);
  fits_update_key_str(outfptr,"TRIMSEC",keyvaluestr,NULL,status);
  fits_update_key_str(outfptr,"DATASEC",keyvaluestr,NULL,status);
  
  /* delete AMPNAME BIAS0001 BIASSEC */
  fits_delete_key(outfptr,"AMPNAME",status);
  fits_delete_key(outfptr,"BIAS0001",status);
  fits_delete_key(outfptr,"BIASSEC",status);
  
  /*update ATV1 ATM1_1*/
  fits_update_key_lng(outfptr,"ATV1",0,NULL,status);
  fits_update_key_lng(outfptr,"ATM1_1",1,NULL,status);
  
  fits_delete_key(outfptr,"DTV1",status);
  fits_delete_key(outfptr,"DTV2",status);
  fits_delete_key(outfptr,"DETSEC",status);
  fits_delete_key(outfptr,"DTM1_1",status);
  fits_delete_key(outfptr,"DTM1_2",status);
  fits_delete_key(outfptr,"DTM2_1",status);
  fits_delete_key(outfptr,"DTM2_2",status);
  /*update CCDSEC AMPSEC ATV2 ATM2_2 
    DETSEC DTV1 DTV2 */
  /*fits_read_key_lng(outfptr,"DTV1",&dtv1,NULL,&status);
    fits_read_key_lng(outfptr,"DTV2",&dtv2,NULL,&status);*/
  fits_read_key_str(outfptr,"CCDSEC",keyvaluestr,NULL,status);
  
  head2range(keyvaluestr,ccdsec);
  /*xccd=(dtv1+ccdsec[0])/2248l;
  yccd=(dtv2+ccdsec[2])/4096l;
  dtv1=xccd*outnaxes[0]+xccd*ccdgaps[0];
  dtv2=yccd*outnaxes[1]+yccd*ccdgaps[1]-(ccdsec[2]-1);
  if (yccd >= 2) dtv2+=(ccdgaps[2]-ccdgaps[1]-1);
  */
  ccdsec[0]=1;
  ccdsec[1]=outnaxes[0];
  range2head(ccdsec,keyvaluestr);
  fits_update_key_str(outfptr,"CCDSEC",keyvaluestr,NULL,status);
  /*fits_update_key_lng(outfptr,"DTV1",dtv1,NULL,status);
    fits_update_key_lng(outfptr,"DTV2",dtv2,NULL,status);*/
  fits_update_key_str(outfptr,"AMPSIZE",keyvaluestr,NULL,status);
  fits_update_key_str(outfptr,"AMPSEC",keyvaluestr,NULL,status);

  /*detsec[0]=dtv1+ccdsec[0];
  detsec[1]=dtv1+ccdsec[1];
  detsec[2]=dtv2+ccdsec[2];
  detsec[3]=dtv2+ccdsec[3];
  range2head(detsec,keyvaluestr);
  fits_update_key_str(outfptr,"DETSEC",keyvaluestr,NULL,status);*/
 
  fits_read_key_lng(outfptr,"ATM2_2",&atm22,NULL,status);
  if (atm22 <0 ){
    temp=ccdsec[3];
    ccdsec[3]=ccdsec[2];
    ccdsec[2]=temp;
    range2head(ccdsec,keyvaluestr);
    fits_update_key_str(outfptr,"AMPSEC",keyvaluestr,NULL,status);
  }  
  return *status;
}

int head2range (const char *head, long *range) {

  char *p,*p2;
  char *head_copy = (char*)malloc(strlen(head) + 1);
  
  if (head_copy != NULL) {
    strcpy(head_copy, head);
    
    p=strtok(head_copy,"[");
    p2=strtok(p,":");
    range[0]=atol(p2);
    p2=strtok(NULL,",");
    range[1]=atol(p2);
    p2=strtok(NULL,":");
    range[2]=atol(p2);
    p2=strtok(NULL,"]");
    range[3]=atol(p2);

    free(head_copy);
  }
  return 0;
}

int range2head (const long *range, char *head) {
  sprintf(head,"[%ld:%ld,%ld:%ld]",range[0],range[1],range[2],range[3]);
  return 0;
}

int getbiassec(int hdupos, int *bias) {
  int ii,start;
  start=(hdupos-2)*4;
  for (ii=0;ii<4;ii++) {
    bias[ii]=biassec[start+ii];
  }
  return 0;
}

int getcoefflist(int hdupos, double *coeff) {
  int ii,start;
  start=(hdupos-2)*3;
  for (ii=0;ii<3;ii++) {
    coeff[ii]=cross_coeff[start+ii];
  }
  return 0;
}

/* Returns 0 if the directory already exists or was created */
int create_directory(const char *path)
{
  struct stat statbuf;

  if (stat(path,&statbuf) != 0) {
    /* does not exist - create it */
    return mkdir(path, 0777);
  }
  else {
    /* check it is a directory */
    if (!S_ISDIR(statbuf.st_mode))
      return -1;
  }

  return 0; /* already exist is a directory */
}

int flag_combine(unsigned char *flag1, unsigned char flag2, char *type) {
  if (strcmp(type,"AND")==0 || strcmp(type,"and")==0) {

    *flag1= (unsigned char)(*flag1 & flag2);
    return 0;
  }
  else if (strcmp(type,"OR")==0 || strcmp(type,"or")==0) {
    *flag1= (unsigned char)(*flag1 | flag2);
    return 0;
  }
  else {
    fprintf(stderr, "Combine type %s is not valid.\n", type);
    return 1;
  }
}

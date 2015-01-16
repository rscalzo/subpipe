/*
 
$Id: cir_wcssubs.c,v 1.14 2005/03/08 20:05:46 jim Exp $
 
*/

#include <stdio.h>
#include <stdlib.h>

#include <cirdr.h>
#include <cir_wcssubs.h>
#include <math.h>

static void fk425pv (double *, double *, double *, double *, double *, 
		     double *);
/*+
 *  Name:
 *      cir_wcsopen
 *
 *  Purpose:
 *      Opens a WCS structure.
 *
 *  Description:
 *      Reads the header of an input images and constructs a WCS structure
 *      from it.
 *
 *  Language:
 *      C
 *
 *  Arguments:
 *      infile = char * (Given)
 *          The name of the input FITS image. This must have a valid WCS
 *          defined in the header
 *      wcs = struct wcsprm ** (Returned)
 *          The WCSLIB structure defined from the header
 *      errmsg = char * (Returned)
 *          Error message generated if something went wrong during the
 *          processing.
 *
 *  Returned values:
 *      Standard cirdr status values (see cirdr.h)
 *
 *  Notes:
 *      None
 *
 *  Dependencies:
 *      cfitsio, wcslib
 *
 *  Authors:
 *      Jim Lewis (CASU, IoA)
 *
 *  Copyright:
 *      Copyright (C) 2002-2005 Cambridge Astronomy Survey Unit.
 *      All Rights Reserved.
 *
+*/

extern int cir_wcsopen(char *infile, struct wcsprm **wcs, char *errmsg) {
    int status,nk,retval,nrej,nwcs;
    fitsfile *iptr = NULL;
    char *hdr = NULL,msg[BUFSIZ];
    struct wcsprm *wwcs = NULL;

    /* Open the file and read the header out of it */

    status = 0;
    (void)fits_open_file(&iptr,infile,READONLY,&status);
    (void)fits_hdr2str(iptr,1,NULL,0,&hdr,&nk,&status);
    if (status != 0) {
	(void)fits_get_errstatus(status,msg);
	(void)sprintf(errmsg,"WCSOPEN: Can't read file header %s -- %s\n",
		      infile,msg);
	closefits(iptr);
	freespace(hdr);
	return(CIR_FATAL);
    }
    closefits(iptr);

    /* Now get the WCS structure */

    retval = wcspih(hdr,nk,0,0,&nrej,&nwcs,&wwcs);
    if (retval != 0) {
	(void)sprintf(errmsg,"WCSOPEN: WCSPIH failed with status: %d\n",
		      retval);
	freespace(hdr);
	wcsvfree(&nwcs,&wwcs);
	return(CIR_FATAL);
    }
    freespace(hdr);

    /* Get the one you want (the first one in this case) */

    *wcs = (struct wcsprm *)cir_calloc(1,sizeof(struct wcsprm));
    (*wcs)->flag = -1;
    if (wcscopy(1,wwcs,*wcs) != 0 || wcsset(*wcs) != 0) {
	(void)sprintf(errmsg,"WCSOPEN: WCSCOPY/WCSSET failed\n");
	wcsvfree(&nwcs,&wwcs);
	wcsfree(*wcs);
	return(CIR_FATAL);
    }

    /* Get out of here */

    wcsvfree(&nwcs,&wwcs);
    return(CIR_OK);
}
 
/*+
 *  Name:
 *      cir_wcsclose
 *
 *  Purpose:
 *      Closes a WCS structure.
 *
 *  Description:
 *      Closes a WCS structure and frees all workspace associated with it
 *
 *  Language:
 *      C
 *
 *  Arguments:
 *      wcs = struct wcsprm * (Given)
 *          The input WCS structure pointer
 *
 *  Returned values:
 *      None
 *
 *  Notes:
 *      This is simply a wrapper routine for WCSLIB routine wcsfree.
 *
 *  Dependencies:
 *      wcslib
 *
 *  Authors:
 *      Jim Lewis (CASU, IoA)
 *
 *  Copyright:
 *      Copyright (C) 2002-2005 Cambridge Astronomy Survey Unit.
 *      All Rights Reserved.
 *
+*/
                                                                               
extern void cir_wcsclose(struct wcsprm *wcs) {
    int retval;

    /* Free the workspace */

    retval = wcsfree(wcs);
    free(wcs);
}

/*+
 *  Name:
 *      cir_xytoradec
 *
 *  Purpose:
 *      Convert X,Y coordinates in pixels to RA and Dec in degrees.
 *
 *  Description:
 *      A WCS structure is used to convert cartesian coordinates into
 *      equatorial.  The structure must exist before entering this routine.
 *
 *  Language:
 *      C
 *
 *  Arguments:
 *      wcs = struct wcsprm * (Given)
 *          The input WCS structure pointer
 *      x = double (Given)
 *          The X coordinate to be converted.
 *      y = double (Given)
 *          The Y coordinate to be converted.
 *      ra = double * (Returned)
 *          The calculated RA
 *      dec = double * (Returned)
 *          The calculated Dec.
 *
 *  Returned values:
 *      None
 *
 *  Notes:
 *      This is simply a wrapper routine for WCSLIB routine wcsp2s. The
 *      equinox of the output coordinates are defined by the wcs structure
 *      when it was created.
 *
 *  Dependencies:
 *      wcslib
 *
 *  Authors:
 *      Jim Lewis (CASU, IoA)
 *
 *  Copyright:
 *      Copyright (C) 2002-2005 Cambridge Astronomy Survey Unit.
 *      All Rights Reserved.
 *
+*/

extern void cir_xytoradec(struct wcsprm *wcs, double x, double y, 
			  double *ra, double *dec) {
    double xy[2],std[2],phi,theta,radec[2];
    int stat,retval;

    /* Load up the information and call the wcslib routine */

    xy[0] = x;
    xy[1] = y;
    retval = wcsp2s(wcs,1,2,xy,std,&phi,&theta,radec,&stat);

    /* Pass it back now */

    *ra = radec[0];
    *dec = radec[1];
}

/*+
 *  Name:
 *      cir_radectoxy
 *
 *  Purpose:
 *      Convert RA and Dec in degrees to X and Y in pixels.
 *
 *  Description:
 *      A WCS structure is used to convert equatorial coordiantes to
 *      cartesian  The structure must exist before entering this routine.
 *
 *  Language:
 *      C
 *
 *  Arguments:
 *      wcs = struct wcsprm * (Given)
 *          The input WCS structure pointer
 *      ra = double (Given)
 *          The input RA
 *      dec = double (Given)
 *          The input Dec.
 *      x = double * (Returned)
 *          The output X coordinate
 *      y = double * (Returned)
 *          The output Y coordinate
 *
 *  Returned values:
 *      None
 *
 *  Notes:
 *      This is simply a wrapper routine for WCSLIB routine wcss2p.
 *
 *  Dependencies:
 *      wcslib
 *
 *  Authors:
 *      Jim Lewis (CASU, IoA)
 *
 *  Copyright:
 *      Copyright (C) 2002-2005 Cambridge Astronomy Survey Unit.
 *      All Rights Reserved.
 *
+*/

extern void cir_radectoxy(struct wcsprm *wcs, double ra, double dec,
			  double *x, double *y) {
    double xy[2],std[2],phi,theta,radec[2];
    int stat,retval;

    /* Load up the information and call the wcslib routine */

    radec[0] = ra;
    radec[1] = dec;
    retval = wcss2p(wcs,1,2,radec,&phi,&theta,std,xy,&stat);

    /* Pass it back now */

    *x = xy[0];
    *y = xy[1];
}

/*+
 *  Name:
 *      cir_radectoxieta
 *
 *  Purpose:
 *      Convert RA and Dec in degrees to standard coordinates (xi and eta)
 *      in degrees.
 *
 *  Description:
 *      A WCS structure is used to convert equatorial coordiantes to
 *      standard coordinates.
 *
 *  Language:
 *      C
 *
 *  Arguments:
 *      wcs = struct wcsprm * (Given)
 *          The input WCS structure pointer
 *      ra = double (Given)
 *          The input RA
 *      dec = double (Given)
 *          The input Dec.
 *      xi = double * (Returned)
 *          The output xi coordinate
 *      eta = double * (Returned)
 *          The output eta coordinate
 *
 *  Returned values:
 *      None
 *
 *  Notes:
 *      This is simply a wrapper routine for WCSLIB routine wcss2p.
 *
 *  Dependencies:
 *      wcstools
 *
 *  Authors:
 *      Jim Lewis (CASU, IoA)
 *
 *  Copyright:
 *      Copyright (C) 2002-2005 Cambridge Astronomy Survey Unit.
 *      All Rights Reserved.
 *
+*/

extern void cir_radectoxieta(struct wcsprm *wcs, double ra, double dec,
			     double *xi, double *eta) {
    double xy[2],std[2],phi,theta,radec[2];
    int stat,retval;

    /* Load up the information and call the wcslib routine */

    radec[0] = ra;
    radec[1] = dec;
    retval = wcss2p(wcs,1,2,radec,&phi,&theta,std,xy,&stat);

    /* Pass it back now */

    *xi = DEGRAD*std[0];
    *eta = DEGRAD*std[1];
}

/*+
 *  Name:
 *      cir_fk425
 *
 *  Purpose:
 *      Convert from FK4 to FK5
 *
 *  Description:
 *      Equatorial coordinates in FK4 are converted to FK5.
 *
 *  Language:
 *      C
 *
 *  Arguments:
 *      ra = double * (Given and Returned)
 *          The input RA.  This value is overwritten by the new value.
 *      dec = double * (Given and Returned)
 *          The input Dec. This value is overwritten by the new value.
 *
 *  Returned values:
 *      None
 *
 *  Notes:
 *      This is simply a wrapper routine for the internal routine fk425pv.
 *
 *  Dependencies:
 *      None
 *
 *  Authors:
 *      Jim Lewis (CASU, IoA)
 *
 *  Copyright:
 *      Copyright (C) 2002-2005 Cambridge Astronomy Survey Unit.
 *      All Rights Reserved.
 *
+*/

extern void cir_fk425(double *ra, double *dec) {
    double rapm,decpm,parallax,rv;

    /* Zero out the stuff you don't care about */

    rapm = (double)0.0;
    decpm = (double)0.0;
    parallax = (double)0.0;
    rv = (double)0.0;

    /* Do the work */

    fk425pv(ra,dec,&rapm,&decpm,&parallax,&rv);
}

/*+
 *  Name:
 *      cir_newwcs
 *
 *  Purpose:
 *      Create a WCS structure given some basic WCS info.
 *
 *  Description:
 *      Takes some WCS information and creates a TAN project WCS structure.
 *
 *  Language:
 *      C
 *
 *  Arguments:
 *      wcs = struct wcsprm ** (Returned)
 *          The output WCS structure
 *      crval1 = double (Given)
 *          The WCS value at the reference point on the X axis.
 *      crval2 = double (Given)
 *          The WCS value at the reference point on the Y axis.
 *      crpix1 = double (Given)
 *          The X axis value of the reference point
 *      crpix2 = double (Given)
 *          The Y axis value of the reference point.
 *      cd1_1 = double (Given)
 *          The 1,1 element of the CD matrix
 *      cd1_2 = double (Given)
 *          The 1,2 element of the CD matrix
 *      cd2_1 = double (Given)
 *          The 2,1 element of the CD matrix
 *      cd2_2 = double (Given)
 *          The 2,2 element of the CD matrix
 *      naxis1 = int (Given)
 *          The size of the proposed X axis.
 *      naxis2 = int (Given)
 *          The size of the proposed Y axis.
 *      errmsg = char * (Returned)
 *          Error message generated if something went wrong during the
 *          processing.
 *
 *  Returned values:
 *      Standard cirdr status values (see cirdr.h)
 *
 *  Notes:
 *      A basic TAN projection header is created. There is no option for any
 *      thing more complicated.  Will rethink if the need arises.
 *
 *  Dependencies:
 *      wcslib
 *
 *  Authors:
 *      Jim Lewis (CASU, IoA)
 *
 *  Copyright:
 *      Copyright (C) 2002-2005 Cambridge Astronomy Survey Unit.
 *      All Rights Reserved.
 *
+*/
 
extern int cir_newwcs(struct wcsprm **wcs, double crval1, double crval2, 
		      double crpix1, double crpix2, double cd1_1, double cd1_2,
		      double cd2_1, double cd2_2, int naxis1, int naxis2,
		      char *errmsg) {
    struct wcsprm *wwcs;
    char *hdrstr;
    int i,nwcs,retval,nrej,j;

    /* Get workspace */

    hdrstr = cir_calloc(2880,sizeof(char *)); 

    /* Create data array headers */

    i = 0;
    sprintf(hdrstr+i*80,"NAXIS   = %20d",2);
    i++;
    sprintf(hdrstr+i*80,"NAXIS1  = %20d",naxis1);
    i++;
    sprintf(hdrstr+i*80,"NAXIS2  = %20d",naxis2);

    /* Create the CTYPE headers */

    i++;
    sprintf(hdrstr+i*80,"CTYPE1  = 'RA---ZPN'");
    i++;
    sprintf(hdrstr+i*80,"CTYPE2  = 'DEC--ZPN'");

    /* Create CRVAL and CRPIX headers */

    i++;
    sprintf(hdrstr+i*80,"CRVAL1  = %20g",crval1);
    i++;
    sprintf(hdrstr+i*80,"CRVAL2  = %20g",crval2);
    i++;
    sprintf(hdrstr+i*80,"CRPIX1  = %20g",crpix1);
    i++;
    sprintf(hdrstr+i*80,"CRPIX2  = %20g",crpix2);
    
    /* Create CD matrix headers */

    i++;
    sprintf(hdrstr+i*80,"CD1_1   = %20g",cd1_1);
    i++;
    sprintf(hdrstr+i*80,"CD1_2   = %20g",cd1_2);
    i++;
    sprintf(hdrstr+i*80,"CD2_1   = %20g",cd2_1);
    i++;
    sprintf(hdrstr+i*80,"CD2_2   = %20g",cd2_2);

    /* Finally add the ZPN parameters to make this equivallent to a TAN.
       Finish with an END hdrstr+i*80 */

    i++;
    sprintf(hdrstr+i*80,"PROJP1  = 1.0");
    i++;
    sprintf(hdrstr+i*80,"PROJP3  = 0.333333");
    i++;
    sprintf(hdrstr+i*80,"PV2_0   = 0.0");
    i++;
    sprintf(hdrstr+i*80,"PV2_1   = 1.0");
    i++;
    sprintf(hdrstr+i*80,"PV2_2   = 0.0");
    i++;
    sprintf(hdrstr+i*80,"PV2_3   = 0.333333");
    i++;
    sprintf(hdrstr+i*80,"END");
    i++;

    /* Clean out any end of string marks */

    for (j = 0; j < 2880; j++) {
        if (hdrstr[j] == '\0') 
            hdrstr[j] = ' ';
    }

    /* Right, now open the WCS structure */

    retval = wcspih(hdrstr,i,0,0,&nrej,&nwcs,&wwcs);
    if (retval != 0) {
	(void)sprintf(errmsg,"NEWWCS: WCSPIH failed with status: %d\n",retval);
	freespace(hdrstr);
	wcsvfree(&nwcs,&wwcs);
	return(CIR_FATAL);
    }
    freespace(hdrstr);

    /* Get the one you want (the first one in this case) */

    *wcs = (struct wcsprm *)cir_calloc(1,sizeof(struct wcsprm));
    (*wcs)->flag = -1;
    if (wcscopy(1,wwcs,*wcs) != 0 || wcsset(*wcs) != 0) {
	(void)sprintf(errmsg,"NEWWCS: WCSCOPY/WCSSET failed\n");
	wcsvfree(&nwcs,&wwcs);
	wcsfree(*wcs);
	return(CIR_FATAL);
    }

    /* Tidy and get out of here */

    wcsvfree(&nwcs,&wwcs);    
    return(CIR_OK);
}

/*+
 *  Name:
 *      cir_crpixshift
 *
 *  Purpose:
 *      Shift the values of crpix by a certain amount
 *
 *  Description:
 *      Reads the values of CRPIX from an input header and then shifts
 *      them by a user specified amount...
 *
 *  Language:
 *      C
 *
 *  Arguments:
 *      fptr = fitsfile * (Given)
 *          The cfitsio pointer to the header with the WCS to be changed
 *      scalefac = double (Given)
 *          The amount to rescale the CRPIX vector by (divided)
 *      shifts = double [] (Given)
 *          The amount to shift CRPIX for each axis (subtracted)
 *      n = int (Given)
 *          The number of axes to be shifted
 *
 *  Returned values:
 *      Standard cirdr status value (see cirdr.h)
 *
 *  Notes:
 *      None
 *
 *  Dependencies:
 *      cfitsio
 *
 *  Authors:
 *      Jim Lewis (CASU, IoA)
 *
 *  Copyright:
 *      Copyright (C) 2002-2005 Cambridge Astronomy Survey Unit.
 *      All Rights Reserved.
 *
+*/
 
extern int cir_crpixshift(fitsfile *fptr, double scalefac, double shifts[],
                          int n) {
    int i,status;
    char comment[BUFSIZ],key[BUFSIZ];
    double val;
 
    /* Loop through and shift the values in the header... */
 
    status = 0;
    for (i = 1; i <= n; i++) {
        sprintf(key,"CRPIX%d",i);
        (void)fits_read_key(fptr,TDOUBLE,key,&val,comment,&status);
        if (status != 0)
            return(CIR_FATAL);
        val = (val - shifts[i-1])/scalefac - 1.0;
        (void)fits_update_key(fptr,TDOUBLE,key,&val,comment,&status);
    }
    return(CIR_OK);
}
 
/*+
 *  Name:
 *      cir_rescalecd
 *
 *  Purpose:
 *      Rescale the CD matrix by a certain amount
 *
 *  Description:
 *      Reads the values of CD matrix from an input header and then multiplies
 *      them by a user specified factor. Does the same for the CRPIX vector
 *
 *  Language:
 *      C
 *
 *  Arguments:
 *      fptr = fitsfile * (Given)
 *          The cfitsio pointer to the header with the WCS to be changed
 *      scalefac = double (Given)
 *          The amount to rescale the CD matrix by (multiplied)
 *
 *  Returned values:
 *      Standard cirdr status value (see cirdr.h)
 *
 *  Notes:
 *      Only goes up to 2d. Couldn't be bothered to to the general case.
 *
 *  Dependencies:
 *      cfitsio
 *
 *  Authors:
 *      Jim Lewis (CASU, IoA)
 *
 *  Copyright:
 *      Copyright (C) 2003-2006 Cambridge Astronomy Survey Unit.
 *      All Rights Reserved.
 *
+*/
 
extern int cir_rescalecd(fitsfile *fptr, double scalefac) {
    int i,j,status;
    char comment[BUFSIZ],key[BUFSIZ];
    double val;
 
    /* Loop through and shift the values in the header... */
 
    status = 0;
    for (i = 1; i <= 2; i++) {
        for (j = 1; j <= 2; j++) {
            sprintf(key,"CD%d_%d",i,j);
            (void)fits_read_key(fptr,TDOUBLE,key,&val,comment,&status);
            if (status != 0)
                return(CIR_FATAL);
            val *= scalefac;
            (void)fits_update_key(fptr,TDOUBLE,key,&val,comment,&status);
        }
    }
    return(CIR_OK);
}

static void fk425pv (double *ra, double *dec, double *rapm, double *decpm,
		     double *parallax, double *rv) {

/* double	*ra;	   Right ascension in degrees (J2000 in, B1950 out) */
/* double	*dec;	   Declination in degrees (J2000 in, B1950 out) */
/* double	*rapm;	   Proper motion in right ascension */
/* double	*decpm;	   Proper motion in declination
			 * In:  ra/dec degrees/Julian year
			 * Out: ra/dec degrees/tropical year */
/* double       *parallax; Parallax (arcsec) */
/* double       *rv;	   Rradial velocity (km/s, +ve = moving away) */

/* This routine converts stars from the old, Bessel-Newcomb, FK4
   system to the new, IAU 1976, FK5, Fricke system, using Yallop's
   implementation (see ref 2) of a matrix method due to Standish
   (see ref 3).  The numerical values of ref 2 are used canonically.

   Notes:

      1)  The proper motions in ra are dra/dt rather than
 	   cos(dec)*dra/dt, and are per year rather than per century.

      2)  Conversion from besselian epoch 1950.0 to Julian epoch
 	   2000.0 only is provided for.  Conversions involving other
 	   epochs will require use of the appropriate precession,
 	   proper motion, and e-terms routines before and/or
 	   after fk425 is called.

      3)  In the FK4 catalogue the proper motions of stars within
 	   10 degrees of the poles do not embody the differential
 	   e-term effect and should, strictly speaking, be handled
 	   in a different manner from stars outside these regions.
 	   However, given the general lack of homogeneity of the star
 	   data available for routine astrometry, the difficulties of
 	   handling positions that may have been determined from
 	   astrometric fields spanning the polar and non-polar regions,
 	   the likelihood that the differential e-terms effect was not
 	   taken into account when allowing for proper motion in past
 	   astrometry, and the undesirability of a discontinuity in
 	   the algorithm, the decision has been made in this routine to
 	   include the effect of differential e-terms on the proper
 	   motions for all stars, whether polar or not.  At epoch 2000,
 	   and measuring on the sky rather than in terms of dra, the
 	   errors resulting from this simplification are less than
 	   1 milliarcsecond in position and 1 milliarcsecond per
 	   century in proper motion.

   References:

      1  "Mean and apparent place computations in the new IAU System.
          I. The transformation of astrometric catalog systems to the
 	  equinox J2000.0." Smith, C.A.; Kaplan, G.H.; Hughes, J.A.;
	  Seidelmann, P.K.; Yallop, B.D.; Hohenkerk, C.Y.
 	  Astronomical Journal vol. 97, Jan. 1989, p. 265-273.

      2  "Mean and apparent place computations in the new IAU System.
	  II. Transformation of mean star places from FK4 B1950.0 to
 	  FK5 J2000.0 using matrices in 6-space."  Yallop, B.D.;
	  Hohenkerk, C.Y.; Smith, C.A.; Kaplan, G.H.; Hughes, J.A.;
	  Seidelmann, P.K.; Astronomical Journal vol. 97, Jan. 1989,
	  p. 274-279.

      3  "Conversion of positions and proper motions from B1950.0 to the
	  IAU system at J2000.0", Standish, E.M.  Astronomy and
	  Astrophysics, vol. 115, no. 1, Nov. 1982, p. 20-22.

   P.T.Wallace   Starlink   20 December 1993
   Doug Mink     Smithsonian Astrophysical Observatory  7 June 1995 */

    double r1950,d1950;		/* B1950.0 ra,dec (rad) */
    double r2000,d2000;		/* J2000.0 ra,dec (rad) */

    /* Miscellaneous */

    double ur,ud,sr,cr,sd,cd,w,wd;
    double x,y,z,xd,yd,zd;
    double rxyz, rxysq, rxy, rxyzsq, spxy, spxyz;
    int	i,j;

    double r0[3],rd0[3];	/* star position and velocity vectors */
    double v1[6],v2[6];		/* combined position and velocity vectors */

    /* Constants */

    double zero = (double) 0.0;
    double vf = 21.095;	/* Km per sec to AU per tropical century */
			/* = 86400 * 36524.2198782 / 149597870 */
    double a[3] = {-1.62557e-6, -0.31919e-6, -0.13843e-6};
    double ad[3] = {1.245e-3,  -1.580e-3,  -0.659e-3};
    double d2pi = 6.283185307179586476925287;	/* two PI */
    double tiny = 1.e-30; /* small number to avoid arithmetic problems */

    /* Convert B1950.0 FK4 star data to J2000.0 FK5 */

    double em[6][6] = {
    {	 0.9999256782,		/* em[0][0] */
	-0.0111820611,		/* em[0][1] */
	-0.0048579477,		/* em[0][2] */
	 0.00000242395018,	/* em[0][3] */
	-0.00000002710663,	/* em[0][4] */
	-0.00000001177656 },	/* em[0][5] */
 
    {	 0.0111820610,		/* em[1][0] */
	 0.9999374784,		/* em[1][1] */
	-0.0000271765,		/* em[1][2] */
	 0.00000002710663,	/* em[1][3] */
	 0.00000242397878,	/* em[1][4] */
	-0.00000000006587 },	/* em[1][5] */
 
    {	 0.0048579479,		/* em[2][0] */
	-0.0000271474,		/* em[2][1] */
	 0.9999881997,		/* em[2][2] */
	 0.00000001177656,	/* em[2][3] */
	-0.00000000006582,	/* em[2][4] */
	 0.00000242410173 },	/* em[2][5] */
 
    {	-0.000551,		/* em[3][0] */
	-0.238565,		/* em[3][1] */
	 0.435739,		/* em[3][2] */
	 0.99994704,		/* em[3][3] */
	-0.01118251,		/* em[3][4] */
	-0.00485767 },		/* em[3][5] */
 
    {	 0.238514,		/* em[4][0] */
	-0.002667,		/* em[4][1] */
	-0.008541,		/* em[4][2] */
	 0.01118251,		/* em[4][3] */
	 0.99995883,		/* em[4][4] */
	-0.00002718 },		/* em[4][5] */
 
    {	-0.435623,		/* em[5][0] */
	 0.012254,		/* em[5][1] */
	 0.002117,		/* em[5][2] */
	 0.00485767,		/* em[5][3] */
	-0.00002714,		/* em[5][4] */
	 1.00000956 }		/* em[5][5] */
    };



    /* Convert B1950 RA and Dec from degrees to radians */

    r1950 = DEGRAD*(*ra);
    d1950 = DEGRAD*(*dec);

    /* Convert B1950 RA and Dec proper motion from degrees/year to arcsec/tc */

    ur = *rapm  * 360000.0;
    ud = *decpm * 360000.0;

    /* Convert direction to Cartesian */

    sr = sin (r1950);
    cr = cos (r1950);
    sd = sin (d1950);
    cd = cos (d1950);
    r0[0] = cr * cd;
    r0[1] = sr * cd;
    r0[2] = sd;

    /* Convert motion to Cartesian */

    w = vf * *rv * *parallax;
    if (ur != zero || ud != zero || (*rv != zero && *parallax != zero)) {
	rd0[0] = (-sr * cd * ur) - (cr * sd * ud) + (w * r0[0]);
	rd0[1] =  (cr * cd * ur) - (sr * sd * ud) + (w * r0[1]);
	rd0[2] = 	                (cd * ud) + (w * r0[2]);
    } else {
	rd0[0] = zero;
	rd0[1] = zero;
	rd0[2] = zero;
    }

    /* Remove e-terms from position and express as position+velocity 6-vector */
    w = (r0[0] * a[0]) + (r0[1] * a[1]) + (r0[2] * a[2]);
    for (i = 0; i < 3; i++)
	v1[i] = r0[i] - a[i] + (w * r0[i]);

    /* Remove e-terms from proper motion and express as 6-vector */

    wd = (r0[0] * ad[0]) + (r0[1] * ad[1]) + (r0[2] * ad[2]);
    for (i = 0; i < 3; i++)
	v1[i+3] = rd0[i] - ad[i] + (wd * r0[i]);

    /* Alternately: Put proper motion in 6-vector without adding e-terms

    for (i = 0; i < 3; i++)
	v1[i+3] = rd0[i]; */

    /* Convert position + velocity vector to FK5 system */

    for (i = 0; i < 6; i++) {
	w = zero;
	for (j = 0; j < 6; j++)
	    w += em[i][j] * v1[j];
	v2[i] = w;
    }

    /* Vector components */

    x = v2[0];
    y = v2[1];
    z = v2[2];
    xd = v2[3];
    yd = v2[4];
    zd = v2[5];

    /* Magnitude of position vector */

    rxysq = x*x + y*y;
    rxy = sqrt (rxysq);
    rxyzsq = rxysq + z*z;
    rxyz = sqrt (rxyzsq);

    spxy = (x * xd) + (y * yd);
    spxyz = spxy + (z * zd);

    /* Convert back to spherical coordinates */

    if (x == zero && y == zero)
	r2000 = zero;
    else {
	r2000 = atan2 (y,x);
	if (r2000 < zero)
	    r2000 = r2000 + d2pi;
    }
    d2000 = atan2 (z,rxy);

    if (rxy > tiny) {
	ur = ((x * yd) - (y * xd)) / rxysq;
	ud = ((zd * rxysq) - (z * spxy)) / (rxyzsq * rxy);
    }

    if (*parallax > tiny) {
	*rv = spxyz / (*parallax * rxyz * vf);
	*parallax = *parallax / rxyz;
    }

    /* Return results */

    *ra = r2000/DEGRAD;
    *dec = d2000/DEGRAD;
    *rapm  = ur / 360000.0;
    *decpm = ud / 360000.0;

    return;
}

/*
 
$Log: cir_wcssubs.c,v $
Revision 1.14  2005/03/08 20:05:46  jim
Removed extraneous declaration to keep compiler happy

Revision 1.13  2005/02/25 10:17:51  jim
Fixed bug in routine that creates a new WCS

Revision 1.12  2004/09/09 10:51:57  jim
Changed malloc to calloc in cir_wcsopen.

Revision 1.11  2004/09/07 14:18:55  jim
Tidied up so that strict compilation doesn't throw up loads of warnings

Revision 1.10  2004/08/19 11:34:25  jim
Added new cir_memory routines for memeory allocation

Revision 1.9  2004/08/04 10:52:58  jim
Fixed bug in cir_wcsopen and cir_wcsclose

Revision 1.8  2004/08/02 11:49:33  jim
New version using wcslib

Revision 1.7  2004/05/05 10:52:45  jim
Modified routine cir_crpixshift to scale the value of crpix as well as shift
it
 
Revision 1.6  2004/01/23 08:48:25  jim
Fixed bug in cir_newwcs where PV2_0 had the wrong value
 
Revision 1.5  2003/07/10 13:45:40  jim
Added PV matrix info into cir_newwcs
 
Revision 1.4  2003/04/16 20:08:54  jim
Added cir_rescalecd and cir_crpixshift
 
Revision 1.3  2002/12/16 10:25:06  jim
Added prologues
 
Revision 1.2  2002/08/26 07:31:51  jim
Added cir_newwcs
 
Revision 1.1.1.1  2002/06/21 09:48:57  jim
Initial import into CVS
 
 
*/

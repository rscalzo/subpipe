/*
 
$Id: cir_wcssubs.h,v 1.6 2004/09/07 14:20:47 jim Exp $
 
*/

#define DEGRAD  0.017453292519943

#include <wcs.h>
#include <wcshdr.h>
#include <fitsio.h>
#include <longnam.h>

extern int cir_wcsopen(char *, struct wcsprm **, char *);
extern void cir_wcsclose(struct wcsprm *); 
extern void cir_xytoradec(struct wcsprm *, double, double, double *, double *);
extern void cir_radectoxy(struct wcsprm *, double, double, double *, double *);
extern void cir_radectoxieta(struct wcsprm *, double, double, double *,
			     double *);
extern void cir_fk425(double *, double *);
extern int cir_newwcs(struct wcsprm **, double, double, double, double, double,
		      double, double, double, int, int, char *);
extern int cir_crpixshift(fitsfile *, double, double [], int);
extern int cir_rescalecd(fitsfile *, double);

/*
 
$Log: cir_wcssubs.h,v $
Revision 1.6  2004/09/07 14:20:47  jim
Tidied up to avoid warning messages

Revision 1.5  2004/08/02 11:48:35  jim
New version with wcslib

Revision 1.4  2004/05/05 10:56:37  jim
Modifed call to cir_crpixshift
 
Revision 1.3  2003/04/16 20:10:21  jim
Added cir_rescalecd and cir_crpixshift
 
Revision 1.2  2002/08/26 07:28:03  jim
Added prototype of cir_newwcs
 
Revision 1.1.1.1  2002/06/21 09:48:57  jim
Initial import into CVS
 
 
*/

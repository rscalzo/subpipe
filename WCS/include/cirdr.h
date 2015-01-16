/*

$Id: cirdr.h,v 1.57 2005/05/31 07:37:56 jim Exp $

*/
#include <stdlib.h>
#include <stdio.h>
#include <fitsio.h>
#include <longnam.h>
#include <wcs.h>
#include <cir_memory.h>

/* Error status codes */

#define CIR_OK    0
#define CIR_WARN  1
#define CIR_FATAL 2

/* Codes for type of statistic */

#define MEANCALC    1
#define MEDIANCALC  2

/* Codes for type of weighting to be done in combinations */

#define NOBPM      0
#define USEBPM     1
#define USEWEIGHT  2

/* Codes for supported rejection algorithms */

#define REJNONE      0
#define REJPOISSON   1
#define REJAVSIGCLIP 2

/* Codes for sky subtraction algorithms */

#define NOSKYSUB        "NONE"
#define OFFSETSKY       "OFFSET_SKY"
#define TIMED_OFFSETSKY "OFFSET_SKY_TIMED"
#define WINDOWSKY       "SLIDING_WINDOW"
#define WINDOWSKY_RESTRICT "SLIDING_WINDOW_RESTRICTED"
#define TILESKY         "TILE_SKY"

/* Codes for reset correction algorithms */

#define NORESET         "NONE"
#define PLANE2M1        "PLANE2M1"
#define EXTN2M1         "EXTN2M1"

/* Codes for debiassing */

#define NODEBIAS         "NONE"
#define FRAMEBIAS        "FRAME"
#define CONSTANTBIAS     "CONSTANT"
#define FRAME_OSBIAS     "FRAME_OS"

/* Codes for 'tarting up' frames */

#define RANOM_FILT   "RESET_ANOM_FILTER"
#define RANOM_SURF   "RESET_ANOM_SURFACE"
#define TWOPT_X      "TWOPOINT_X"
#define TWOPT_Y      "TWOPOINT_Y"
#define DESTR_X      "DESTRIPE_X"
#define DESTR_Y      "DESTRIPE_Y"

/* Codes for defringing algorithms */

#define FRINGE_AUTO   1
#define FRINGE_MANUAL 2

/* Types of masks */

#define MASK_BPM 1
#define MASK_CPM 2
#define MASK_OPM 3

/* Codes for calibration frame types */

#define CALFLAT    "FLAT"
#define CALFRINGE  "FRINGE"
#define CALLIN     "LIN"
#define CALBIAS    "BIAS"
#define CALSKY     "SKY"
#define CALCPM     "CPM"
#define CALDARK    "DARK"
#define CALBPM     "BPM"

/* Useful macro */

#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) < (B) ? (A) : (B))

/* Macros that deal with tidying up */

#define freespace(_p) if (_p != NULL) {free(_p); _p = NULL;}
#define closefits(_p) status = 0; if (_p != NULL) {(void)fits_close_file(_p,&status); _p = NULL;}
#define deletefits(_p) status = 0; if (_p != NULL) {(void)fits_delete_file(_p,&status); _p = NULL;}
#define closefile(_p) if (_p != NULL) {fclose(_p); _p = NULL;}
#define freewcs(_p) if (_p != NULL) {cir_wcsclose(_p); _p = NULL;}

/* Oft used constants */

#define M_PI            3.14159265358979323846  /* pi */
#define M_SQRT2         1.41421356237309504880  /* sqrt(2) */
#define M_LN2           0.69314718055994530942  /* log_e 2 */

/* All exportable CIRDR routines */

extern int cir_apm(char *, char *, char *, char *, int, int, float, 
                   int, int, char *);
extern int cir_bpmcreate(char **, char *, float, int, int, char *);
extern int cir_catcoord(char *, char *, char *);
extern int cir_ccdproc(char *, char *, char *, char *, float, char *, char *,
                       float, float, char *, char *);
extern int cir_chanback(char *, char *, char *);
extern int cir_chanstats(char **, int, char *, char *, char *);
extern int cir_compress_fits(char *, char *, char *);
extern int cir_compress_fits_mef_f2i(char *, double, char *, char *);
extern int cir_coverage(char *, int, float *, float *, float *, float *, 
			char *);
extern int cir_crosstalk(char **, int, char *, char *, char *);
extern int cir_defringe(char *, char *, char *, char *, float, float, float,
			int, float, int, int, int, float, float, int, char *);
extern int cir_dmean(double *, unsigned char *, int, double *, char *);
extern int cir_dmed(double *, unsigned char *, int, double *, char *);
extern int cir_fback(char *, int, int, char *);
extern void cir_fk425(double *, double *);
extern int cir_getstds(char *, char *, char *, char *, char *, float, int, int,
		       char *);
extern int cir_imadd(char *, char *, char *);
extern int cir_imaddk(char *, float, char *);
extern int cir_imavsig(char *, char *, float *, float *, char *);
extern int cir_imcombine(char *, char *, char *, char *, int, int, int, int, 
                         int, int, int, int, int, float, float, char *, 
                         char *, char *, char *);
extern int cir_imcombine_lite(char **, int, int, int, int, float, char *, 
			      char *, char *);
extern int cir_imcore(char *, char *, int, float, int, float, int, char *, 
		      char *, int, int, char *);
extern int cir_imdiv(char *, char *, char *);
extern int cir_imdivk(char *, float, char *);
extern int cir_immult(char *, char *, char *);
extern int cir_immultk(char *, float, char *);
extern int cir_imqmedsig(char *, char *, float, int, float, float, float *, 
			 float *, char *);
extern int cir_imrepl(char *, float, float, char *);
extern int cir_imslice(char *, char *, int, char *);
extern int cir_imsqrt(char *, char *);
extern int cir_imsubt(char *, char *, char *);
extern int cir_imsubtk(char *, float, char *);
extern int cir_imstack(char *, char *, char *, int, int, int, int, int, int, 
                       int, int, float, float, char *, char *, char *, char *);
extern int cir_imstat(char *, char *, int, float *, char *);
extern int cir_interleaf(char *, char *, char *, char *, int, int, int, 
			 char *);
extern int cir_lincoefs_postreset(char **, int, char *, int, char *, float,
				  float, char *, char *);
extern int cir_lincoefs_postreset_allpix(char **, int, char *, int, char *, 
					 float, float, char *, char *);
extern int cir_lincorrect_postreset(char *, char *, char *, char *, float,
				    float, char *);
extern int cir_lincorrect_postreset_allpix(char *, char *, char *, char *, 
					   float, float, char *);
extern int cir_makemask(char *, char *, float, int, char *);
extern int cir_matchstds(char *, char *, float, int, int, char *, int *, 
                         char *);
extern int cir_matchxy(char *, char *, float, int, int, float *, float *,
		       int *, char *);
extern int cir_mean(float *, unsigned char *, int, float *, char *);
extern int cir_meansig(float *, unsigned char *, int, float *, float *, 
		       char *);
extern int cir_meansigcut(float *, unsigned char *, int, float, float, float *,
			  float *, char *);
extern int cir_med(float *, unsigned char *, int, float *, char *);
extern int cir_med_zs(float *, unsigned char *, float *, float *, int, float *,
		      char *);
extern int cir_mkconf(char *, char *, float, int, int, char *);
extern int cir_mkconf_2(char **, char *, float, int, int, char *);
extern int cir_mkconf_3(char *, char *, char *, char *);
extern int cir_mscale_subt(char *, char *, float, char *, char *);
extern int cir_offon(char *, char *, char **, char **, int, float *, char *);
extern int cir_persist(char *, char *, long, long, char *);
extern int cir_pixstats(char **, int, char *, char *);
extern int cir_platesol(char *, char *, int, int, int, char *);
extern int cir_pmcombine(char *, char *, char *, char *);
extern int cir_prep_fringe_file(char *, int, int, int, int, char *, char *, 
				char *);
extern void cir_prochist(char *);
extern int cir_qblkmed(char *, char *, int, float, char *);
extern int cir_quadback(char *, char *);
extern int cir_qmedsig(float *, unsigned char *, int, float, int, float, float,
		       float *, float *, char *);
extern int cir_radec2xy(char *, double, double, double *, double *, char *);
extern int cir_stage1(char *, char *, char *, char *, char *, float, char *, 
		      char *);
extern void cir_stamp(char *, char *);
extern int cir_sumbpm(unsigned char *, int, int *, char *);
extern int cir_tartup(char *, char *);
extern int cir_tabledump(char *, char *, char **, int, char *);
extern int cir_testbias(char *, char *, char *, int, float, int *, char *);
extern int cir_uncompress(char *, char *, char *);
extern int cir_update_hdr(char *, char *, char *, char *, char *, char *);
extern int cir_wcsoffset(char *, char *, int, char *);
extern int cir_wmean(float *, float *, int, float *, float *, char *);
extern int cir_xy2radec(char *, double, double, double *, double *, char *);
extern int cir_xyarrays_match(float *, float *, int, float *, float *, int,
			      float, int *, int *, char *);

/* WFCAM specific routines */

extern int cir_wfcam_ustep_pt1(char **, int, char *, char *, char *, char *,
			       char *, char *, char *, float, int, int, float,
			       int, float, int, int, char *);
extern int cir_wfcam_mskysubt(char **, int, char *, float, char *);


/*

$Log: cirdr.h,v $
Revision 1.57  2005/05/31 07:37:56  jim
added cir_xyarrays_match

Revision 1.56  2005/05/25 09:06:38  jim
Added exptime to argument list for cir_imcombine_lite

Revision 1.55  2005/05/23 08:25:45  jim
Modified argument list for cir_qblkmed

Revision 1.54  2005/05/13 09:27:00  jim
Added cir_imqmedsig

Revision 1.53  2005/04/19 12:37:06  jim
added cir_uncompress

Revision 1.52  2005/03/21 16:14:57  jim
Added cir_lincoef_postreset_allpix

Revision 1.51  2005/03/14 21:12:48  jim
added cir_wfcam_mskysubt

Revision 1.50  2005/03/10 22:48:53  jim
Added cir_wfcam_ustep_pt1

Revision 1.49  2005/03/08 20:06:45  jim
Added cir_mkconf_3

Revision 1.48  2005/03/07 20:38:06  jim
Added definition of CALBPM and cir_bpmcreate

Revision 1.47  2005/03/03 22:15:15  jim
Added cir_imcombine_lite and cir_compress_mef_f2i

Revision 1.46  2005/02/25 10:13:14  jim
Added cir_offon

Revision 1.45  2005/02/08 10:32:38  jim
Added cir_chanback

Revision 1.44  2005/02/04 10:30:08  jim
Added cir_lincoefs_postreset_allpix

Revision 1.43  2004/11/26 19:58:27  jim
Modified parameters to cir_qblkmed

Revision 1.42  2004/11/22 21:01:00  jim
Added cir_quadback and cir_qblkmed

Revision 1.41  2004/10/19 15:34:03  jim
Added cir_mkconf_2

Revision 1.40  2004/09/21 13:58:15  jim
Added entry for linearity routines

Revision 1.39  2004/09/09 10:50:31  jim
Added cir_imavsig

Revision 1.38  2004/09/07 14:20:47  jim
Tidied up to avoid warning messages

Revision 1.37  2004/09/03 10:49:33  jim
added cir_imslice

Revision 1.36  2004/09/02 12:56:54  jim
Added new arithmetic routines

Revision 1.35  2004/08/19 11:40:20  jim
Added reference to cir_memory.h

Revision 1.34  2004/08/02 11:48:14  jim
New definition of freewcs

Revision 1.33  2004/07/26 11:09:06  jim
Added cir_chanstats

Revision 1.32  2004/07/23 06:39:45  jim
Modified entry for cir_crosstalk

Revision 1.31  2004/06/09 14:40:40  jim
Added cir_pixstats

Revision 1.30  2004/05/18 13:30:18  jim
Added cir_compress_fits and removed cir_compress_simple_fits

Revision 1.29  2004/05/05 10:56:19  jim
Defined macro 'deletefits'

Revision 1.28  2004/04/05 11:27:28  jim
Modified for new version of cir_imcore

Revision 1.27  2004/03/30 13:36:11  jim
Added cir_crosstalk

Revision 1.26  2004/03/09 12:40:25  jim
Added cir_compress_simple_fits

Revision 1.25  2004/01/30 14:25:51  jim
Added cache parameter to cir_getstds

Revision 1.24  2003/11/27 10:26:44  jim
Added new sky subtraction algorithm WINDOWSKY_RESTRICT

Revision 1.23  2003/11/13 11:48:24  jim
Added cir_imaddk routine

Revision 1.22  2003/11/06 16:34:17  jim
Added CALDARK constant

Revision 1.21  2003/10/05 20:56:01  jim
Added cir_tabledump

Revision 1.20  2003/10/04 00:28:30  jim
Added cir_update_hdr

Revision 1.19  2003/09/30 20:28:55  jim
Added cir_radec2xy

Revision 1.18  2003/09/16 11:37:15  jim
Added cir_xy2radec

Revision 1.17  2003/09/11 11:58:57  jim
Added in revised statistical routines

Revision 1.16  2003/09/01 10:34:27  jim
Modified entry for cir_getstds

Revision 1.15  2003/07/05 09:05:43  jim
Added nbsize parameter to cir_imcore.c

Revision 1.14  2003/05/06 07:58:24  jim
Added cir_persist

Revision 1.11  2003/03/14 13:01:06  jim
Added calibration frame types

Revision 1.10  2003/02/03 09:10:32  jim
Added cir_prochist

Revision 1.9  2003/02/03 09:08:51  jim
Added cir_stamp

Revision 1.8  2002/12/12 12:09:29  jim
Added to calling sequence for cir_coverage

Revision 1.7  2002/11/11 09:31:00  jim
Added TILESKY option for sky subtraction

Revision 1.6  2002/08/02 08:29:10  jim
Added entry for cir_coverage

Revision 1.5  2002/07/19 10:23:56  jim
Added entry for cir_testbias

Revision 1.4  2002/07/08 12:22:43  jim
Removed masktype parameter from argument list to cir_imstat

Revision 1.3  2002/07/04 12:36:06  jim
Changed prototype for cir_imcombine to add maskloc parameter

Revision 1.2  2002/07/04 12:13:12  jim
Added MASK_OPM global

Revision 1.1.1.1  2002/06/21 09:48:57  jim
Initial import into CVS


*/

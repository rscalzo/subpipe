/*

$Id: cir_stamp.c,v 1.4 2004/08/19 11:34:25 jim Exp $

*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include <cirdr.h>

#define SZ_HIST 70

/*+
 *  Name:
 *      cir_stamp
 *
 *  Purpose:
 *      Write a time stamp to a FITS image's history
 *
 *  Description:
 *      Two history cards are added to the header of a given FITS image.
 *      The first is the current date and time.  The second is a given
 *      timestamp for the calling routine.  This should be something like
 *      the CVS ID string.
 *
 *  Arguments:
 *      fname = char * (Given)
 *          The full name for the input FITS image
 *      stamp = char * (Given)
 *          Some sort of timestamp to be written to the history record.
 *          This could be a CVS ID or something of the sort.
 *
 *  Returned values:
 *      None
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

extern void cir_stamp (char *fname, char *stamp) {
    fitsfile *iptr;
    struct tm *tp;
    time_t time_p;
    int status;
    char datestr[BUFSIZ];

    /* Get time info */

    time(&time_p);
    tp = localtime(&time_p);
    strftime(datestr,BUFSIZ,"%Y%m%d %H:%M:%S",tp);

    /* Open the file and write the info to the header */

    status = 0;
    (void)fits_open_file(&iptr,fname,READWRITE,&status);
    if (status != 0) 
        return;
    (void)fits_write_history(iptr,datestr,&status);
    sprintf(datestr,"   %s",stamp);
    (void)fits_write_history(iptr,datestr,&status);
    (void)fits_close_file(iptr,&status);
}

/*+
 *  Name:
 *      cir_prochist
 *
 *  Purpose:
 *      Reorganise all the history records in a single FITS header 
 *
 *  Description:
 *      The history cards from a single FITS header are copied into memory
 *      and then deleted.  Then they are recopied onto the bottom of the
 *      header.
 *
 *  Arguments:
 *      fname = char * (Given)
 *          The full name for the input FITS image
 *
 *  Returned values:
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

extern void cir_prochist(char *fname) {
    fitsfile *iptr;
    int nkeys,status,i,nhist;
    char *histstr,*hstr,val[SZ_HIST],com[SZ_HIST],key[SZ_HIST];

    /* Open the FITS file */

    status = 0;
    (void)fits_open_file(&iptr,fname,READWRITE,&status);
    if (status != 0)
        return;

    /* Get the size of the header block and allocate enough memory
       for that number of history cards. */

    (void)fits_get_hdrspace(iptr,&nkeys,NULL,&status);
    histstr = cir_malloc(nkeys*SZ_HIST*sizeof(char *));

    /* Run through the header and mark all of the history cards. Then
       delete them.*/

    nhist = 0;
    i = 1;
    while (1) {
        (void)fits_read_keyn(iptr,i,key,val,com,&status);
        if (status != 0) {
            status = 0;
            break;
        }
        if (!strcmp(key,"HISTORY")) {
            nhist++;
            hstr = histstr + (nhist - 1)*SZ_HIST;
            (void)strcpy(hstr,com);
            (void)fits_delete_record(iptr,i,&status);
        } else {
            i++;
        }
    }

    /* Now write the history cards all to the end of the file */

    for (i = 0; i < nhist; i++)
        (void)fits_write_history(iptr,histstr+i*SZ_HIST,&status);

    /* Tidy and exit */

    free(histstr);
    (void)fits_close_file(iptr,&status);
}

/*

$Log: cir_stamp.c,v $
Revision 1.4  2004/08/19 11:34:25  jim
Added new cir_memory routines for memeory allocation

Revision 1.3  2003/11/07 13:05:43  jim
Removed some debugging info

Revision 1.2  2003/11/07 12:58:28  jim
Modified cir_prochist so that the declaration for the FITS keyword is now
much bigger than 8 characters.  This is to take ESO keywords into account

Revision 1.1  2003/02/03 09:07:29  jim
First entry into CVS


*/

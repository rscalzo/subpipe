

import numpy as np
from bisect import bisect_left


def close_match_radec(ra1,dec1,ra2,dec2,ep,allow=1,silent=0,box=0):

#+
# NAME:
# close_match_radec
#
# PURPOSE:
# this will find close matches between 2 sets of points (ra1,dec1)
# and (ra2,dec2) (note ra1,ra2,dec1,dec2 are all arrays) in ra dec space.
# 
# CALLING SEQUENCE:
# close_match_radec(ra1,dec1,ra2,dec2,m1,m2,ep,allow,miss1,silent=1,box=1)
#
# INPUTS:
# ra1,dec1: the ra dec of the first set
# ra2,dec2: the ra dec of the second set
# ep:  this is the error which defines a close match. A pair is considered
# a match if their spherical separation is <= ep.
# allow: how many matches in the (ra2,dec2) space will you allow for
# each (ra1,dec1)
#
# OUTPUTS:
# m1,m2: the indices of the matches in each space. That is  
# (ra1(m1),dec1(m1)) matches (ra2(m2),dec2(m2))
# miss1: this gives the index of the things in x1 NOT found to match (optional)
#
# OPTIONAL KEYWORD PARAMETERS:
# silent: no printed output
# box: old square box matching instead of new circular matching
#
# NOTES:
# It sorts the ra2 list so that it can do a binary search for each ra1.
# It then carves out the sublist of ra2 where it is off by less than ep.
# It then only searches that sublist for points where dec1 and dec2 differ
# by less than ep. 
# PROCEDURES CALLED:
# binary_search, rem_dup
# REVISION HISTORY:
# written by David Johnston -University of Michigan June 97
#
#   Tim McKay August 97
# 	bug fixed, matcharr extended to "long" to accomodate ROTSE I images
#   Tim McKay 6/23/98
#	A<ered to be an ra dec match, with appropriate scaling of ra range...
#   Tim McKay 7/8/99
#	A<ered to order matches by distance, rather than just ra dec distance
#   Erin Scott Sheldon 08-Mar-2001
#       Made code human readable. Fixed bug where miss1 was not returned
#       when no matches found.
#   E.S.S.  Fixed bug where some matches could be returned outside of the
#           requested distance cut.  10-June-2004
#   FY IDL to Python
#-
 
# first sort out the allowed errors in ra and dec.....
 
    epdec=ep
 
    n1=len(ra1)
    n2=len(ra2)

    matcharr=np.zeros((n1,allow),dtype='int32')	#the main book-keeping device for 
    matcharr[:]=-1		#matches -initialized to -1
    ind=np.arange(n2)
    sor=np.argsort(ra2)  			#sort ra2
    ra2sort=ra2[sor]
    dec2sort=dec2[sor]
    ind=ind[sor]			#keep track of index
    runi=0
    endra2=ra2sort[n2-1]
    for i in range(n1):      #the main top level loop over ra1
     
        epra=ep/np.cos(dec1[i]*0.01745329)
        ra1minus = ra1[i]-epra     #sets the range of good ones
        ra1plus  = ra1[i]+epra
        in1=find_lt(ra2sort,ra1minus) 
                                #searched for the first good one
        if in1 == -1: 
            if ra1minus < endra2: in1=0 
                                #in case ra1minus smaller than all ra2 
                                #but still may be some matches

        if in1 !=-1: 

        ## OK, we have the first match in the sorted list. Look at 
        ## objects that come after until we find no more matches

            in2 = in1
            jj = in2+1
        while jj < n2 :
            if ra2sort[in2+1] < ra1plus :
                in2=in2+1
                jj=jj+1
            else: jj=n2
         
        if n2 ==1 : in2 =0

        ## while loop carved out sublist to check

        if in1 <= in2 :  
            dec2check=dec2sort[in1:in2+1] #the sublist to check
            ra2check=ra2sort[in1:in2+1]

            ## look in box
            decoffby=abs(dec2check-dec1[i])
            raoffby=abs(ra2check-ra1[i])
            good=[goodind for goodind in range(len(decoffby)) if (decoffby[goodind] < epdec and raoffby[goodind] <epra)]
            good=np.array(good)+in1
            ngood=len(good)
                        #the selection made here
            if ngood != 0 :  

                ######################################################
                ## Above was only a square cut
                ## Get those that are actually within the
                ## requested radius. E.S.S. 10-June-2004
                ######################################################
                if (not box) : 
                    offby=sphdist(ra1[i],dec1[i],ra2sort[good],dec2sort[good])
                    good_offby = np.where(offby <= ep)[0]
                    ngood=len(good_offby)
                else : 
                    good_offby = np.arange(ngood)
                    offby = raoffby
                    ## ngood is the same
                 
                if ngood != 0 :

                    good = good[good_offby]
                    offby = offby[good_offby]
                    if ngood > allow :  

                        ## Sort by closeness
                        good=good[np.argsort(offby)] 

                        ## not more than you are allowed by 'allow' 
                        ngood=allow 
                        good=good[0:allow]                         

                    matcharr[i,0:ngood]=good
                                #finally the matches are entered in 
                    runi=runi+ngood #a running total    
                     

    if not silent: print 'total put in bytarr',runi
    matches=np.where(matcharr != -1)
    this=len(matches)
    if this == 0:
        if not silent: print 'no matches found'
        miss1=np.arange(n1)
        m1=-1
        m2=-1
        return m1,m2,miss1

    m1=matches[0]              #a neat trick to extract them correctly 
    m2=matcharr[matches[0],matches[1]]           #from the matcharr matrix
    if not silent: print len(m1),' matches'
    m2=ind[m2]                     #remember, must unsort
    dif=list(set(m1))
    # RS 2012/02/10:  Added this to prevent occasional failure; the code below
    # seems to assume that m1 will be sorted, when in fact sets have no order
    # and so we have to sort it ourselves.
    dif.sort()

    if not silent: print len(dif),' different matches'
    
    if len(dif) < n1 : 
        miss1=range(n1)
        for difind in dif.__reversed__():
            del miss1[difind]
        if not silent: print len(miss1),'  misses'
    else : 
        miss1=[]  
        if not silent: print 'no misses'
    return m1,m2,miss1


def find_lt(a, x):
    'Find rightmost value less than x'
    i = bisect_left(a, x)
    return i-1
    

def sphdist (ra1, dec1, ra2, dec2):
    """measures the spherical distance in degrees
        The input has to be in degrees too
        """
    dec1_r = np.deg2rad(dec1)
    dec2_r = np.deg2rad(dec2)
    dist=2*np.rad2deg(np.arcsin(np.sqrt((np.sin((dec1_r-dec2_r)/2))**2+np.cos(dec1_r)*np.cos(dec2_r) *(np.sin((np.deg2rad(ra1-ra2))/2))**2)))
    return dist

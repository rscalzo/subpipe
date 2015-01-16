#!/usr/bin/env perl

# FINDSTARS: Run source detection using Sextractor, pull
# out all the stellar objects

#$snpath = $ENV{ 'SNPATH' } || die " SNPATH ENV VARIABLE NOT FOUND\n";
$snpath="/export/miner1/skymap/brian/SOFTWARE/SNBIN";

# Command line help
if ($#ARGV<0){
    die "USAGE: runsexnew.pl [image]\n".
	"Options:\n".
	"        -outdir dir       directory to write output and intermediate files\n".
	"        -nsigma nsigma    SExtractor detection threshold, nsigma above bkgd\n".
	"        -fwhm fwhm        fwhm in pixels of PSF\n".           
	"        -sat sat          saturation level (ADU) of image\n".
	"        -gain gain        gain (e-/ADU) of image\n".
	"        -threshadu        SExtractor detection threshold, in ADU above back- overrides nsigma)\n".
	"        -zp               Mag offset to apply to\n".
        "        -photo            Use photographic film params\n";
} 

### 0. Read in parameters and make sure everything's kosher.

$image = $ARGV[0];
$DETECTOR="CCD";

die "ERROR: $image not found!\n" unless(-e $image);
$zp=0;
for ($i=1; $i<=$#ARGV; $i++){
    if ($ARGV[$i] eq "-outdir")        { $outdir=$ARGV[++$i]; }
    elsif ($ARGV[$i] eq "-nsigma")     { $nsigma=$ARGV[++$i]; }
    elsif ($ARGV[$i] eq "-fwhm")       { $fwhm=$ARGV[++$i]; }
    elsif ($ARGV[$i] eq "-sat")        { $sat=$ARGV[++$i]; }
    elsif ($ARGV[$i] eq "-gain")       { $gain=$ARGV[++$i]; }
    elsif ($ARGV[$i] eq "-noclean")    { $noclean=1; }
    elsif ($ARGV[$i] eq "-noiter")     { $noiter=1; }
    elsif ($ARGV[$i] eq "-elmax")      { $elongmax=$ARGV[++$i]; }
    elsif ($ARGV[$i] eq "-clmin")      { $classmin=$ARGV[++$i]; }
    elsif ($ARGV[$i] eq "-flagmax")    { $flagmax=$ARGV[++$i]; }
    elsif ($ARGV[$i] eq "-threshadu")  { $threshadu=$ARGV[++$i]; }
    elsif ($ARGV[$i] eq "-zp")         { $zp=$ARGV[++$i]; }
    elsif ($ARGV[$i] eq "-photo")      { $DETECTOR="PHOTO"; }
    else {
	die "ERROR: bad argument <$ARGV[$i]>\n";
    }
}

unless ($gain){
    # Try the header?
    chomp($out = `$snpath/gethead GAIN $image`);
    
    if ($out) {
	$gain=$out;
	print "No gain specified, using GAIN=$gain from header\n";
    } else {
	$gain=1;
	print "WARNING: NO gain supplied, settomg gain=$gain\n";
    }
}

unless ($fwhm)        { $fwhm=3;       print "Setting fwhm to $fwhm pixels\n"; }
unless ($nsigma)      { $nsigma=7;   print "Setting nsigma to $nsigma\n"; }
unless ($sat)         { $sat=30000;    print "Setting saturation to $sat ADU\n"; }
unless ($elongmax)    { $elongmax=1.5; print "Setting elongmax to $elongmax\n"; }
unless ($classmin)    { $classmin=0.3; print "Setting s/g classmin to $classmin \n"; }
unless ($flagmax)     { $flagmax=4;    print "Setting flagmax to $flagmax\n"; }


### Files

($im=$image) =~ s/.*\///;
($impath=$image) =~ s/$ref//;

$outdir = "." unless ($outdir);
$stars  = "${outdir}/${im}.stars";



###  Run SExtractor once ###
$sexout=Sextract($image,$outdir,$fwhm,$sat,$gain,$nsigma,$threshadu);

unless (($sexout =~ /SUCCESS/) && (-e $stars)){
    print STDERR "ERROR: Sextract failed for $image\n";
    exit(0);
}


    





######################################################################
##################### SUBROUTINES ####################################
######################################################################

sub Sextract {
   my ($image,$outdir,$fwhm,$sat,$gain,$sigma,$thresh) = @_;

   # zap the path from the image name
   ($imname=$image) =~ s/.*\///;

   $sexpar = "${outdir}/${imname}.sexpar";
   $defaultconv = "${outdir}/${imname}.conv";
   $defaultnnw = "${outdir}/${imname}.nnw";
   $sexfile  = "${outdir}/${imname}.default_sex";
   $starfile = "${outdir}/${imname}.stars";
   
   # Create default Sextractor files
   makesexfiles($defaultnnw,$defaultconv,$sexpar);

   # RS 2011/04/29:  Modified some of these filenames to correspond
   # to kernels actually existing in Brian's software directory!
   # We should check in any such kernels as part of the pipeline.
   if ($fwhm < 2) {
       system("cp $snpath/gauss_1.5_3x3.conv $defaultconv");
   }
   elsif ($fwhm < 2.5) {
       system("cp $snpath/gauss_2.0_3x3.conv $defaultconv");
   }
   elsif ($fwhm < 3) {
       system("cp $snpath/gauss_2.5_5x5.conv $defaultconv");
   }
   elsif ($fwhm < 3.5) {
       system("cp $snpath/gauss_3.0_5x5.conv $defaultconv");
   }
   elsif ($fwhm < 4) {
       system("cp $snpath/gauss_3.0_7x7.conv $defaultconv");
   }
   elsif ($fwhm < 4.5) {
       system("cp $snpath/gauss_4.0_7x7.conv $defaultconv");
   }
   else{
       system("cp $snpath/gauss_5.0_9x9.conv $defaultconv");
   }

   if ($thresh) {
       $threshtype="ABSOLUTE";
       $sigma=$thresh;
   }
   else {
       $threshtype="RELATIVE";
   }


   ############## create sextractor config file ################
   open (F,">$sexfile");
   print  F <<EOT;

# Default configuration file for SExtractor V1.2
# EB 18/08/97
# (*) indicates parameters which can be omitted from this config file.

#-------------------------------- Catalog ------------------------------------

CATALOG_NAME	$starfile	# name of the output catalog
CATALOG_TYPE	FITS_1.0	# "ASCII_HEAD","ASCII","FITS_1.0" or "FITS_LDAC"

PARAMETERS_NAME	$sexpar	        # name of the file containing catalog contents

#------------------------------- Extraction ----------------------------------

DETECT_TYPE	$DETECTOR	# "CCD" or "PHOTO" (*)
#DETECT_IMAGE	SAME		# "SAME" or <image filename>
FLAG_IMAGE	flag.fits	# filename for an input FLAG-image
DETECT_MINAREA	3		# minimum number of pixels above threshold
DETECT_THRESH	$sigma		# <sigmas> or <threshold>,<ZP> in mag.arcsec-2
ANALYSIS_THRESH	$sigma		# <sigmas> or <threshold>,<ZP> in mag.arcsec-2
THRESH_TYPE     $threshtype

FILTER		Y		# apply filter for detection ("Y" or "N")?
FILTER_NAME	$defaultconv	# name of the file containing the filter

DEBLEND_NTHRESH	32		# Number of deblending sub-thresholds
DEBLEND_MINCONT	0.005		# Minimum contrast parameter for deblending

CLEAN		Y		# Clean spurious detections? (Y or N)?
CLEAN_PARAM	1.0		# Cleaning efficiency

#BLANK		Y		# Blank detected objects (Y or N)?

#------------------------------ Photometry -----------------------------------

PHOT_APERTURES	5		# MAG_APER aperture diameter(s) in pixels
PHOT_AUTOPARAMS	2.5, 3.5	# MAG_AUTO parameters: <Kron_fact>,<min_radius>

SATUR_LEVEL	$sat        	# level (in ADUs) at which arises saturation

MAG_ZEROPOINT	$zp		# magnitude zero-point
MAG_GAMMA	4.0		# gamma of emulsion (for photographic scans)
GAIN		$gain       	# detector gain in e-/ADU.
PIXEL_SCALE	1.0      	# size of pixel in arcsec (0=use FITS WCS info).

#------------------------- Star/Galaxy Separation ----------------------------

SEEING_FWHM	$fwhm      	# stellar FWHM in arcsec
STARNNW_NAME	$defaultnnw	# Neural-Network_Weight table filename

#------------------------------ Background -----------------------------------

BACK_SIZE	16		# Background mesh: <size> or <width>,<height>
BACK_FILTERSIZE	3		# Background filter: <size> or <width>,<height>

BACKPHOTO_TYPE	GLOBAL		# can be "GLOBAL" or "LOCAL" (*)
BACKPHOTO_THICK	24		# thickness of the background LOCAL annulus (*)

#------------------------------ Check Image ----------------------------------

CHECKIMAGE_TYPE	NONE            # can be one of "NONE", "BACKGROUND",
				# "MINIBACKGROUND", "-BACKGROUND", "OBJECTS",
				# "-OBJECTS", "SEGMENTATION", "APERTURES",
				# or "FILTERED" (*)
CHECKIMAGE_NAME	$checkim	# Filename for the check-image (*)

#--------------------- Memory (change with caution!) -------------------------

MEMORY_OBJSTACK	2000		# number of objects in stack
MEMORY_PIXSTACK	1000000		# number of pixels in stack
MEMORY_BUFSIZE	4096		# number of lines in buffer

#---------------- Scanning parameters (change with caution!) -----------------

#SCAN_ISOAPRATIO	0.6		# maximum isoph. to apert ratio allowed (*)

#----------------------------- Miscellaneous ---------------------------------

VERBOSE_TYPE	NORMAL		# can be "QUIET", "NORMAL" or "FULL" (*)

#------------------------------- New Stuff -----------------------------------

# Surprise!!

EOT

  close(F);

   ################### done making config file #######################


   printf STDERR "Doing source detection on $image...";
   @sexoutput=`/usr/local/bin/sex $image -c $sexfile`;
   my $errorflag=0;
   foreach(@sexoutput){
       if ($_ =~ /Background/){
	   chomp;
	   @line = split(/ +/);
	   $background = $line[2];
	   $rms = $line[4];
	   $thresh = $line[7];
       }
       
       if ($_ =~ /Error|ERROR/){
	   printf STDERR "ERROR: $_/n Sextractor didn't work, exiting...\n";
	   exit(0);
       }
   }

   if (!(-e $starfile)){
       printf STDERR "ERROR: Starfile $catalogue not created, exiting\n";
       exit(0);
   }
   printf STDERR "done!\n";

   system "rm -f $defaultconv"    if -e $defaultconv;
   system "rm -f $defaultnnw"     if -e $defaultnnw;
   system "rm -f $sexfile"        if -e $sexfile;
   system "rm -f $sexpar"         if -e $sexpar;
   
   printf "Sextracted objects in $starfile\n";
   return "back: $background rms: $rms thresh: $thresh SUCCESS! ";
   
}

sub cleanstars{
# input: starlistname FLAG_MAX ELONG_MAX CLASSIFIER_MIN STARSUSED_MIN STARSUSED_MAX
# if number stars < STARSUSED_MIN, original stars file is NOT overwritten,
# otherwise, original stars files is copied to starlistname.old, and cleaned 
# list is written into starlistname.
# if STARSUSED_MAX=0 then no upper limit
# output: @cleanstarlist, 
 
 my $starlistname=$_[0];
 my $FLAG_MAX=$_[1];
 my $ELONG_MAX=$_[2];
 my $CLASSIFIER_MIN=$_[3];
 my $STARSUSED_MIN=$_[4];
 my $STARSUSED_MAX=$_[5];
 my $dummy=$ENV{"PWD"};


 chomp(@starlist=`cat $starlistname`);
  
 undef my @cleanstarlist;

 if ($verbose) {printf STDERR "$#starlist stars in list\n";}

 my $meanfwhm = 0;
 my $nstars = 0;
 my $ngoodstars = 0;
 foreach (@starlist){
     chomp($_);
     @star=split(" ",$_);
     if ($star[0] eq "#") {
	 push(@cleanstarlist,$_);
     } elsif (@star>0) {
	 my $num = $star[0];
	 my $x=$star[1];
	 my $y=$star[2];
	 my $mag=$star[3]+25;
	 
	 my $flag = $star[4];
	 my $elong = $star[8];
	 my $fwhm = $star[9];
	 my $classifier = $star[10];
	 $nstars++;

	 # check flags.  if it passes, save it
	 # also, keep track of the mean FWHM

	 if (($flag<=$FLAG_MAX) && ($elong<=$ELONG_MAX) && 
	     ($classifier>=$CLASSIFIER_MIN) && ($fwhm > 0) &&
	     (($STARSUSED_MAX==0) || ($ngoodstars<$STARSUSED_MAX))){
	     
	     $ngoodstars++;
	     push(@cleanstarlist,$_);
	     $meanfwhm = (1./$ngoodstars)*(($ngoodstars-1)*$meanfwhm + $fwhm);
	     push(@fwhm,$fwhm);
	 }
     }
 }

 # calculate sigma-clipped mean fwhm

 foreach $fwhm (@fwhm) { $sum = ($fwhm - $meanfwhm)*($fwhm - $meanfwhm); }

 if ($ngoodstars > 1){
     $sigma = sqrt( $sum / ($ngoodstars - 1));  $n=0;
     foreach $fwhm (@fwhm) {
	 if (($fwhm < ($meanfwhm+3*$sigma)) && ($fwhm > $meanfwhm-3*$sigma)){
	     $n++;
	     $clippedmeanfwhm = (1./$n)*(($n-1)*$clippedmeanfwhm + $fwhm);
	 }
     }
 } else {
     printf "too few stars to calculate mean fwhm\n";
     $clippedmeanfwhm = $meanfwhm;
 }

 if ($ngoodstars>=$STARSUSED_MIN) {
     $outstr = join("\n",@cleanstarlist);
     open(F,">$starlistname");
     print F $outstr;
     close F;
     printf STDERR "Cleaned $starlistname - $ngoodstars/$nstars kept\n";
 } else {
     printf STDERR "Only $ngoodstars remain - leaving $starlistname unchanged\n";
 }

 return $clippedmeanfwhm;

}


sub wcsmatch{
    my ($xref_ptr,$yref_ptr,$xshft_ptr,$yshft_ptr,$xwcs_ptr,$ywcs_ptr,$dex_ptr,$nx,$ny,$matchlist,$searchrad) = @_;

    $searchrad2=$searchrad*$searchrad;

    # Do matching

    undef my @matches;undef my @nomatch;undef my @ambiguousmatch;
    $Ngood=0;$Nbad=0;$Nambiguous=0;$Ntot=0;


    # Loop through (wcs) shft objects

    for($i=0; $i<@$xwcs_ptr; $i++){    
	my $xw = $$xwcs_ptr[$i];
	my $yw = $$ywcs_ptr[$i];
	
	$Ntot++;
	
	my $Nfound=0;
	my $bestdistance2=$searchrad2;my $bestindex=-1;


	my $startx = int($xw-$searchrad);
	$startx = 0 if ($startx < 0);
	$startx = $nx if ($startx > $nx);
	my $jmin=$$dex_ptr[$startx];

	for($j=$jmin; $j<@$xref_ptr; $j++){
	    
	    my $xr = $$xref_ptr[$j];
	    my $yr = $$yref_ptr[$j];
	    
	    # Don't search further than you need...
	    last if ($xr>$xw+$searchrad);
	    
	    $distance2=($xr-$xw)*($xr-$xw) + ($yr-$yw)*($yr-$yw);

	    if ($distance2 < $searchrad2){
		$Nfound++;
		if ($distance2<$bestdistance2){
		    $bestindex=$j;
		    $bestdistance2=$distance2;	
		}
	    }
	    
	}
	
	# Save xref, yref and xshft, yshft of good matches
	if ($Nfound>0){
	    if ($Nfound==1){
		$Ngood++;
		$matchstr .= "$$xshft_ptr[$i] $$yshft_ptr[$i] ".
		    "$$xref_ptr[$bestindex] $$yref_ptr[$bestindex]\n";

	    } else { 
		$Nambiguous++; 
	    }
	} else { 
	    $Nbad++; 
	}
	
	if (($Ntot % 500 == 0) && ($Ntot>0)) {
	    printf STDERR sprintf("%5d: %5d(%3.1f\%) matched %5d(%3.1f\%) unmatched\n",$Ntot,$Ngood,$Ngood/$Ntot*100.,$Nbad ,$Nbad/$Ntot*100);
	}
    }
    print STDERR sprintf("total: %5d(%3.1f\%) matched %5d(%3.1f\%) unmatched\n",$Ngood,$Ngood/$Ntot*100,$Nbad ,$Nbad/$Ntot*100);
    print STDERR sprintf("total matches: %5d\n",$Ngood);

    return $matchstr;
}

sub pointer2sortedlist{
    my ($nx,@xlist)=@_;
    # create pointers to the sorted lists: makes matching faster...
    # $cmpindex[$xx] is the index of the first measurement, for which x>=$xx
    undef my @cmpindex;
    my $Ncmpobj=@xlist;

    # first, set all pointers to 0
    for(my $x=0;$x<=$nx;$x++){$cmpindex[$x]=0;}

    my $xlast=-1;
    my $lastgoodindex=0;
    for($i=0;$i<$Ncmpobj;$i++){
	my $x= int($xlist[$i]);
	next if ($x < 0.0); # skip negative values. This is taken care of by setting $cmpindex[0]=0 later on
	next if ($x > $nx);
	$lastgoodindex=$i;
	next if ($x == $xlast);
	for (my $xx=$xlast+1;$xx<=$x;$xx++){
	    $cmpindex[$xx]=$i;
	}
	$xlast=$x;
	
    }

    $cmpindex[0]=0; # make the first index always point to the first entry. Allows negative X values!
 
    # fillup the rest of cmpindex
    $lastgoodindex++ if ($lastgoodindex+1<$Ncmpobj);
    for (my $xx=$xlast+1;$xx<=$nx;$xx++){
	$cmpindex[$xx]=$lastgoodindex;
    }
    
    return(@cmpindex);
}

sub xsort{
    my ($x1) = $a =~ /^\s*(\S+)/;
    my ($x2) = $b =~ /^\s*(\S+)/;
    $x1 <=> $x2;
}

sub numerically { $a <=> $b; }


sub getparamval {
    # look for keyword in parameter list
    # if it doesn't exist, return the default

    my ($keyword,$default,@plist)=@_;
    undef my $pval;
    
    unless ($#plist<0) {
	$pval = getparamfromlist($keyword,@plist);	
	if ($pval){ return $pval; }
    }

    return $default;
}
    


sub getparamfromlist{
    my ($keyword,@list)=@_;
    my $returnstring="";

    foreach $line (@list)  {
	if ($line =~ /$keyword/) {
	    my @pair = split (/\s+/, $line);
	    $Nline=@pair;
	    if ($Nline>0) {
		if ($pair[0] eq $keyword) {
		    $returnstring=$pair[1];
		    last;
		}
	    }
	}
    }
    return $returnstring;
}	


sub do_xform {
    my ($matchlist,$map,$halfpix,$order,$maxrms)=@_;
    my $success=0;
    
    # Go from order of polynomial to number of terms at that
    # order, which is how xform wants it.
    if    ($order==1){ $order=3; }
    elsif ($order==2){ $order=6; }
    elsif ($order==3){ $order=10;}
    else { die "ERROR: don't know how to fit order $order!\n"; }

    # If halfpix is set, do NOT pass halfpix flag to xform
    # Confusing, I know, sorry...
    if ($halfpix){ $halfflag = ""; }
    else { $halfflag = "-halfpix";}

    $xformcmd = "xform $matchlist $map $halfflag $order";
    print "$xformcmd\n";
    my @xformoutput;
#    if ($Verbose) {  @xformoutput=`$xformcmd 2>&1`; }  
#    else { @xformoutput=`$xformcmd`; }
    @xformoutput=`$xformcmd 2>&1`; 


    ### get the rms of the astrometry
    my $rms=-1;
    for ($i=0;$i<scalar(@xformoutput);$i++){
	print "XFORM: $xformoutput[$i]";
	if ($xformoutput[$i]=~/rms =\s+(\S+)/){
	    chomp($rms=$1);
	    last;
	}
    }

    printf(stderr "RMS: $rms MAX RMS: $maxrms\n");	
    
    if (($rms > $maxrms) || ($rms eq "nan") || ($rms < 0)){
	printf(stderr "ERROR: rms $rms is not acceptable!\n");
	$success=0;
    } else {
	printf(stderr "ASTRO RMS IS GOOD!\n");
	$success=1;
    }
    return($success);
}

sub do_remap{
    my $shftim = $_[0];  # image to get remapped
    my $remapim =$_[1];  # output image of remapping
    my $mapping =$_[2];  # geometric transformation from xform
    my $nx      =$_[3]; 
    my $ny      =$_[4]; 
    my $bitpix  =$_[5];  # bits per pixel
    my $sat     =$_[6];  # saturation level (ADU)
    my $maskval =$_[7];  # value of masked pixels
    my $noiseim =$_[8];  # noise image created if filename is passed
    my $nbscl   =$_[9];  # BSCALE of noise image
    my $gain    =$_[10]; # e-/ADU, need this to get noise right
    my $rdnoise =$_[11]; # read noise in e-

    my $cmd = "remap -src $shftimage -im $remapimage -map $mapping -nx ".
	"$nx -ny $ny -bitpix $bitpix -satval $sat -maskval $maskval";
    
    if ($noiseim){
	$cmd .= " -noise $noiseim -nbscl $nbscl -gain $gain -rd $rdnoise";
    }
	    
    #  Do remapping
    print "$cmd\n";
    my @remapoutput = `$cmd 2>&1`;

    # Did we make it?
    die "ERROR: remap failed!" if("@remapoutput"=~/ERROR/);
    die "ERROR: failed to create image" if (!(-e $remapim));

    # Guess so...
    return "Remapping successful";
}


sub copywcs {
    my ($wcs_fits,$file2wcs) = @_; 

    my $result=0;
    #remove multiple slashes because iraf sucks
    while($wcs_fits =~ s/\/\//\//g){;}
    while($file2wcs =~ s/\/\//\//g) {;}
    
    my $cl_file   = "cl".$$.".file";
    unlink $cl_file if -e $cl_file; 
    open(CLFILE, ">$cl_file");  
    my $pwd = $ENV{'PWD'};
    
    my @split1;my @split2;
    undef @split1;
    undef @split2;
    @split1 = split(/\//,$wcs_fits);
    @split2 = split(/\//,$file2wcs);
    
    my $i;
    my $min_len;
    if($#split1 < $#split2)  { $min_len =  $#split1;}
    else {$min_len =  $#split2;}
    my $common_path;
    
    ### search for common path (try to reduce iraf filename!!!) 
    #no path, use current working direcotry!   

    if($#split1 <= 0  && $#split2 <= 0 ) {
	$common_path = ".";
    }
    else {
	for($i=1;$i<$min_len;$i++) {
	    if($split1[$i] eq $split2[$i]) {
		$common_path = $common_path."/".$split1[$i];
	    }
	    else {
		last;
	    }
	}
    }
    
    my $wcs_fits_cpy = $wcs_fits;
    my $file2wcs_cpy = $file2wcs;

    if($common_path eq "") {  $common_path = "."; }
    if($common_path ne ".") {
	$wcs_fits =~ s/$common_path//;
	$file2wcs =~ s/$common_path//;
	$wcs_fits =~ s/^\//\.\//;
	$file2wcs =~ s/^\//\.\//;
    }
    else {
	$common_path = $pwd;
    }
    
    print "************************comm: $common_path wcs_fits:$wcs_fits file2wcs:$file2wcs\n";
       
    printf CLFILE "\n";
    printf CLFILE "print hello\n";
    printf CLFILE "cd $common_path\n";
    printf CLFILE "wcscopy $file2wcs $wcs_fits\n";
    printf CLFILE "hedit.add = yes\n";
    printf CLFILE "hedit.addonly = no\n";
    printf CLFILE "hedit.delete = no\n";
    printf CLFILE "hedit.verify = no\n";
    printf CLFILE "hedit.show = yes\n";
    printf CLFILE "hedit.update = yes\n";
    printf CLFILE "hedit $file2wcs HISTORY 'IRAF: wcscopy $file2wcs $wcs_fits'\n";
    
    printf CLFILE "logout\n";
    close CLFILE;
    printf STDERR "Copying IRAF-specify WCS data FROM:$wcs_fits  TO:$file2wcs using cl script $cl_file\n\n";
    my @clout = &cl_run2($cl_file);

    unlink $cl_file if -e $cl_file;
    
    foreach $line (@clout){
	if ((grep /ERROR|Error|fail|exception/, $line)){
	    printf STDERR "ERROR COPYWCS : something went wrong with copywcs during IRAF-stage(wcscopy or hedit)!!!!\n";
	    printf STDERR "OFFENDING LINE: $line\n";
	    return(0);
	}	
    }
    
#### NOTE: IRAF's WCSCOPY does NOT seem to copy over CRPIX, but actually deletes it, so I use some wcstools to copy those.
    #return the names to original ones...
    $wcs_fits = $wcs_fits_cpy;
    $file2wcs = $file2wcs_cpy;
    printf STDERR "Copying standard WCS data FROM:$wcs_fits  TO:$file2wcs\n\n"; 
    #get standard-WCS info
    #(my $wcs_header1) = "ASTROCAT EWCSXAS EWCSYAS WCSASTRM WCSDIM EQUINOX CTYPE1 CRVAL1  CRPIX1  CTYPE2 CRVAL2  CRPIX2  CD1_1 CD2_1 CD1_2 CD2_2";
    # Kick out WCSASTRM, since it gives segmentation fault!
    (my $wcs_header1) = "ASTROCAT EWCSXAS EWCSYAS WCSDIM EQUINOX CTYPE1 CRVAL1  CRPIX1  CTYPE2 CRVAL2  CRPIX2  CD1_1 CD2_1 CD1_2 CD2_2";
    
    printf STDERR "Executing:gethead -e $wcs_fits $wcs_header1 \n";
    (my $wcs_info1) = `gethead -te $wcs_fits $wcs_header1`; # separate by tabs
    $wcs_info1 =~ s/\s+$//g; # remove leading and trailing whitespace
    $wcs_info1 =~ s/^\s+//g;
    
    (my @fields1) = split(/\t/,$wcs_info1); # split on tabs
    my $i;
    undef $wcs_info1;
    for($i=0;$i<=$#fields1;$i++) {
	$fields1[$i] =~ s/\s+$//g;
	$fields1[$i] =~ s/^\s+//g;
	# add single-quotes around the 'value' 
	$fields1[$i] =~ s/\=/\=\'/;
	$fields1[$i] = $fields1[$i]."'";
	$wcs_info1 = $wcs_info1." ".$fields1[$i];
    }
    
    printf STDERR "Executing:sethead -hv $wcs_info1  $file2wcs\n";
    (my @output) = `sethead -hv $wcs_info1  $file2wcs`;
    unless((grep /success/, @output)){
	printf STDERR "ERROR in wcstool's  sethead, when copying standard WCS...exiting...\n\n";
	return(0);
    }
#    printf STDERR "copywcs: SUCCESS\n";
    return(1);
}


sub cl_run2 {
   my ($iraf_script) = @_;
   my $cur_dir = $ENV{'PWD'};
   my $iraf_home = $ENV{'HOME'}."/iraf";

   chdir("$iraf_home");

   printf "Running IRAF...\n";

   if (!(-e "$cur_dir$,/$iraf_script")) {
       die "$cur_dir$,/$iraf_script does not exist\n";
   }

   my @clout = `cl < $cur_dir$,/$iraf_script`;
   chdir("$cur_dir");
   
   my $myflag = 0;
   my $cline;
   for $cline (@clout) {
      $myflag = 1 if $cline =~ /Welcome to IRAF/;
      
      printf "$cline" if $myflag;
   }

   unlink $iraf_script if -e $iraf_script;
   return @clout;
}


###### WARNING!
###### WARNING!
###### WARNING!
### This is a bit ugly.  SExtractor needs a few files in order to
### run properly (neural network weights, convolution kernel, etc.)
### Since these can be located in different places depending on where 
### SExtractor has been installed, I have chosen to have wcsalign
### automatically create them.  Do NOT edit or muck with the
### following subroutine unless you are a SExtractor guru. You 
### have been warned,

sub makesexfiles{
    
    ($nnwfile, $convfile, $outformatfile) = @_;

    ############# 1. Neural Network Weights ####################
    unlink $nnwfile;
    open(NNW,">$nnwfile") || die "ERROR: couldn't open $nnwfile!\n";

    print NNW <<EOT;
NNW
# Neural Network Weights for the SExtractor star/galaxy classifier (V1.3)
# inputs:	9 for profile parameters + 1 for seeing.
# outputs:	``Stellarity index'' (0.0 to 1.0)
# Seeing FWHM range: from 0.025 to 5.5'' (images must have 1.5 < FWHM < 5 pixels)
# Optimized for Moffat profiles with 2<= beta <= 4.

 3 10 10  1

-1.56604e+00 -2.48265e+00 -1.44564e+00 -1.24675e+00 -9.44913e-01 -5.22453e-01  4.61342e-02  8.31957e-01  2.15505e+00  2.64769e-01
 3.03477e+00  2.69561e+00  3.16188e+00  3.34497e+00  3.51885e+00  3.65570e+00  3.74856e+00  3.84541e+00  4.22811e+00  3.27734e+00

-3.22480e-01 -2.12804e+00  6.50750e-01 -1.11242e+00 -1.40683e+00 -1.55944e+00 -1.84558e+00 -1.18946e-01  5.52395e-01 -4.36564e-01 -5.30052e+00
 4.62594e-01 -3.29127e+00  1.10950e+00 -6.01857e-01  1.29492e-01  1.42290e+00  2.90741e+00  2.44058e+00 -9.19118e-01  8.42851e-01 -4.69824e+00
-2.57424e+00  8.96469e-01  8.34775e-01  2.18845e+00  2.46526e+00  8.60878e-02 -6.88080e-01 -1.33623e-02  9.30403e-02  1.64942e+00 -1.01231e+00
 4.81041e+00  1.53747e+00 -1.12216e+00 -3.16008e+00 -1.67404e+00 -1.75767e+00 -1.29310e+00  5.59549e-01  8.08468e-01 -1.01592e-02 -7.54052e+00
 1.01933e+01 -2.09484e+01 -1.07426e+00  9.87912e-01  6.05210e-01 -6.04535e-02 -5.87826e-01 -7.94117e-01 -4.89190e-01 -8.12710e-02 -2.07067e+01
-5.31793e+00  7.94240e+00 -4.64165e+00 -4.37436e+00 -1.55417e+00  7.54368e-01  1.09608e+00  1.45967e+00  1.62946e+00 -1.01301e+00  1.13514e-01
 2.20336e-01  1.70056e+00 -5.20105e-01 -4.28330e-01  1.57258e-03 -3.36502e-01 -8.18568e-02 -7.16163e+00  8.23195e+00 -1.71561e-02 -1.13749e+01
 3.75075e+00  7.25399e+00 -1.75325e+00 -2.68814e+00 -3.71128e+00 -4.62933e+00 -2.13747e+00 -1.89186e-01  1.29122e+00 -7.49380e-01  6.71712e-01
-8.41923e-01  4.64997e+00  5.65808e-01 -3.08277e-01 -1.01687e+00  1.73127e-01 -8.92130e-01  1.89044e+00 -2.75543e-01 -7.72828e-01  5.36745e-01
-3.65598e+00  7.56997e+00 -3.76373e+00 -1.74542e+00 -1.37540e-01 -5.55400e-01 -1.59195e-01  1.27910e-01  1.91906e+00  1.42119e+00 -4.35502e+00

-1.70059e+00 -3.65695e+00  1.22367e+00 -5.74367e-01 -3.29571e+00  2.46316e+00  5.22353e+00  2.42038e+00  1.22919e+00 -9.22250e-01 -2.32028e+00


 0.00000e+00 
 1.00000e+00 

EOT

close NNW;


    ############### 2. Convultion filter #########################
    unlink $convfile;
    open(CF,">$convfile") || die "ERROR: couldn't open $convfile!\n";
    print CF <<EOT;
CONV NORM
# 3x3 ``all-ground'' convolution mask with FWHM = 2 pixels.
1 2 1
2 4 2
1 2 1

EOT
close CF;


    ################ 3.Output format file #######################
    unlink $outformatfile;
    open(OF,">$outformatfile") || die "ERROR: couldn't open $outformatfile!\n";
    print OF <<EOT;
NUMBER
X_IMAGE
Y_IMAGE
ALPHA_J2000
DELTA_J2000
MAG_BEST
FLAGS
A_IMAGE
B_IMAGE
THETA_IMAGE
ELONGATION
FWHM_IMAGE
CLASS_STAR
EOT
close OF;


}

#!/usr/bin/env perl

$brianwcs = $ENV{ 'BRIANWCS' } || die "Whoops!  Please setenv BRIANWCS and try again";

require "skysubroutines.pl";

if ( $#ARGV < 0 ) {
    print( "WCS.pl: image1 image2 ... image for WCS to be applied\n scale is in arcseconds per pixel\n scaleuncertainty is uncertainty in arcsec per pixel of scale (0.005 is default)\n");
    print("\n\t-scale \n\t-tol scaleuncertainty \n\t-h headertype \n\t-m Mosaic size \n\t-o iraforder \n\t-stars starmatches \n\t-match #starstomatchatonce \n\t-amiss constant_tol_in_fit");
    print("\n\t-cbmiss scale_tol_in_fit\n\t-headertype RA DEC format (H=hours DH=decimal hours R=radians D=decimal degrees) (H is default)");
    print("\n\t-iraforder (2-5) determine auto is default\n\t-zpB (calculate zpB for data (GSC))\n\t-zpR calculate zpR for data (UNSNO)");
    print("\n\t-usewcs (use existing WCS as starting point to get RA, DEC and scale)\n\t-tycho (use TYCHO Catalog)\n\t-sao (USE SAOcatalog)\n\t-macho (USE macho LMCcatalog)");
    print("\n\t-macs (USE MACS LMCcatalog)\n\t-localcat (USE a local user-defined catalog)\n\t-mag (bright limit)\n\t-interac (make interactive)\n\t-roterr [45] allowed rotation error in degrees");
    print("\n\t-ucac2 (USE ucac2 catalog)\n\t-2mass (use 2MASS Catalog - only on miner right now)");
    print("\n\t-donotapply Match, but do not apply fit to header");
    print("\n\t-thresh threshold for SExtractor detection");
    print("\n\t-mask maskname \n\n");
    die;
}
($img) = @ARGV;
$donotapply=0;
$interac="";
$sao=0;
$tycho=0;
$ucac1=0;
$twomass=0;
$ucac2=0;
$roterr="-roterr 2";
$macho=0;
$macs=0;
$localcat="";
$miss=3.0;
$zpR=0;
$zpB=0;
$tolinit=0.005;
$headertype="H";
$Mosaic=1;
$iraforder=0;
$starmatchesinit=8;
$maxatonceinit=50;
$amiss=10;
$cbmiss=0.01;
$m1=0.;
$m2=22;
$force=0;
$usewcs=0;
$scaleinit="ND";
$thresh=10.;
$clean=0.2;
$cleansig=3;
$maskname="";
$usemask=0;

#get list of files
foreach(@ARGV) {
    if(/^\-/) {
	last;
    }
    else {
	$files[$#files+1]=$ARGV[$#files+1];
    }
}
	

for ($i = $#files+1; $i <= $#ARGV; $i++) {
    if ($ARGV[$i] eq '-o') {
	$i++;
	$iraforder=$ARGV[$i];
    }
    elsif ($ARGV[$i] eq '-tol') {
	$i++;
	$tolinit=$ARGV[$i];
    }
    elsif ($ARGV[$i] eq '-interac') {
	$interac="-interac";
    }
    elsif ($ARGV[$i] eq '-donotapply') {
	$donotapply=1;
    }
    elsif ($ARGV[$i] eq '-ucac1') {
	$ucac1=1;
    }
    elsif ($ARGV[$i] eq '-2mass') {
	$twomass=1;
    }
    elsif ($ARGV[$i] eq '-ucac2') {
	$ucac2=1;
    }
    elsif ($ARGV[$i] eq '-macho') {
	$macho=1;
    }
    elsif ($ARGV[$i] eq '-macs') {
	$macs=1;
    }
    elsif ($ARGV[$i] eq '-localcat') {
	$i++;
	$localcat=$ARGV[$i];
    }
    elsif ($ARGV[$i] eq '-sao') {
	$sao=1;
    }
    elsif ($ARGV[$i] eq '-tycho') {
	$tycho=1;
    }
    elsif ($ARGV[$i] eq '-mag') {
	$i++;
	$m2=$ARGV[$i];
    }
    elsif ($ARGV[$i] eq '-h') {
	$i++;
	$headertype=$ARGV[$i];
    }
    elsif ($ARGV[$i] eq '-m') {
	$i++;
	$Mosaic=$ARGV[$i];
    }
    elsif ($ARGV[$i] eq '-force') {
	$force=1;
    }
    elsif ($ARGV[$i] eq '-thresh') {
	$i++;
	$thresh=$ARGV[$i];
    }
    elsif ($ARGV[$i] eq '-match') {
	$i++;
	$maxatonceinit=$ARGV[$i];
    }
    elsif ($ARGV[$i] eq '-stars') {
	$i++;
	$starmatchesinit=$ARGV[$i];
    }
    elsif ($ARGV[$i] eq '-scale') {
	$i++;
	$scaleinit=$ARGV[$i];
    }
    elsif ($ARGV[$i] eq '-amiss') {
	$i++;
	$amiss=$ARGV[$i];
    }
    elsif ($ARGV[$i] eq '-zpB') {
	$zpB=1;
    }
    elsif ($ARGV[$i] eq '-zpR') {
	$zpR=1;
    }
    elsif ($ARGV[$i] eq '-clean') {
	$i++;
	$clean=$ARGV[$i];
    }
    elsif ($ARGV[$i] eq '-cleansig') {
	$i++;
	$cleansig=$ARGV[$i];
    }
    elsif ($ARGV[$i] eq '-usewcs') {
	$usewcs=1;
    }
    elsif ($ARGV[$i] eq '-cbmiss') {
	$i++;
	$cbmiss=$ARGV[$i];
    }
    elsif ($ARGV[$i] eq '-roterr') {
	$i++;
	$roterr="-roterr $ARGV[$i]";
    }
    elsif ($i==1 && ($ARGV[$i] >0 && $ARGV[$i] < 10000)) {
	$scale=$ARGV[$i]*1;
	#for compatability with previous versions allow first thing to be scale
    }
    elsif ($ARGV[$i] eq '-mask') {
	$i++;
	$maskname=$ARGV[$i];
	$usemask=1;
    }
    else {
	die "command $ARGV[$i], not understood\n";
    }
}

foreach(@files) {

    # initialise variables which we clobber
    #
    #
    $scale=$scaleinit;
    $tol=$tolinit;
    $starmatches=$starmatchesinit;
    $maxatonce=$maxatonceinit;

    $img=$_;
    $_=$img;
    if (/(\S+).fits/) {$img=$1;}

#no star list for image..Run Sexextract to find stars
    $stars1="$img.fits.xieta";
    $stars2="$img.fits.stars.txt";
    $matchlist="$img.fits.wcsmatch";


    if (! -e $stars2 || $force) {
	if ($usemask==1) {
	    print "$brianwcs/bin/runsex.pl $img.fits $thresh -mask $maskname\n";
	    system("$brianwcs/bin/runsex.pl $img.fits $thresh -mask $maskname");
	    system("mv $img.fits.stars $stars2")
	}
	else {
	    print "$brianwcs/bin/runsex.pl $img.fits $thresh\n";
	    system("$brianwcs/bin/runsex.pl $img.fits $thresh");
	    system("mv $img.fits.stars $stars2")
	}
	#copy over flag 0 objects
	#open(F,"<$img.fits.stars");
	#open(OUT,">$stars2");
	#while (<F>) {
	#   $line=$_;
	#   ($XX[$i],$YY[$i],$flag,$fwhm[$i++])=/^\s*\S+\s+(\S+)\s+(\S+)\s+\S+\s+(\S+)\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+)/;
	#   if ($flag ne "0" ) {
	#	$i--;
	#    } else {
	#	printf(OUT $line);
	#    }   
	#}
	#close(OUT);
	#close(F);	
    }
    
    open(F,"<$stars2");
    $i=0;
    while(<F>) {
	($XX[$i],$YY[$i],$flag,$fwhm[$i++])=/^\s*\S+\s+(\S+)\s+(\S+)\s+\S+\s+(\S+)\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+)/;
	if ($flag ne "0" ) {
	    $i--;
	}
    }
    close(F);
    $Nstars=$i;
    
#adjust match parameters
    print "Number of good stars: $Nstars \n";
    if ($Nstars < 60) {
	$starmatches = 5;
    }

#Get header info
    

    $_ = `gethead $img.fits RA DEC NAXIS1 NAXIS2 CTYPE1 CTYPE2`;
    ($RA,$DEC,$xaxis,$yaxis,$ctype1,$ctype2) = split;
    if ($usewcs==0) {
	
	if ($RA eq "___") {
	    $_ = `gethead $img.fits MEANRA MEANDEC`;
	    ($RA,$DEC) = split;
	}
	if ($RA eq "___" || $DEC eq "___" ) {
	    die "No RA/DEC in header\n";
	}
	
	
	if ($headertype eq "H" ) {
	    ($RAr,$DECr)=radecstring2rad($RA,$DEC);
	    printf ("Read from Header:");
	}
	if ($headertype eq "DH" ) {
	    $RAr=$RA/24*2*3.1415926535;
	    $DECr=$DEC/360*2*3.1415926535;
	    ($RA,$DEC)=radecrad2string($RAr,$DECr);
	    printf ("Read from Header: ");
	}
	
	if ($headertype eq "R") {
	    $RAr=$RA;
	    $DECr=$DEC;
	    ($RA,$DEC)=radecrad2string($RAr,$DECr);
	    printf ("Read from Header: ");
	}
	
	if ($headertype eq "D") {
	    $RA=$RA/15.;
	    $DEC=$DEC;
	    $RAr=$RA/24*2*3.1415926535;
	    $DECr=$DEC/360*2*3.1415926535;
	    printf ("Read from Header: ");
	}
	$_=$ctype1;
	if (/RA/ && $force !=1 && $usewcs !=1) {
	    $_=$ctype2;
	    if (/DEC/){
		print  "WCS already present in image. Use -force to force update. Or -usewcs\n";
		next;
	    }
	}
    }
    else { #use WCS to get center RA and DEC
	$midx=$xaxis/2;
	$midy=$yaxis/2;
	$_=`xy2sky $img.fits $midx $midy`;
	($RA,$DEC)=/^(\S+)\s+(\S+)/;
	($RAr,$DECr)=radecstring2rad($RA,$DEC);
	$midx+=100; #offset by 100 pixels to get scale
	$midy+=100;
	$_=`xy2sky $img.fits $midx $midy`;
	($RA1,$DEC1)=/^(\S+)\s+(\S+)/;
	($RAr1,$DECr1)=radecstring2rad($RA1,$DEC1);
	($dtheta,$dphi)=eq2st($RAr,$DECr,$RAr1,$DECr1);
	if($scale eq "ND") {
	    $scale=sqrt($dtheta*$dtheta+$dphi*$dphi)/sqrt(20000);
	}
    }
    printf ("RA DEC read from Header %s %s\n",$RA,$DEC);
#if scale is undefined, try the header...
    if ($scale eq "ND") {
	$_=`gethead $img.fits SECPIX SECPIX1 PIXSCAL1`;
	($secpix,$secpix1,$pixscal1)=/(\S+)\s+(\S+)\s+(\S+)/;
	if ($secpix >0) {
	    $scale =$secpix;
	    print "SCALE Read from Header...SECPIX=$scale arcsec/pixel\n";
	}
	elsif ($secpix1 > 0) {
	    $scale =$secpix1;
	    print "SCALE Read from Header...SECPIX1=$scale arcsec/pixel\n";
	}
	elsif ($pixscal1 > 0) {
	    $scale =$pixscal1;
	    print "SCALE Read from Header...PIXSCAL1=$scale arcsec/pixel\n";
	}
	else { #no info..do our best to figure it out.
	    $scale=0.0;
	}
    }
    printf ("RA:$RA DEC:$DEC Scale:%5.3f arcsec/pixel\n",$scale);
    $size1=$xaxis*$scale*1.2*$Mosaic;
    $size2=$yaxis*$scale*1.2*$Mosaic;
    if ($size1 > $size2) {$gscsize=$size1/2;}
    else {$gscsize=$size2/2};
    
    
    $OKmatch=sqrt(5*$scale*5*$scale+$miss*$miss);
    
    
#no gsc list Extract list
# if forced then remove existing img.fits.catalog
    if ($force) {
	unlink ("$img.fits.catalog");
    }
    if (! -e "$img.fits.catalog") {
	$usnocat="usa2";
	if ($ucac1) {
	    print "$brianwcs/bin/readucac1.pl $img.fits.catalog $RA $DEC $gscsize\n";
	    system("$brianwcs/bin/readucac1.pl $img.fits.catalog $RA $DEC $gscsize");
	}
	if ($twomass) {
	    print "$brianwcs/bin/read2mass.pl $img.fits.catalog $RA $DEC $gscsize\n";
	    system ("$brianwcs/bin/read2mass.pl $img.fits.catalog $RA $DEC $gscsize");
	}
	elsif ($ucac2) {
	    print "$brianwcs/bin/readucac2.pl $img.fits.catalog $RA $DEC $gscsize\n";
	    system("$brianwcs/bin/readucac2.pl $img.fits.catalog $RA $DEC $gscsize >/dev/null");
	}
	elsif ($macho) {
	    print "$brianwcs/bin/readmacho.pl $img.fits.catalog $RA $DEC $gscsize\n";
	    system("$brianwcs/bin/readmacho.pl $img.fits.catalog $RA $DEC $gscsize");
	}
	elsif ($macs) {
	    print "$brianwcs/bin/readmacs.pl $img.fits.catalog $RA $DEC $gscsize\n";
	    system("$brianwcs/bin/readmacs.pl $img.fits.catalog $RA $DEC $gscsize");
	}
	elsif ($localcat ne "") {
	    if (! -e "$localcat") {
		die "Local catalog $localcat does not exist\n";
	    }
	    else {
		system("cp $localcat $img.fits.catalog");
	    }
	}
	elsif ($tycho) {
	    if ( $ENV{'TY2_PATH'} eq "") { 
		die "TY2_PATH not found\n";
	    }
	    print "scat -m $m1 $m2 -c ty2 -r $gscsize -n 1000 $RA $DEC 2000.0 | sort -k 5 > $img.fits.catalog\n";
	    system("scat -m $m1 $m2 -c ty2 -r $gscsize -n 1000 $RA $DEC 2000.0 | sort -k 5 > $img.fits.catalog");
	}
	elsif ($sao) {
	    if ( $ENV{'SAO_PATH'} eq "") { 
		die "SAO_PATH not found\n";
	    }
	    print "scat -m $m1 $m2 -c sao -r $gscsize -n 1000 $RA $DEC 2000.0 > $img.fits.catalog\n";
	    system("scat -m $m1 $m2 -c sao -r $gscsize -n 1000 $RA $DEC 2000.0 > $img.fits.catalog");
	}
	else {
	    if ( $ENV{'TMC_PATH'} eq "") { 
		die "TMC_PATH variable not set - cannot use this as the default catalog - nothing else asked for\n";
	    }
	    print "scat -m $m1 $m2 -c tmc -r $gscsize -n 1000 $RA $DEC 2000.0 | sort -k 5 >> $img.fits.catalog\n";
	    system("scat -m $m1 $m2 -c tmc -r $gscsize -n 1000 $RA $DEC 2000.0 | sort -k 5 >> $img.fits.catalog");
	}
    }
    
#make a magnitude sorted list on tangent plane
    open(F,">$img.fits.xieta");
    open(S,"$img.fits.catalog"); #open catalog
    $count=0;
    while(<S>) {
	if (! /#/) {
	    ($id1,$RA1,$DEC1,$mag)=/^\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/;
	    if ($RA1 ne "") {
		($RA1r,$DEC1r)=radecstring2rad($RA1,$DEC1); #go to radians	
		($x,$y)=eq2st($RAr,$DECr,$RA1r,$DEC1r); #tangent plane
		$count++;
		$RA1r *=180/3.1415926535; #go to degrees
		$DEC1r *=180/3.1415926535; #go to degrees
		print F "$id1 $x $y $mag $RA1r $DEC1r\n"; #write to .xieta
	    }
	}
    }
    close (F);
    close (S);

    if (-e "$img.fits.wcsmatch") {
	$wc=`wc -l $img.fits.wcsmatch`;
	if ($wc < 4 || $force) { 
	    unlink "$img.fits.wcsmatch";
	}
    }
    
    if (! -e "$img.fits.wcsmatch" || $force) {
	print "$starmatches $maxatonce\n";
	if ($scale < 1) { $scaletmp=1;}
	else {$scaletmp=$scale;}
	$amiss*=$scaletmp;
	$cbmiss*=$scaletmp;
	
	
	
	
	$scale="-scale $scale";
	
	$tol="-tol $tol";
	
	$starmatches="-ntri $starmatches";
	
	$maxatonce="-nmatch $maxatonce";
	
    #match using starmatch
	print "$brianwcs/bin/starmatch_new $stars1 $stars2 $scale $tol $starmatches $maxatonce $roterr -outids\n";
	$xform = `$brianwcs/bin/starmatch_new $stars1 $stars2 $scale $tol $starmatches $maxatonce $roterr -outids`;
	print"->$xform<-\n";
	@xf = split( / /, $xform );
	$aa=$xf[0];
	$bb=$xf[1];
	$cc=$xf[2];
	$dd=$xf[3];
	$ee=$xf[4];
	$ff=$xf[5];
	$_=$xf[5];
	($xf[5])=/(\S+)/;
	
	if ( ( ( $xf[ 0 ] == 0 ) && ( $xf[ 1 ] == 0 ) && ( $xf[ 2 ] == 0 ) &&
	       ( $xf[ 0 ] == 0 ) && ( $xf[ 4 ] == 0 ) && ( $xf[ 5 ] == 0 ) ) ) {
	    print "Transform is invalid! \n";
	    next;
	#die;
	} 
	else {
	    print( "\nx = $xf[ 0 ]\*x0 \+ $xf[ 1 ]\*y0 \+ $xf[ 2 ]\n" );
	    print( "y = $xf[ 3 ]\*x0 \+ $xf[ 4 ]\*y0 \+ $xf[ 5 ]\n" );
	}
	
	
	
#match stars and find all that are within 2 pixels
	open( MCH, ">$matchlist" );
	open( S2, "<${stars2}" );
	$i=0;
	
	print "$brianwcs/bin/finalmatch $stars1 $stars2 $matchlist $xf[0] $xf[1] $xf[2] $xf[3] $xf[4] $xf[5] $OKmatch 1000 1 1\n";
	system ("$brianwcs/bin/finalmatch $stars1 $stars2 $matchlist $xf[0] $xf[1] $xf[2] $xf[3] $xf[4] $xf[5] $OKmatch 1000 1 1");
	printf("Final Matching of stars done\n");
    }
#
#figure out appropriate IRAF order if it is set to 0.
#
    
    $_=`wc $matchlist`;
    ($wc)=/(\S+)/;
    if ($iraforder != 0) {
	if ($wc < 5 && $iraforder > 2) { $iraforder = 2;}
	elsif ($wc < 12 && $iraforder > 3) {$iraforder=3;}
	elsif ($wc < 25 && $iraforder > 4) {$iraforder=4;}
    }
    else {
	$_=`wc $matchlist`;
	($wc)=/(\S+)/;
	if ($wc < 5) { $iraforder = 2;}
	elsif ($wc < 12) {$iraforder=3;}
	elsif ($wc < 100) {$iraforder=4;}
	else {$iraforder=5;}
    }
    
    
# derive WCS  and write into header
#use IRAF imcoords.ccmap to fit and write FITS WCS header info
#print "updateWCS.pl $img.fits -o $iraforder $interac\n";
#system ("updateWCS.pl $img.fits -o $iraforder $interac");
    
    unless ($donotapply) {
   
    #clean WCS up and put higher order fit onto it
	system("$brianwcs/bin/updateWCS.pl $img -o $iraforder -clean $clean -sig $cleansig $interac");
    }
    
#
# compute zps if so desired..
#
    if ($zpR || $zpB) {
	open (ZP, ">$img.ZP");
	if ($zpB) {
	    #make list of appropriate stars
	    open (C, "<$img.fits.catalog");
	    open (F, ">tmp.RA");
	    $i=0;
	    while(<C>) {
		$i++;
		($ra,$dec,$magB[$i])=/^\s*\S+\s+(\S+)\s+(\S+)\s+(\S+)/;
		print F "$ra $dec\n";
	    }
	    close (F);
	    close (C);
	    # find these stars' xy position in image
	    open(S, "sky2xy $img.fits \@tmp.RA |");
	    open (F, ">tmp.coo");
	    $j=0;
	    while (<S>) {
		($x,$y)=/\-\>\s*(\S+)\s+(\S+)/;
		if ($x > 0 && $x < 50000 && $y > 0 && $y < 50000) {
		    print F "1 $x $y\n";
		}
		else {
		    print F "1 0 0\n";
		}
		
	    }
	    close (F);
	    close (S);
	    
	    # get stars' mags
	    open (G, "getfluxbk $img.fits tmp.coo -aper 15 -sat 50000 |"); 
	    $j=0;
	    $k=0;
	    while (<G>) {
		$j++;
		($c,$back)=/\s*\S+\s+\S+\s+(\S+)\s+(\S+)/;
		if ($c > 0) {
		    $mag[$k++]=$magB[$j]+2.5*log($c)/2.30259;
		}
	    }
	    close(G);
	    ($ZPBlue)=median(@mag);
	    
	    printf (ZP "ZPBlue:%9.2f N:%d\n",$ZPBlue,$#mag);
	    
	}
	if ($zpR) {
	    #make list of appropriate stars
	    open (C, "<$img.fits.catalog");
	    open (F, ">tmp.RA");
	    $i=0;
	    while(<C>) {
		$i++;
		($ra,$dec,$dummy,$magR[$i])=/\s*\S+\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/;
		if ($magR[$i] == 0) { 
		    $i--;
		    next;
		}
		else {
		    print F "$ra $dec\n";
		}
	    }
	    close (F);
	    close(C);
	    # find these stars' xy position in image
	    print "/sky2xy $img.fits \@tmp.RA |\n";
	    open(S, "/sky2xy $img.fits \@tmp.RA |");
	    open (F, ">tmp.coo");
	    $j=0;
	    while (<S>) {
		($x,$y)=/\-\>\s*(\S+)\s+(\S+)/;
		if ($x > 0 && $x < 50000 && $y > 0 && $y < 50000) {
		    print F "1 $x $y\n";
		    for ($i=1;$i<=$Nstars;$i++) { #get FWHM of real stars
			if ((abs($x-$XX[$i])+abs($y-$YY[$i])) < 10) { 
			    $seeinggood[$j++]=$fwhm[$i];
			    $i=$Nstars+1;
			} 
		    }
		}
		else {
		    print F "1 0 0\n";
		}
	    }
	    close (F);
	    close (S);
	    
	    # get stars' mags
	    open (G, "getfluxbk $img.fits tmp.coo -aper 15 -sat 40000 |"); 
	    $j=0;
	    $k=0;
	    while (<G>) {
		$j++;
		($c,$back)=/\s*\S+\s+\S+\s+(\S+)\s+(\S+)/;
		if ($c > 0) {
		    $magr[$k]=$magR[$j]+2.5*log($c)/2.30259;
		    $k++;
		}
		
	    }
	    close(G);
	    ($ZPRed)=median(@magr);
	    print "$#seeinggood";
	    ($seeing)=median(@seeinggood);
	    
	    printf (ZP "ZPRed: %9.2f (%6.2f) (%6.2f) N:%d %d\n",$ZPRed,$seeing,1.,$#magr,$#seeinggood);
	}
	close (ZP);
    }
}    
    
    sub eq2st { # Compute std coord offsets (arcsec) given RA and Dec and plate centre in rad.
	local( $plate_centre_ra, $plate_centre_dec, $obj_ra, $obj_dec ) = @_;
	local $ARCSEC_PER_RADIAN = 206264.8062470964;
	
	local $div = ( sin( $obj_dec ) * sin( $plate_centre_dec ) +
		       cos( $obj_dec ) * cos( $plate_centre_dec ) *
		       cos( $obj_ra - $plate_centre_ra ) );
	
	# Compute standard coords and convert to arcsec
	local $xi_obj = cos( $obj_dec ) * sin( $obj_ra - $plate_centre_ra ) *
	    $ARCSEC_PER_RADIAN / $div;
	
	local $eta_obj = ( sin( $obj_dec ) * cos( $plate_centre_dec ) -
			   cos( $obj_dec ) * sin( $plate_centre_dec ) *
			   cos( $obj_ra - $plate_centre_ra ) ) *
			   $ARCSEC_PER_RADIAN / $div;
	
	( $xi_obj, $eta_obj );
    }
    
    sub st2eq { # Compute RA and Dec given std coord offsets (arcsec) and plate centre
	local( $plate_centre_ra, $plate_centre_dec, $xi, $eta ) = @_;
	local $ARCSEC_PER_RADIAN = 206264.8062470964;
	local $TWOPI = 2.0 * 3.14159265358979323846;
	
	local $object_xi = $xi / $ARCSEC_PER_RADIAN;
	local $object_eta = $eta / $ARCSEC_PER_RADIAN;
	
	# Convert to RA and Dec
	local $numerator = $object_xi;
	
	local $denominator = cos( $plate_centre_dec ) - 
	    $object_eta * sin( $plate_centre_dec );
	local $ra = atan2( $numerator, $denominator ) + $plate_centre_ra;
	if ( $ra < 0.0 ) { $ra = $ra + $TWOPI; }
	
	$numerator = cos( $ra - $plate_centre_ra ) *
	    ( cos( $plate_centre_dec ) * $object_eta + 
	      sin( $plate_centre_dec ) );
	
	$denominator = cos( $plate_centre_dec ) -
	    $object_eta * sin( $plate_centre_dec );
	$dec = atan2( $numerator, $denominator );
	
	( $ra, $dec );
    }
    
    
    sub separation { # Return angular separation, degrees.
	local ( $ra1, $dec1, $ra2, $dec2 ) = @_;
	local $hours = 3.819718634205; # Hours in a radian.
	local $degrees = 57.2957795130823; # Degrees in a radian.
	local $x1 = cos( $ra1 / $hours ) * cos( $dec1 / $degrees );
	local $y1 = sin( $ra1 / $hours ) * cos( $dec1 / $degrees );
	local $z1 = sin( $dec1 / $degrees );
	local $x2 = cos( $ra2 / $hours ) * cos( $dec2 / $degrees );
	local $y2 = sin( $ra2 / $hours ) * cos( $dec2 / $degrees );
	local $z2 = sin( $dec2 / $degrees );
	local $theta = acos( $x1*$x2 + $y1*$y2 + $z1*$z2 ) * $degrees;
	$theta;
    }   
    
    sub acos {
	local ( $inp ) = @_;
	while( $inp > 1.0 ) { $inp -= 1.0; }
	while( $inp < -1.0 ) { $inp += 1.0; }
	local $outp = atan2( sqrt( 1.0 - ( $inp * $inp ) ), $inp );
    }
    
    sub radecstring2rad {
	local( $inra, $indec ) = @_;
	
	local $ARCSEC_PER_RADIAN = 206264.8062470964;
	local( $raH, $raM, $raS ) = split( /:/, $inra );
	local( $decD, $decM, $decS ) = split( /:/, $indec );
	
	if ( $decD =~ /(-\d+)/ ) {
	    $dec = -(abs($decD) + $decM / 60.0 + $decS / 3600.0);
	} else {
	    $dec = $decD + $decM / 60.0 + $decS / 3600.0;
	}
	$ra = $raH + $raM / 60.0 + $raS / 3600.0;
	
	$decrad = $dec * 3600.0 / $ARCSEC_PER_RADIAN;
	$rarad = $ra * 15.0 * 3600.0 / $ARCSEC_PER_RADIAN;
	
	( $rarad, $decrad );
    }
    
    sub radecrad2string {
	local( $rarad, $decrad ) = @_;
	
	local $ARCSEC_PER_RADIAN = 206264.8062470964;
	local $ra = $rarad * $ARCSEC_PER_RADIAN / 3600.0 / 15.0;
	local $dec = $decrad * $ARCSEC_PER_RADIAN / 3600.0;
	
	local $raH = int( $ra );
	local $raM = int( ( $ra - int( $ra ) ) * 60.0 );
	local $raS = ( ( $ra - int( $ra ) ) - $raM / 60.0 ) * 3600.0;
	$raS = int( $raS * 100.0 + 0.5 ) / 100.0;
	
	local $decD = int( $dec );
	local $decM = int( abs( $dec - int( $dec ) ) * 60.0 );
	local $decS = ( abs( $dec - int( $dec ) ) - $decM / 60.0 ) * 3600.0;
	$decS = int( $decS * 100.0 + 0.5 ) / 100.0;
	
	
	if ($raS >= 60) {
	    $raS-=60;
	    $raM+=1;
	}
	if ($raM >= 60) {
	    $raM-=60;
	    $raH+=1;
	}
	if ($raH >=24) {
	    $raH-=24;
	}
	
	if ($decS >= 60) {
	    $decS-=60;
	    $decM+=1;
	}
	if ($decM >= 60) {
	    $decM-=60;
	    $decD+=1;
	}
	
	
	
	print "$raH $raM $raS $decD $decM $decS\n";
	
	
	if ( $raH < 10 ) { $raH = '0' . $raH; }
	if ( $raM < 10 ) { $raM = '0' . $raM; }
	if ( $raS < 10 ) { $raS = '0' . $raS; }
	
	$_=$decD;
	if (! /-0/) {
	    if ( ( $decD < 10 ) && ( $decD >= 0 ) ) { $decD = '0' . $decD; }
	    elsif ( ( $decD > -10 ) && ( $decD < 0 ) ) { $decD = '-0' . abs( $decD ); }
	}
	else {
	    $decD = '-00';
	}
	if ( $decM < 10 ) { $decM = '0' . $decM; }
	if ( $decS < 10 ) { $decS = '0' . $decS; }
	
	local $rastr = join( ':', $raH, $raM, $raS );
	local $decstr = join( ':', $decD, $decM, $decS );
	
	print "$rastr $decstr\n";
	( $rastr, $decstr );
    }
    
    sub getoffsets {
	local ( $racen, $deccen, $inra, $indec ) = @_;
	print "Getting offsets between $racen $deccen and $inra $indec\n";
	local $ARCSEC_PER_RADIAN = 206264.8062470964;
	
	local ( $racenrad, $deccenrad ) = radecstring2rad( $racen, $deccen );
	local ( $rarad, $decrad ) = radecstring2rad( $inra, $indec );
	
	local $racenas = $racenrad * $ARCSEC_PER_RADIAN;
	local $deccenas = $deccenrad * $ARCSEC_PER_RADIAN;
	local $raas = $rarad * $ARCSEC_PER_RADIAN;
	local $decas = $decrad * $ARCSEC_PER_RADIAN;
	
	local $offNS = ( $deccenas - $decas );
	if ( $offNS >= 0 ) { $offNS = int( abs( $offNS ) * 1000.0 + 0.5 ) / 1000.0 . "\" N"; }
	elsif ( $offNS < 0 ) { $offNS = int( abs( $offNS ) * 1000.0 + 0.5 ) / 1000.0 . "\" S"; }
	
	local $offEW = ( $racenas - $raas ) * cos( ( $deccenrad + $decrad ) / 2 );
	if ( $offEW >= 0 ) { $offEW = int( abs( $offEW ) * 1000.0 + 0.5 ) / 1000.0 . "\" E"; }
	elsif ( $offEW < 0 ) { $offEW = int( abs( $offEW ) * 1000.0 + 0.5 ) / 1000.0 . "\" W"; }
	
	( $offEW, $offNS );
    }
    
    
    sub radec2tex {
	local( $inra, $indec ) = @_;
	
	local( $raH, $raM, $raS ) = split( /:/, $inra );
	local( $decD, $decM, $decS ) = split( /:/, $indec );
	
	$raS = int( $raS * 1000.0 + 0.5 ) / 1000.0;
	$decS = int( $decS * 1000.0 + 0.5 ) / 1000.0;
	
	local $outra = "\$\\alpha= $raH\^h\~$raM\^m\~$raS\^s\$";
	local $outdec = "\$\\delta= $decD^\\circ\~$decM\^\{\\prime\}\~$decS\^\{\\prime\\prime\}\$";
	( $outra, $outdec );
    }
    
    
    
    

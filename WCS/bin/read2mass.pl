#!/usr/bin/env perl

if ( $#ARGV !=3 ) { die "Usage: vizquery.pl catalog outname RA DEC box (arcsec)\n";
 }

($fileout,$RA,$DEC,$size)=@ARGV;

$DEG2RAD=3.1415926535/180.;
$maxJ=$maxH=$maxK=0.15; #
open (OUT, ">$fileout") || die "cannot open $fileout for writing\n";

require "skysubroutines.pl";

($ra,$dec)=radecstring2rad($RA,$DEC);
$racent=$ra*180/3.1415926535;
$deccent=$dec*180/3.1415926535;
$size=$size/60;

# RS 2012/04/17:  fixing sign-related crash
my $sgn = (($deccent > 0) ? "+" : "");

# RS 2012/05/10:  added -get flag to vizquery after advice from CDS help desk
print "vizquery -get -mime=tsv -source=2MASS -c=$racent$sgn$deccent -c.eq=j2000 -c.rm=$size\n";
open (CAT,"vizquery -get -mime=tsv -source=2MASS -c=$racent$sgn$deccent -c.eq=j2000 -c.rm=$size | sort -nk 4 |") || die "cannot open vizquery\n";

$start=0;
$i=1;
while(<CAT>) {
    if (/^#/) { 
	next;
    }
    if ($start==1) {
	if (/^\s*$/) {last;} #quit with a blank line
	($rat[$i],$dect[$i],$id,$J[$i],$Je[$i],$H[$i],$He[$i],$K[$i],$Ke[$i])=split(/\t/);
	if (($Je[$i] < $maxJ && $Je[$i] > 0) || ($He[$i] < $maxH && $He[$i] > 0)  || ($Ke[$i] < $maxK && $Ke[$i] > 0) ) { #has to have SNR=1 in at least one band
	    ($ra,$dec)=radecrad2string($rat[$i]*$DEG2RAD,$dect[$i]*$DEG2RAD);
	    printf(OUT "%5d %12s %12s %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n",$i,$ra,$dec,$J[$i],$Je[$i],$H[$i],$He[$i],$K[$i],$Ke[$i]);
	    $i++;
	}
    }
    if (/--------/) { $start=1;}
}
close(OUT);
    


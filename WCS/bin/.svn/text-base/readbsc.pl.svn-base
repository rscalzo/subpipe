#!/usr/bin/env perl

#
#reads modified version of the BrightStarCatalog
#
if ( $#ARGV !=3 ) { die "Usage: readbsc.pl outname RA DEC box (degrees)\nReads Bright Star Catalog and output catalog\n";
 }

($fileout,$RA,$DEC,$size)=@ARGV;

open (OUT, ">$fileout") || die "cannot open $fileout for writing\n";

$m1=0;
$m2=7;
$epoch=2000.0;

for ($i = 4; $i <= $#ARGV; $i++) {
    die "Your ${i}th parameter $ARGV[$i] not understood\n";
}


require "skysubroutines.pl";


($ra,$dec)=radecstring2rad($RA,$DEC);
$racent=$ra;
$deccent=$dec;


open (CAT,"/data/brian/BrightStar.dat\n") || die "cannot open Bright Star Catalog";
$id=0;
while(<CAT>) {
    @data=split(/\s+/);
    $id++;
    $ras=$data[0];
    $decs=$data[1];
    $mag=$data[2];
    ($myrar,$mydecr)=radecstring2rad($ras,$decs);
    ($a,$b)=eq2st($racent,$deccent,$myrar,$mydecr);
    $dist=sqrt($a*$a+$b*$b)/3600;
    if ($dist < $size && abs($racent-$myrar)< 3.14159/2) {
	printf (OUT "%6d %15s %15s %6.2f %.2f\n",$id,$ras,$decs,$mag,$dist);
    }
}
close (OUT);

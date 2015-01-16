#!/usr/bin/env perl

if ( $#ARGV !=3 ) { die "Usage: readucac1.pl outname RA DEC box (arcsec)\n";
 }

($fileout,$RA,$DEC,$size)=@ARGV;

open (OUT, ">$fileout") || die "cannot open $fileout for writing\n";

for ($i = 4; $i <= $#ARGV; $i++) {
    if ($ARGV[$i] eq '-s') {
	$i++;
	$seeing=$ARGV[$i];
    } elsif ($ARGV[$i] eq '-stars') {
	$i++;
	$needthismanystars=$ARGV[$i];
    }
}

require "skysubroutines.pl";


($ra,$dec)=radecstring2rad($RA,$DEC);
$racent=$ra;
$deccent=$dec;

$low=-1*$size;
$hi=$size;

($ra2,$dec2)=st2eq($racent,$deccent,$hi,$hi);
print "$deccent $dec1 $dec2 $hi $low \n";
($ra1,$dec1)=st2eq($racent,$deccent,$low,$low);


$ra1*=12./3.1415926535;
$ra2*=12./3.1415926535;
$dec1*=180/3.1415926535;
$dec2*=180/3.1415926535;

print "$ra1 $ra2 $dec1 $dec2\n";

open (CAT, "| /data/u1/unix/catread");
print CAT "b\n";
print CAT "$ra1 $ra2\n";
print CAT "$dec1 $dec2\n";
print CAT "x\n";

close (CAT);


open (CAT,"sort -nk 3 fort.9 |");
<CAT>;
<CAT>;
<CAT>;
<CAT>;
while(<CAT>) {
    ($sra,$spd,$mag)=/\s*(\S+)\s+(\S+)\s+(\S+)/;
    $s++;
    $ra=$sra/206264.806/1000.;  #in RA RAD
    $dec=($spd-324000000)/206264.806/1000.; #DEC RAD
    ($xi,$eta)=eq2st($ra,$dec,$racent,$deccent);
    ($ra,$dec)=radecrad2string($ra,$dec);
    $dist=sqrt($xi*$xi+$eta*$eta);
    printf (OUT "%6d %15s %15s %6.2f %6.2f %10.3f\n",$s,$ra,$dec,$mag/100.,0,$dist);
}
close (OUT);


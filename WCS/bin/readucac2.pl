#!/usr/bin/env perl

if ( $#ARGV !=3 ) { die "Usage: readucac1.pl outname RA DEC box (arcsec)\n";
 }

($fileout,$RA,$DEC,$size)=@ARGV;

open (OUT, ">$fileout") || die "cannot open $fileout for writing\n";

$m1=5;
$m2=18;
$epoch=2000.0;

for ($i = 4; $i <= $#ARGV; $i++) {
    if ($ARGV[$i] eq '-s') {
	$i++;
	$seeing=$ARGV[$i];
    } elsif ($ARGV[$i] eq '-m1') {
	$i++;
	$m1=$ARGV[$i];
    }
    elsif ($ARGV[$i] eq '-m2') {
	$i++;
	$m2=$ARGV[$i];
    }
    elsif ($ARGV[$i] eq '-epoch') {
	$i++;
	$epoch=$ARGV[$i];
    }
    else {
	die "Your ${i}th parameter $ARGV[$i] not understood\n";
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
if (-e "tmp.ucac2") { unlink "tmp.ucac2";}
open (CAT, "| u2access");
print CAT "/home/brian/CATALOGS/UCAC2/\n"; #path for original UCAC2 files = 
print CAT "3\n"; #fmt = 3 = updated RA, Dec (hms) photom..
print CAT "tmp.ucac2\n"; #file name output table or "s" = screen dump
print CAT "$m1 $m2\n"; #select magnitude range (UCAC mags, 5 - 18 = all):
print CAT "$epoch\n"; #update positions for proper motions to new epoch
print CAT "3\n"; #no sort = default use mag instead
print CAT "\n"; #-- no byte flip required
print CAT "\n"; #-- no byte flip required
print CAT "r\n"; #== new box: r=range, c=center, q=quit r
print CAT "$ra1 $ra2\n"; # RA
print CAT "$dec1 $dec2\n"; #DEC
print CAT "q\n";  #== new box: r=range, c=center, q=quit r
close (CAT);

open (CAT,"tmp.ucac2\n") || die "cannot open tmp.ucac2";
<CAT>;
<CAT>;
<CAT>;
<CAT>;
<CAT>;
<CAT>;
<CAT>;
<CAT>;
<CAT>;
$id=0;
while(<CAT>) {
    @data=split(/\s+/);
    $id++;
    $ra="$data[0]:$data[1]:$data[2]";
    $dec="$data[3]:$data[4]:$data[5]";
    $rmag=$data[13];
    $jmag=$data[14];
    $hmag=$data[15];
    $kmag=$data[16];
    printf (OUT "%6d %15s %15s %6.2f %6.2f %6.3f %6.3f\n",$id,$ra,$dec,$rmag,$jmag,$hmag,$kmag);
#    printf ("%6d %15s %15s %6.2f %6.2f %6.3f %6.3f\n",$id,$ra,$dec,$rmag,$jmag,$hmag,$kmag);
}
close (OUT);
if (-e "tmp.ucac2") { unlink "tmp.ucac2";}


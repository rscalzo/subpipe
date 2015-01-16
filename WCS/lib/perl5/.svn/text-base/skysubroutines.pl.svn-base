   
sub tworeg { 

#x2=a*x1+b*y1+c

local(*x1,*y1,*x2,*yes)=@_;

local $i;
local $sum=0;
local $sumx=0;
local $sumy=0.;
local $sumxy=0.;
local $sumu=0.;
local $sumx2=0.;
local $sumy2=0.;
local $sumux=0.;
local $sumuy=0.;
local $delta=0;

for ($i=0; $i<=$#x1;$i++) {
    if ($yes[$i]==0) {
	$sum ++;
	$sumx += $x1[$i];
	$sumy += $y1[$i];
	$sumu += $x2[$i];
	$sumx2 += $x1[$i]*$x1[$i];
	$sumxy += $x1[$i]*$y1[$i];
	$sumy2 += $y1[$i]*$y1[$i];
	$sumux += $x1[$i]*$x2[$i];
	$sumuy += $y1[$i]*$x2[$i];
    }
}
$delta = $sum*($sumx2*$sumy2 - $sumxy*$sumxy) + $sumx*($sumxy*$sumy - $sumx*$sumy2) + $sumy*($sumx*$sumxy-$sumx2*$sumy);

local $a = ($sumu*($sumx2*$sumy2 - $sumxy*$sumxy) + $sumx*($sumxy*$sumuy - $sumux*$sumy2) + $sumy*($sumux*$sumxy - $sumx2*$sumuy))/$delta;
local $b = ($sum*($sumux*$sumy2 - $sumxy*$sumuy) + $sumu*($sumxy*$sumy - $sumx*$sumy2) + $sumy*($sumuy*$sumx - $sumy*$sumux))/$delta;
local $c = ($sum*($sumx2*$sumuy - $sumxy*$sumux) + $sumx*($sumux*$sumy - $sumx*$sumuy) + $sumu*($sumx*$sumxy - $sumx2*$sumy))/$delta;

    return ($a,$b,$c);
}

sub onereg { 
    
#x=a*y+b;
    
local(*x,*y,*yes)=@_;

local $i;
local $sum=0;
local $sumx=0;
local $sumy=0.;
local $sumxy=0.;
local $sumu=0.;
local $sumx2=0.;
local $sumy2=0.;
local $sumux=0.;
local $sumuy=0.;
local $delta=0;

for ($i=0; $i<=$#x;$i++) {
    if ($yes[$i]==0) {
	$sum ++;
	$sumx += $x[$i];
	$sumy += $y[$i];
	$sumx2 += $x[$i]*$x[$i];
	$sumxy += $x[$i]*$y[$i];
	$sumy2 += $y[$i]*$y[$i];
    }
}
$b = ( $sumy * $sumx2 - $sumx * $sumxy ) / ($sum * $sumx2 - $sumx * $sumx );
$s = - ( $sumy * $sumx -$sum * $sumxy ) / ($sum * $sumx2 - $sumx * $sumx );
return ($b,$s);
}
sub getmedsig {
   
    local @median = @_;
    local $mid=int($#median/2)+1;
    local $lowsig=int($#median/6)+1;
    local $hisig=int($#median/6*5)+1;
    local (@sorted);
    @sorted = sort sortbynumber @median;
    return ($sorted[$mid],$sorted[$lowsig],$sorted[$hisig]);
}  


sub median1 {
    local @median= @_;
    local $mid=int($#median/2)+1;
    local (@sorted);
    @sorted = sort sortbynumber @median;
    
    if ($#median/2 eq int($#median/2)) {
	return ($sorted[$mid]);
    }
    else {
	return (($sorted[$mid]+$sorted[$mid-1])/2);
    }
}  



sub medianfraction {
    local (*median,$fract)= @_;
    local $mid=int($#median*$fract)+1;
    local (@sorted);
    @sorted = sort sortbynumber @median;
    
#    print "Mid $mid Median $#median Fract $fract\n";
    return ($sorted[$mid]);
}  

sub swap {
    local ($a,$b)=@_;
    local $tmp;

    $tmp=$a;
    $a=$b;
    $b=$tmp;
    return ($a,$b);
}

sub minmax {
(*x,*yes)=@_;

local $i;
local $min=1e30,$max=-1e30;
for ($i=0;$i<=$#x;$i++) {
#    print "minmax: $i $yes[$i] $x[$i] $min $max\n";
    unless ($yes[$i] == 1) {
	if ($x[$i] > $max) {$max=$x[$i];}
	if ($x[$i] < $min) {$min=$x[$i];}
    }
}
return($min,$max);
}



sub max {
    local @median= @_;
    local $mid=$#median;
    local (@sorted);
    @sorted = sort sortbynumber @median;
    return ($sorted[$mid]);
}  

sub min {
    local @median= @_;
    local $mid=0;
    local (@sorted);
    @sorted = sort sortbynumber @median;
    return ($sorted[$mid]);
}  


    
sub sortbynumber { $a<=>$b;}    

sub median {
    local @median = @_;
    local $i;

#    print "median list @median\n";

    local $mid=int($#median/2)+1;

    open (F, ">tmp");
    for ($i=0;$i <=$#median;$i++) {
	if ($median[$i] ne "") {
	    print F "$median[$i]\n";
	}
    }
    close (F);
    $_=`sort -n tmp | awk \'\{s++;if (s == $mid) print\}\'`;
    
    ($med)=/(\S+)/;
    unlink "tmp" if -e "tmp";
    ($med);
}
	
sub stats {
    local (*data,*YES)=@_;

    
    local $n=$#data+1;
    local $sum=0;
    local $sumx=0;
    local $sumxx=0;
    local $i;
    local $ave=0.;


    for ($i=0;$i<=$#data;$i++) {
	if ($YES[$i] == 0) {
	    $sum+=$data[$i];
	}
	else {
	    $n--;
	}
    }
    if ($n >0) {    $ave=$sum/$n;}
    else {$ave=0;}
    
    
    $sum=0;
    for ($i=0;$i<=$#data;$i++) {
	if ($YES[$i] == 0) {
	    $sum+=($data[$i]-$ave)*($data[$i]-$ave);
	}
    }
    if ($n > 0) {
	$std=sqrt($sum/$n);
    }
    else {
	$std=0;
    }
    return($ave,$std);

}
    
	




sub acos {

    local ( $inp ) = @_;

    while( $inp > 1.0 ) { $inp -= 1.0; }

    while( $inp < -1.0 ) { $inp += 1.0; }

    local $outp = atan2( sqrt( 1.0 - ( $inp * $inp ) ), $inp );

}



#converts RADEC to Radians



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



#converts Radians to RADEC string



sub radecrad2string {

    local( $rarad, $decrad ) = @_;



    local $ARCSEC_PER_RADIAN = 206264.8062470964;

    local $ra = $rarad * $ARCSEC_PER_RADIAN / 3600.0 / 15.0;

    local $dec = $decrad * $ARCSEC_PER_RADIAN / 3600.0;


    local $raH = int( $ra );

    local $raM = int( ( $ra - int( $ra ) ) * 60.0 );

    local $raS = ( ( $ra - int( $ra ) ) - $raM / 60.0 ) * 3600.0;

    $raS = int( $raS * 100.0 + 0.5 ) / 100.0;

    if ($raS >= 60) {
	$raS = 0;
	$raM +=1;
    }
    if ($raM >= 60) {
	$raM=0;
	$raH+=1;
    }
    if ($raH >= 24) {
	$raH=0;
    }
    
    if ($dec < 0) {
	$dec=abs($dec);
	$sign="-";
    }
    else {
	$sign="+";
    }
    local $decD = int( $dec );

    local $decM = int( abs( $dec - int( $dec ) ) * 60.0 );

    local $decS = ( abs( $dec - int( $dec ) ) - $decM / 60.0 ) * 3600.0;

    $decS = int( $decS * 100.0 + 0.5 ) / 100.0;

    if ($decS >= 60) {
	$decS = 0;
	$decM +=1;
    }
    if ($decM >= 60) {
	$decM=0;
	$decD+=1; #SK 130706 bug fix: TMJ found this was "$decH"
    }
    #replaced below on advice of TMJ 20Jul07 code change..
    return ( sprintf( "%02u:%02u:%06.3f", $raH, $raM, $raS ), $sign . sprintf( "%02u:%02u:%05.2f", $decD, $decM, $decS ) ) ;
}







sub eq2st { # Compute std coord offsets (arcsec) given RA and Dec and plate centre in rad.

    my( $plate_centre_ra, $plate_centre_dec, $obj_ra, $obj_dec ) = @_;

    my $ARCSEC_PER_RADIAN = 206264.8062470964;



    my $div = ( sin( $obj_dec ) * sin( $plate_centre_dec ) +

		  cos( $obj_dec ) * cos( $plate_centre_dec ) *

		  cos( $obj_ra - $plate_centre_ra ) );

 

    # Compute standard coords and convert to arcsec

    my $xi_obj = cos( $obj_dec ) * sin( $obj_ra - $plate_centre_ra ) *

	$ARCSEC_PER_RADIAN / $div;

 

    my $eta_obj = ( sin( $obj_dec ) * cos( $plate_centre_dec ) -

		      cos( $obj_dec ) * sin( $plate_centre_dec ) *

		      cos( $obj_ra - $plate_centre_ra ) ) *

			  $ARCSEC_PER_RADIAN / $div;



    ( $xi_obj, $eta_obj );

}



sub st2eq { # Compute RA and Dec given std coord offsets (arcsec) and plate centre

    my( $plate_centre_ra, $plate_centre_dec, $xi, $eta ) = @_;

    my $ARCSEC_PER_RADIAN = 206264.8062470964;

    my $TWOPI = 2.0 * 3.14159265358979323846;

    

    my $object_xi = $xi / $ARCSEC_PER_RADIAN;

    my $object_eta = $eta / $ARCSEC_PER_RADIAN;

 

    # Convert to RA and Dec

    my $numerator = $object_xi;

 

    my $denominator = cos( $plate_centre_dec ) - 

	$object_eta * sin( $plate_centre_dec );

    my  $ra = atan2( $numerator, $denominator ) + $plate_centre_ra;

    if ( $ra < 0.0 ) { $ra = $ra + $TWOPI; }

 

    $numerator = cos( $ra - $plate_centre_ra ) *

	( cos( $plate_centre_dec ) * $object_eta + 

	 sin( $plate_centre_dec ) );

    

    $denominator = cos( $plate_centre_dec ) -

	$object_eta * sin( $plate_centre_dec );

    $dec = atan2( $numerator, $denominator );



    ( $ra, $dec );

}

sub linint { #linerar interpolate 
local(*x,*y,$newX,)=@_;
local $i,$true=0;
$i=0;
while ($newX > $x[$i] && $i<$#x) {
    $i++;
}

local $x1,$x2,$y1,$y2,$newY;

if ($newX > $x[$#x] || $newX < $x[0])  {
    $newY=0;
}
else {

    $x1=$x[$i-1];
    $x2=$x[$i];
    $y1=$y[$i-1];
    $y2=$y[$i];

    
    $newY=($newX-$x1)*($y2-$y1)/($x2-$x1)+$y1;
}

return($newY);
}

sub date2jd { # 
    local($m,$day,$year) = @_;
    local $A,$B,$C,$D,$jd,$y;
    $y=$year;
    if ($m <=2) {
	$y=$year-1;
	$m=$m+12;
    }
    $B=0;
    if ($y >= 1582 ) {
	$A=int($y/100.);
	$B=2-$A+int($A/4.);
    }
 $C=int(365.25*$y);
    $D=int(30.6001*($m+1));
    $jd=$B+$C+$D+$day+1720994.5;
    return ($jd);
}

sub jd2date {
    local ($jd) = @_;
    local $A,$B,$C,$D,$E,$F,$G,$I,$day,$month,$year;
    $jd=$jd+.5;
    $I=int($jd);
    $F=$jd-$I;
    $B=$I;
    if ($I>2299160) {
	$A=int(($I-1867216.25)/(36524.25));
	$B=$I+1+$A-int($A/4);
    }

    $C=$B+1524.;
    $D=int(($C-122.1)/365.25);
    $E=int(365.25*$D);
    $G=int(($C-$E)/30.6001);
    $day=$C-$E+$F-int(30.6001*$G);
    $month=$G-1;
    if ($month > 12.5) {
	$month= $month-12.;
    }
    $year=$D-4715;
    if ($month>2.5) {
	$year=$D-4716;
    }
    return ($month,$day,$year);
}


1; #return true;

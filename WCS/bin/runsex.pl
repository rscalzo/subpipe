#!/usr/bin/env perl

if ( $#ARGV <1 ) { die "Usage: runsex.pl image sigma -sat (saturated pixels 30000) -zp (zeropoint) -nobad -mask maskname\n ";
 }
($img,$sigma)=@ARGV;

# RS 2011/11/02:  No help for it, sextractor needs to know where these files
# are stored.  Probably the least evil way is to set a BRIANWCS environment
# variable to point to the base of the tree.
$brianwcs = $ENV{ 'BRIANWCS' } || die "Whoops!  Please setenv BRIANWCS and try again";

$sat=60000;
$zp=25;
$nobad=0;
$weight_image="None";
$weight_type="None";

for ($i = 2; $i <= $#ARGV; $i++) {
    if ($ARGV[$i] eq '-sat') {
	$i++;
	$sat=$ARGV[$i];
    }
    if ($ARGV[$i] eq '-zp') {
	$i++;
	$zp=$ARGV[$i];
    }
    if ($ARGV[$i] eq '-nobad') {
	$nobad=1;
    }
    if ($ARGV[$i] eq '-mask') {
	$i++;
	$weight_image=$ARGV[$i];
	$weight_type="MAP_WEIGHT";
    }
}

# The configuration files must then be in $BRIANWCS/etc/.
(system("cp $brianwcs/etc/{daofind.param,default.conv,default.nnw} .") == 0)
   || die "Can't find sextractor config files in $BRIANWCS/etc";

open (F,">default.sex");
print  F <<EOT;

# Default configuration file for SExtractor V1.2
# EB 18/08/97
# (*) indicates parameters which can be omitted from this config file.

#-------------------------------- Catalog ------------------------------------

CATALOG_NAME	test.cat	# name of the output catalog
CATALOG_TYPE	ASCII_HEAD	# "ASCII_HEAD","ASCII","FITS_1.0" or "FITS_LDAC"

PARAMETERS_NAME	daofind.param	# name of the file containing catalog contents

#------------------------------- Extraction ----------------------------------

DETECT_TYPE	CCD		# "CCD" or "PHOTO" (*)
#DETECT_IMAGE	SAME		# "SAME" or <image filename>
FLAG_IMAGE	flag.fits	# filename for an input FLAG-image
DETECT_MINAREA	5		# minimum number of pixels above threshold
DETECT_THRESH	$sigma		# <sigmas> or <threshold>,<ZP> in mag.arcsec-2
ANALYSIS_THRESH	$sigma		# <sigmas> or <threshold>,<ZP> in mag.arcsec-2

FILTER		Y		# apply filter for detection ("Y" or "N")?
FILTER_NAME	default.conv	# name of the file containing the filter

DEBLEND_NTHRESH	32		# Number of deblending sub-thresholds
DEBLEND_MINCONT	0.005		# Minimum contrast parameter for deblending

CLEAN		Y		# Clean spurious detections? (Y or N)?
CLEAN_PARAM	1.0		# Cleaning efficiency

#BLANK		Y		# Blank detected objects (Y or N)?

#------------------------------ Photometry -----------------------------------

PHOT_APERTURES	5		# MAG_APER aperture diameter(s) in pixels
PHOT_AUTOPARAMS	2.5, 3.5	# MAG_AUTO parameters: <Kron_fact>,<min_radius>

SATUR_LEVEL	$sat		# level (in ADUs) at which arises saturation

MAG_ZEROPOINT	$zp		# magnitude zero-point
MAG_GAMMA	4.0		# gamma of emulsion (for photographic scans)
GAIN		0.0		# detector gain in e-/ADU.
PIXEL_SCALE	1.0		# size of pixel in arcsec (0=use FITS WCS info).

#------------------------- Star/Galaxy Separation ----------------------------

SEEING_FWHM	1.2		# stellar FWHM in arcsec
STARNNW_NAME	default.nnw	# Neural-Network_Weight table filename

#------------------------------ Background -----------------------------------

BACK_SIZE	64		# Background mesh: <size> or <width>,<height>
BACK_FILTERSIZE	3		# Background filter: <size> or <width>,<height>

BACKPHOTO_TYPE	GLOBAL		# can be "GLOBAL" or "LOCAL" (*)
BACKPHOTO_THICK	24		# thickness of the background LOCAL annulus (*)

#------------------------------ Check Image ----------------------------------

CHECKIMAGE_TYPE	NONE		# can be one of "NONE", "BACKGROUND",
				# "MINIBACKGROUND", "-BACKGROUND", "OBJECTS",
				# "-OBJECTS", "SEGMENTATION", "APERTURES",
				# or "FILTERED" (*)
CHECKIMAGE_NAME	check.fits	# Filename for the check-image (*)

#--------------------- Memory (change with caution!) -------------------------

MEMORY_OBJSTACK	2000		# number of objects in stack
MEMORY_PIXSTACK	100000		# number of pixels in stack
MEMORY_BUFSIZE	512		# number of lines in buffer

#---------------- Scanning parameters (change with caution!) -----------------

#SCAN_ISOAPRATIO	0.6		# maximum isoph. to apert ratio allowed (*)

#----------------------------- Miscellaneous ---------------------------------

VERBOSE_TYPE	NORMAL		# can be "QUIET", "NORMAL" or "FULL" (*)

#------------------------------- New Stuff -----------------------------------

WEIGHT_IMAGE    $weight_image
WEIGHT_TYPE     $weight_type
WEIGHT_THRESH   0.0001
# Surprise!!

EOT

    close(F);
 system("sex $img -VERBOSE_TYPE QUIET");

if ($nobad) {
    system("sort -n -k 4 test.cat | awk \'\{ if (\$5 < 1) print\}\' > $img.stars"); 
}
else {system("sort -n -k 4 test.cat > $img.stars");
  }
system("rm -f daofind.param default.conv default.nnw default.sex test.cat");



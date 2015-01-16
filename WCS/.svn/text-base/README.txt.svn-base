------------------------------------------------------------------------------
** README on building Brian's WCS code
** RS 2012/08/28 Updated for Nicolas.  Includes information about central
   shared package build.
** RS 2011/11/18 Pitched cirdr.  Brian only needed three include files and
   three source files, so I've just included these in the distribution.
   This solves several problems with linking various libraries which I kept
   running into every time I tried to build.
** RS 2011/11/02 Original draft
------------------------------------------------------------------------------


I.  The short version

After checking out the code from svn, it should be in the following tree:

   WCS/
      bin/ (includes executable python and perl scripts)
      src/ (includes source code for brian_fitwcs and brian_xi_eta)
      lib/ (includes skysubroutines.pl in subdirectory perl5/)
      etc/ (includes configuration files needed for the code to run)

To set up the code, do the following from within the WCS directory:

   setenv BRIANWCS `pwd`
   setenv PATH $BRIANWCS/bin:$PATH
   setenv EXTROOT /export/maipenrai/skymap/skymap/external-packages/
   cd src
   make install

To use it, from within a directory containing the output of fits64_to_32
(with 32 subdirectories each with a single-CCD extension file of the form,
e.g., 17/Skymapper_805424129_00000_2011-08-25T23:01:01_17.fits), run:

   SM-WCS.py Skymapper_805424129_00000_2011-08-25T23:01:01.fits

If all goes well, at the end all of the extensions should have TNX-format
WCS headers with residuals on the order of 0.1 arcsec.


II. Detailed list of dependencies

Here's a list of the dependencies of Brian's code, including all local
installs of externally-maintained packages.  Before you start, make sure you
have installed all of the necessary external dependencies.

   external packages (install these first if you don't have 'em):
      CFITSIO           http://heasarc.gsfc.nasa.gov/fitsio/
      WCStools          http://tdc-www.harvard.edu/wcstools/
      WCSlib 3.5        http://www.atnf.csiro.au/people/mcalabre/WCS/
      cdsclient-3.6     http://cdsarc.u-strasbg.fr/doc/cdsclient.html
      UCAC2 u2access    http://badc.lamost.org/archives/UCAC2/sw/
      pyRAF             http://www.stsci.edu/resources/software_hardware/pyraf
   WCS.pl -- calls:
      skysubroutines.pl (20 subroutines, self-contained)
      read*.pl (for reading different astronomical catalogs)
      packages from SVN/Code/COMMONSOURCE/:
         starmatch_new (requires starmatch.c, indexx.c, hunt.c, nrutil.c)
         finalmatch (requires finalmatch.c, nrutil.c)
      WCStools binaries:  gethead, sethead, scat, sky2xy, xy2sky
      cdsclient executables:  vizquery, wwwget
      UCAC2 u2access
   brian_fitwcs:  binary written in C by Brian
      (original in /export/miner1/skymap/brian/SOFTWARE/WCS/IoA/cirdr/src/)
      depends on:
         CIRDR (libcirdr.a with included *.c, *.h, see this same path)
         CFITSIO (libcfitsio.a)
         WCSLIB (libwcs.a)
      Brian's relevant code (also in this path):
         brian_fitwcs.c invert_c.c minimrq.c select.c nrutil.c
         brian_xi_eta.c
   updateWCS.py -- calls:
      skysubroutines.pl (20 subroutines, self-contained)
      WCStools binaries:  gethead, sethead, wcshead, delwcs
      pyRAF (task iraf.imcoords.ccmap())

Make sure you know where all the necessary libraries are and that they live in
your environment.  The pipeline uses an environment variable called EXTROOT
to cover installs of external packages; thus pipeline C programs will look for
include files in $EXTROOT/include/ when compiling and for libraries in
$EXTROOT/lib/ when linking, and will install in $EXTROOT/bin/ by default when
done.  Put external packages here if they're not globally installed.
The pipeline assumes that all external binaries are somewhere in your $PATH;
subpipe_env.csh adds $EXTROOT/bin to $PATH.  If any of these are missing,
edit subpipe_env.csh to make sure they're covered.  I suggest putting $EXTROOT
in a sensible place like $HOME/local/.

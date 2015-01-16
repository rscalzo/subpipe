
#Define enviroment variables
#should really be included in user configuration

#calibration directories
setenv STAP_cal

#calibration names
setenv STAP_ff $STAP_cal/ff/ #ext/filename_filter_ext.fits
setenv STAP_fr $STAP_cal/fr/ #ext/filename_filter_ext.fits

#image directories
setenv STAP_temp      #/filename.fits
setenv STAP_raw       #ext/filename_ext.fits
setenv STAP_new       #ext/filename_ext.fits
setenv STAP_ref       #ext/filename_ext.fits
setenv STAP_diff      #ext/filename_ext.fits
setenv STAP_arch      #date/ext/filename_ext.fits

#product directories
setenv STAP_thumb     #date/filename_ext_HHMMSS+DDMMSS.jpg

#log directory
setenv STAP_logs      #date/processid_YYYYMMDDHHMMSS.log

#configuration root
setenv STAP_etc

#configuration files
setenv STAP_sexconf STAP_etc/stap.sex
setenv STAP_sexpar  STAP_etc/stap.par

#path for running version of pipeline
setenv STAP_source

#control files
setenv STAP_rawlist        #text
setenv STAP_newlist

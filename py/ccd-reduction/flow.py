'''

GOAL:
- to reduce a night's worth of ccd data through flatfielding

REQUIREMENTS FOR THE DATA:

In its current version, this program assumes that you have the following calibration files
- flatfield
- bias
- darks

In the near future, we will enable you to indicate that you took only darks, or only bias frames.

Be sure to move any files that are bod (for whatever reason) into a subdirectory (e.g. junk/).


EXAMPLE:

* Once the bad files are removed from the directory, type

    python ~/github/siena-observatory/ccd-reduction/flow.py

* Then, fix filenames so science image names end in .fits instead of .fit

    cd PROCESSED
    python ~/github/siena-observatory/ccd-reduction/fixname.py

* And fix image headers to provide the information that scamp needs

    python ~/github/siena-observatory/ccd-reduction/fixheader.py


REQUIRED MODULES:
- ccdproc
- astropy
- numpy (is that worth listing?)


################################################################
NOTES:
#### keep in separate routine ######
# For solving WCS and making mosaics:
# - sextractor
# - scamp
# - swarp
#####################################



References:
https://media.readthedocs.org/pdf/ccdproc/latest/ccdproc.pdf


testing with data from

/Users/rfinn/Dropbox/Siena/observatory/images/reduced/2018-12-04

and sky flats from

/Users/rfinn/Dropbox/Siena/observatory/images/reduced/2018-12-05


TO DO:
- add option to subtract bias
- add option to subtract dark
- use dark with matching exposure time if it exists


################################################################

'''

import ccdproc
from ccdproc import ImageFileCollection
from astropy.nddata import CCDData
import os
from astropy.io import fits
import numpy as np
import astropy.units as u
import argparse
# move any bad files to junk



#####################################################
# set the following variables to 1
# when starting the reduction of a night
# If you need to repeat a later stage, you
# can skip the earlier steps (cosmic ray rejection
# in particular is time consuming.
#####################################################



parser = argparse.ArgumentParser(description ='Process images through flatfielding')
parser.add_argument('--nobias', dest = 'nobias', default = False, action = 'store_true', help = 'Skip bias subtraction')
parser.add_argument('--nodark', dest = 'nodark', default = False, action = 'store_true', help = 'Skip dark subtraction')
args = parser.parse_args()

zerocombine = 0
runzapcosmic = 0
darkcombine = 0
flatcombine = 0
process_science = 1
cleanup = 0

# replace these with real values
# for Siena SBIG STL-11000M CCD
stl_rdnoise = 13. # e-, rms
stl_gain = .8#e-/ADU unbinned
stl_gain = 1.6#e-/ADU binned


gain = stl_gain*u.electron/u.adu
rdnoise = stl_rdnoise*u.electron

# check theses using 1/28/19 data
ccdkeyword={'light':'Light Frame','dark':'Dark Frame','flat':'Flat Field','bias':'Bias Frame'}



def zapcosmic(file_list):
    i=0
    nfiles = len(file_list)
    for f in file_list:
        print ('ZAPPING COSMIC RAYS FOR FILE %i OF %i'%(i,nfiles))
        with fits.open(f) as hdu1:
            print ('working on ',f)

            # convert data to CCDData format and save header
            ccd = CCDData(hdu1[0].data, unit=u.adu)
            header = hdu1[0].header
            crimage = ccdproc.cosmicray_lacosmic(ccd, gain = float(gain.value), readnoise = float(rdnoise.value))
            header['HISTORY'] = '= Cosmic rays rejected using ccdproc.cosmicray_lacosmic '
            fits.writeto('z'+f,crimage,header)
            hdu1.close()
        i += 1
        print ('\n')

# make an image collection of all the files in the data directory
# z is prefix for galaxies that already have cosmic ray rejection
ic = ImageFileCollection(os.getcwd(),keywords='*',glob_include='*.fit*')


# combine bias frames

if zerocombine:
    # select all files with imagetyp=='bias'
    bias_files = ic.files_filtered(imagetyp = ccdkeyword['bias'])
    # feed list into ccdproc.combine, output bias

    master_bias = ccdproc.combine(bias_files,method='average',sigma_clip=True,unit=u.adu)
    gaincorrected_master_bias = ccdproc.gain_correct(master_bias,float(gain))
    print('writing fits file for master bias')
    gaincorrected_master_bias.write('bias-combined.fits',overwrite=True)
else:
    print('not combining zeros')
    print('\t reading in bias-combined.fits instead')
    hdu1 = fits.open('bias-combined.fits')
    header = hdu1[0].header
    gaincorrected_master_bias = CCDData(hdu1[0].data, unit=u.electron, meta=header)
    hdu1.close()
###################################################
# COMBINE DARKS
###################################################
# identify darks and combine darks with longest exposure time

# change to use the dark with the closest exposure time

exptimes = np.array(ic.values('exptime'))

image_types = np.array(ic.values('imagetyp'))

if darkcombine:

    # select all files with imagetyp=='dark'
    # want to read in one set of long exposure dark frames, like 120 s
    # observers should take a set of darks that correspond to longest exposure time
    # e.g. 120s

    set_exptime=set(exptimes[image_types == ccdkeyword['dark']])
    dark_exptimes_set = np.array(list(set_exptime),'f')
    max_exposure = max(exptimes[image_types == ccdkeyword['dark']])
    '''
    * subtract bias from combined dark frames

    * skipping for now because bias has more counts that the darks - go figure!

    
    gaincorrected_master_dark_bias = ccdproc.subtract_bias(gaincorrected_master_dark, gaincorrected_master_bias)

    print('writing fits file for master dark')
    gaincorrected_master_dark.write('dark-combined.fits',overwrite=True)
    '''
    # combine darks of a given exposure time
    for expt in set_exptime:
        # skip max exposure time - already combined

        dark_files = ic.files_filtered(imagetyp = ccdkeyword['dark'], exptime = expt)
        dark = ccdproc.combine(dark_files,method='average',sigma_clip=True,unit=u.adu)
        gaincorrected_dark = ccdproc.gain_correct(dark,float(gain))
        '''
        * subtract bias from combined dark frames
    
        gaincorrected_dark_bias = ccdproc.subtract_bias(gaincorrected_dark, gaincorrected_master_bias)
        '''
        print('writing fits file for dark t = ',int(expt))
        gaincorrected_dark.write('dark-combined-'+str(int(expt))+'.fits',overwrite=True)


else:
    print('not combining darks')
    print('\t reading in dark-combined-*.fits instead')
    darkcollection = ImageFileCollection(os.getcwd(),keywords='*',glob_include='dark-combined*.fit*')
    dark_exptimes = np.array(ic.values('exptime'))
    dimage_types = np.array(ic.values('imagetyp'))

    dark_exptimes_set=np.array(list(set(dark_exptimes[image_types == ccdkeyword['dark']])),'f')



print('runzapcosmic = ',runzapcosmic)
if runzapcosmic:

###################################################
# COSMIC RAY REJECTION
###################################################
# run LA cosmic on flat field images

    flat_files = ic.files_filtered(imagetyp = ccdkeyword['flat'])

    zapcosmic(flat_files)

    # run LA cosmic on science images

    sci_files = ic.files_filtered(imagetyp = ccdkeyword['light'])
    zapcosmic(sci_files)


# make an image collection of all the files in the data directory
# z is prefix for galaxies that already have cosmic ray rejection
icz = ImageFileCollection(os.getcwd(),keywords='*',glob_include='z*.fit*')

filters = icz.values('filter')
image_types = np.array(icz.values('imagetyp'))


#####################################################
###### PROCESS FLAT FIELD IMAGES
#####################################################


# create a set of unique filters
all_filters = set(np.array(filters))
# gets rid of 0.0, which is what the dark 'filter' will come through as
try:
    all_filters.remove('0.0')
except KeyError:
    print('no bias frames found after zapping')


# read in available dark files
# dark-bias-
if flatcombine:     
    # loop through filters, combine all the flats in that filter
    for filt in all_filters:
        # get of all flats in this filter
        flat_files = icz.files_filtered(imagetyp = ccdkeyword['flat'], filter = filt)
        #print(filt,flat_files)
        if len(flat_files) == 0:
            print('WARNING: no flat images for filter ',filt)
            print('skipping this filter')
            continue
        # combine images into one flat
        print('combining flats for filter ',filt)
        my_master_flat = ccdproc.combine(flat_files,method='median',unit=u.adu)
        gaincorrected_master_flat = ccdproc.gain_correct(my_master_flat,gain)
        # subtract bias
        # subtract the scaled dark
        print('subtracting bias and scaled dark from combined flat for filter ',filt)
        #master_flat_bias_dark = ccdproc.ccd_process(gaincorrected_master_flat, readnoise=rdnoise, master_bias = gaincorrected_master_bias, dark_frame=gaincorrected_master_dark, exposure_key='exposure', exposure_unit=u.second, dark_scale=True, gain_corrected=True)

        #################################################
        ######## FIND DARK WITH CLOSEST EXPOSURE TIME
        #################################################
        # find dark with closest exposure time
        flat_exptime = my_master_flat.header['EXPTIME']
        # find dark with closest exposure time
        delta_t = abs(flat_exptime - dark_exptimes_set)
        closest_dark = dark_exptimes_set[delta_t == min(delta_t)]
        # open the appropriate dark
        hdu1 = fits.open('dark-combined-'+str(int(closest_dark[0]))+'.fits')
        gaincorrected_dark = CCDData(hdu1[0].data, unit=u.electron, meta=header)
        hdu1.close()

                
        #master_flat_dark = ccdproc.ccd_process(gaincorrected_master_flat, readnoise=rdnoise, dark_frame=gaincorrected_dark, exposure_key='exposure', exposure_unit=u.second, dark_scale=True, gain_corrected=True)
        master_flat_dark = ccdproc.subtract_dark(gaincorrected_master_flat, gaincorrected_dark, exposure_time='exposure', exposure_unit=u.second,scale=False)

        print('writing out combined flat for filter ',filt)
        # write output    
        master_flat_dark.write('flat-'+filt+'.fits',overwrite=True)
else:
    print('skipping flat combine')

#####################################################
###### PROCESS SCIENCE IMAGES
#####################################################
if process_science:
    print('\n Processing science frames!!!')
    for filt in all_filters:
        sci_files = icz.files_filtered(imagetyp = ccdkeyword['light'], filter = filt) 
        master_flat = 'flat-'+filt+'.fits'
        if not(os.path.exists(master_flat)):
            print('WARNING: no flat found for filter ',filt)
            print('skipping images in this filter')
            continue
        hdu = fits.open(master_flat)
        gaincorrected_master_flat = CCDData(hdu[0].data, unit=u.electron)
        hdu.close()
        # write out image with prefix 'fdb' for flat, dark, bias corrected
        for f in sci_files:
            with fits.open(f) as hdu1:
                print ('working on ',f)
                print ('imagetype = ',hdu1[0].header['IMAGETYP'])
                # convert data to CCDData format and save header
                header = hdu1[0].header
                #ccd = CCDData(hdu1[0].data, unit=u.adu,meta={'exposure':header['exposure']})
                ccd = CCDData(hdu1[0].data, unit=u.adu,meta=header)

                #newccd = ccdproc.ccd_process(ccd,error=True, gain=gain, readnoise=rdnoise, master_bias = gaincorrected_master_bias, dark_frame=gaincorrected_master_dark, exposure_key='exposure', exposure_unit=u.second, dark_scale=True,master_flat = gaincorrected_master_flat, gain_corrected=True)


                #################################################
                ######## FIND DARK WITH CLOSEST EXPOSURE TIME
                #################################################
                # find dark with closest exposure time
                sci_exptime = header['EXPTIME']
                # find dark with closest exposure time
                delta_t = abs(sci_exptime - dark_exptimes_set)
                closest_dark = dark_exptimes_set[delta_t == min(delta_t)]
                # open the appropriate dark
                hdu1 = fits.open('dark-combined-'+str(int(closest_dark[0]))+'.fits')
                gaincorrected_dark = CCDData(hdu1[0].data, unit=u.electron, meta=header)
                hdu1.close()
                
                newccd = ccdproc.ccd_process(ccd,error=True, gain=gain, readnoise=rdnoise, dark_frame=gaincorrected_dark, exposure_key='exposure', exposure_unit=u.second,master_flat = gaincorrected_master_flat, gain_corrected=True)

                header['HISTORY'] = '= Processing by ccdproc: dark, flat '
                #fits.writeto('fdb'+f,newccd,header, overwrite=True)
                fits.writeto('fd'+f,newccd,header, overwrite=True)
                hdu1.close()

else:
    print('skipping science frames')    
'''
# then move on to running scamp and swarp to solve WCS
# and make mosaics

# let's start 01-28-2019 data 

'''
#####################################################
###### CLEANING UP
#####################################################

if cleanup:
    # moved processed a subdirectory
    dirname = 'PROCESSED'
    prefix = 'fdb'
    if not(os.path.exists(dirname)):
        os.mkdir(dirname)
    os.system('mv '+prefix+'*.fit* '+dirname+'/.')

    # move zapped to a subdirectory
    dirname = 'ZAPPED'
    prefix = 'z'
    if not(os.path.exists(dirname)):
        os.mkdir(dirname)
    os.system('mv '+prefix+'*.fit* '+dirname+'/.')


    # move originals to subdirectory
    dirname = 'ORIGINALS'
    prefix =''
    if not(os.path.exists(dirname)):
        os.mkdir(dirname)
    os.system('mv '+prefix+'*.fit* '+dirname+'/.')
    



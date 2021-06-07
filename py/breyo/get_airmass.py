#!/usr/bin/env python

'''
Determine the airmass of an image and add it to the image header.



example usage:

python ~/github/breyo/py/breyo/fixcoords.py --filestring zfbp-M39 

This will update all files zfbp-M39*.fits


NOTES: 

'''

import argparse
from astropy.io import fits
import glob
import sys

import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz

if __name__ == '__main__':
    ###########################
    ##### SET UP ARGPARSE
    ###########################

    parser = argparse.ArgumentParser(description ='Program to replace RA and DEC in STL11000M header with user provided input')
    parser.add_argument('--filestring', dest = 'filestring', default=None, help = 'filestring for glob. e.g. to get *M39*.fits, enter M39')
    parser.add_argument('--airmasskeyword', dest = 'seczkeyword', default = 'AIRMASS', help = 'header keyword for airmass.  Default is AIRMASS')    

    args = parser.parse_args()
    matchstring = '*'+args.filestring+'*.fits'
    #print(matchstring)
    #print(args.ra)
    #print(args.dec)
    latitude = 42 + 39/60 +9/3600
    longitude = -1*(73+45/60 + 22/3600)
    breyo = EarthLocation(lat=latitude*u.deg,lon=longitude*u.deg,height=200*u.m)

    if args.filestring is not None:
        files = glob.glob(matchstring)
        print('got ',len(files),' files to update')
        for f in files:
            data,header = fits.getdata(f,header=True)
            # get RA and DEC
            ra = header['CRVAL1']
            dec = header['CRVAL2']
            target_coord = SkyCoord(ra,dec,unit='deg')
            tobs = header['JD']
            t = Time(tobs, format='jd')
            
            frame_tobs = AltAz(obstime=t,location=breyo)
            target_altaz = target_coord.transform_to(frame_tobs)

            
            
            header.set('AIRMASS',target_altaz.secz)
            fits.writeto(f,data,header,overwrite=True)
            
    else:
        print('no filestring provided')
        sys.exit()


#!/usr/bin/env python

'''
replace RA and DEC with provided values

example usage:


'''

import argparse
from astropy.io import fits
import glob
import sys

if __name__ == '__main__':
    ###########################
    ##### SET UP ARGPARSE
    ###########################

    parser = argparse.ArgumentParser(description ='Program to replace RA and DEC in STL11000M header with user provided input')
    parser.add_argument('--filestring', dest = 'filestring', default=None, help = 'filestring for glob. e.g. *M39*.fits')
    parser.add_argument('--ra', dest = 'ra', default = None, help = "RA of object. e.g. '21 32 12'")
    parser.add_argument('--dec', dest = 'dec', default = None, help = 'DEC of object. e.g. "+48 26 00')    
    parser.add_argument('--rakeyword', dest = 'rakey', default = 'OBJCTRA', help = "header keyword for RA.  Default is OBJCTRA.")
    parser.add_argument('--deckeyword', dest = 'deckey', default = 'OBJCTDEC', help = 'header keyword for DEC.  Default is OBJCTDEC')    

    args = parser.parse_args()

    if args.filestring is not None:
        files = glob.glob(args.filestring)
        print('got ',len(files),' files to update')
        if args.ra is not None:
            for f in files:
                data,header = fits.getdata(f,header=True)
                header.set(args.rakeyword,args.ra)
                fits.writeto(f,data,header,overwrite=True)
        if args.dec is not None:
            for f in files:
                data,header = fits.getdata(f,header=True)
                header.set(args.deckeyword,args.dec)
                fits.writeto(f,data,header,overwrite=True)
            
    else:
        print('no filestring provided')
        sys.exit()


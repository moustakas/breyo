#!/usr/bin/env python


import os

import argparse
import glob
import sys

if __name__ == '__main__':
    ###########################
    ##### SET UP ARGPARSE
    ###########################

    parser = argparse.ArgumentParser(description ='Program to change prefix of filename')
    parser.add_argument('--orig', dest = 'orig', default=None, help = 'original prefix that you want to change')
    parser.add_argument('--new', dest = 'new', default = None, help = "new prefix")

    args = parser.parse_args()
    print(args.orig,args.new)
    matchstring = args.orig+'*.fits'
    #print(matchstring)
    #print(args.ra)
    #print(args.dec)    
    
    if args.orig is not None:
        files = glob.glob(matchstring)
        print('got ',len(files),' files to update')
        
        for f in files:
            outfile = f.replace(args.orig,args.new)
            os.rename(f,outfile)
            
    else:
        print('no filestring provided')
        sys.exit()


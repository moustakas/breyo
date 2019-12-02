'''

EXAMPLE:

python ~/github/siena-observatory/ccd-reduction/changeheader.py --filestring 'blank1' --keyword 'IMAGETYP' --newvalue 'Flat Field'


'''


from astropy.io import fits
import argparse
import glob

parser = argparse.ArgumentParser(description ='change header field in a file')
parser.add_argument('--filestring', dest = 'filestring', default = 'h', help = 'string to use to get input files (default = "h" which grabs all files "h*o00.fits")')

parser.add_argument('--keyword', dest = 'keyword', help = 'Header keyword to update')
parser.add_argument('--newvalue', dest = 'new', help = 'new value for header keyword')

args = parser.parse_args()

files = glob.glob(args.filestring+'*')

for f in files:
    t = fits.open(f)
    t[0].header[args.keyword] = args.new
    print('updating file ',f)
    print('image size = ',t[0].header['naxis1'],t[0].header['naxis2'])
    print('')
    t.writeto(f, overwrite = True)

    




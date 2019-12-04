"""
breyo.reduce
============

Reduction code.

"""
import os
import glob

import numpy as np

def make_masterbias(overwrite=False):
    
    masterbiasfile = os.path.join(reduxdir, 'master-bias.fits')
    if not os.path.isfile(masterbiasfile) or overwrite:
        print('Building the master bias frame.')
        allbiasfiles = glob.glob('{}/bias-????.fit'.format(datadir))
        # build the median stack
        allbias = np.array([fits.getdata(biasfile) for biasfile in allbiasfiles])
        masterbias = np.median(allbias, axis=0)
        # read one header and write out
        hdr = fits.getheader(allbiasfiles[0])
        hdr['NBIAS'] = len(allbias)
        print('Writing {}'.format(masterbiasfile))
        fits.writeto(masterbiasfile, masterbias.astype('f4'), header=hdr, overwrite=True)
    else:
        print('Reading {}'.format(masterbiasfile))
        masterbias = fits.getdata(masterbiasfile)
    return masterbias

def make_masterflat(masterbias, overwrite=False):
    allmasterflat = []
    for filt in ('b', 'v', 'r'):
        masterflatfile = os.path.join(reduxdir, 'master-flat-{}.fits'.format(filt))
        if not os.path.isfile(masterflatfile) or overwrite:
            print('Working on filter {}'.format(filt))
            # Get the filenames.
            allflatfiles = glob.glob('{}/skyflat-????{}.fit'.format(datadir, filt))
            # Read the data, one image at a time and normalize by the mean counts.
            allflat = []
            for flatfile in allflatfiles:
                print('  Reading {}'.format(flatfile))
                thisflat = fits.getdata(flatfile) * 1.0 # convert to float
                thisflat -= masterbias
                thisflat /= np.median(thisflat)
                allflat.append(thisflat)
            # median stack
            thismasterflat = np.median(allflat, axis=0)
            # read one header and write out
            hdr = fits.getheader(allflatfiles[0])
            hdr['NFLAT'] = len(allflat)
            print('Writing {}'.format(masterflatfile))
            fits.writeto(masterflatfile, thismasterflat.astype('f4'), header=hdr, overwrite=True)
        else:
            print('Reading {}'.format(masterflatfile))
            thismasterflat = fits.getdata(masterflatfile)
        allmasterflat.append(thismasterflat)
    return allmasterflat

# Commented out IPython magic to ensure Python compatibility.
# %time masterbias = make_masterbias(overwrite=do_masterbias)

# Commented out IPython magic to ensure Python compatibility.
# %time masterflat_b, masterflat_v, masterflat_r = make_masterflat(masterbias, overwrite=do_masterflats)

"""### Reduce the data!

Reduce all the science data by subtracting the bias and dividing by the flat-field (in the appropriate filter).
"""

def reduce_alldata(bias, bflat, vflat, rflat, prefix='fb', overwrite=False):
    t0 = time.time()
    allimgfiles = glob.glob('{}/NGC6830-*.fit'.format(datadir))
    print('Found {} images in {}'.format(len(allimgfiles), datadir))
    for imgfile in allimgfiles:
        imgoutfile = os.path.join(reduxdir, '{}-{}'.format(prefix, os.path.basename(imgfile).replace('.fit', '.fits')))
        if not os.path.isfile(imgoutfile) or overwrite:
            img = fits.getdata(imgfile) * 1.0 # convert to float
            # Subtract the bias and then divide by the appropriate flat-field
            img -= bias
            hdr = fits.getheader(imgfile)
            filt = hdr['FILTER'].lower().strip()
            if filt == 'blue':
                img /= bflat
            elif filt == 'green':
                img /= vflat
            elif filt == 'red':
                img /= rflat
            else:
                print('Unrecognized filter {}!'.format(filt))
            # Trim the first 100 columns, which are funky.
            img = img[:, 100:]
            # Finally multiply by the gain and divide by the exposure time.
            img = img * hdr['EGAIN'] / hdr['EXPTIME'] # [e-/sec]
            print('Writing {}'.format(imgoutfile))
            fits.writeto(imgoutfile, img.astype('f4'), header=hdr, overwrite=True)
    print('All done in {:.2f} sec'.format(time.time() - t0))

def get_sky_background(img, verbose=True):
    from astropy.stats import SigmaClip
    from photutils import Background2D, MedianBackground
    sigma_clip = SigmaClip(sigma=3.)
    bkg = Background2D(img, (100, 100), filter_size=(5, 5), sigma_clip=sigma_clip, 
                       bkg_estimator=MedianBackground())
    if verbose:
        print('Sky background median = {:.3f}, rms = {:.3f} ADU.'.format(bkg.background_median, bkg.background_rms_median))
    return bkg

def find_stars(image, imgfile, fwhm=3.0, nsigma=3, sigma=None, 
               verbose=True, overwrite=False):
    from astropy.table import Table
    
    starsfile = os.path.join(reduxdir, 'stars-{}'.format(os.path.basename(imgfile)))
    if not os.path.isfile(starsfile) or overwrite:
        from photutils import DAOStarFinder
        if sigma is None:
            sigma = np.std(image)

        daofind = DAOStarFinder(fwhm=fwhm, threshold=nsigma * sigma)
        srcs = daofind(image)
        # reverse-sort by flux 
        srcs.sort('flux')
        srcs.reverse()
        if verbose:
            print('Found {} sources'.format(len(srcs)))

        print('Writing {} stars to {}'.format(len(srcs), starsfile))
        srcs.write(starsfile, overwrite=True)
    else:
        srcs = Table.read(starsfile)
        print('Read {} stars from {}'.format(len(srcs), starsfile))
    return srcs

def get_astrometry(imgfile, srcs=None, api_key=None, prefix='w', overwrite=False):
    from astropy.io import fits
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from astroquery.astrometry_net import AstrometryNet  

    wcsfile = os.path.join(reduxdir, '{}{}'.format(prefix, os.path.basename(imgfile)))
    if not os.path.isfile(wcsfile) or overwrite:
        img, hdr = fits.getdata(imgfile, header=True)

        # Initialize the API.
        ast = AstrometryNet()
        if api_key:
            ast.api_key = api_key
        #ast.show_allowed_settings()

        # Get the initial position center based on the header.
        c = SkyCoord(hdr['OBJCTRA']+hdr['OBJCTDEC'], unit=(u.hourangle, u.deg))
        print('Initial RA, Dec = {:.5f}, {:.5f}'.format(c.ra.value, c.dec.value))

        # Query the astrometry.net engine!
        t0 = time.time()
        wcshdr = ast.solve_from_source_list(
            srcs['xcentroid'], srcs['ycentroid'], hdr['naxis1'], hdr['naxis2'],
            center_ra=c.ra.value, center_dec=c.dec.value, radius=15/60.0, 
            scale_type='ev', scale_est=0.8, scale_err=10, scale_units='arcsecperpix',
            crpix_center=True)
        print('Total time = {:.3f} min'.format((time.time() - t0)/60))

        # update the original header
        for key in wcshdr.keys():
            if key not in hdr and key != 'COMMENT' and key != 'HISTORY':
                hdr[key] = wcshdr[key]

        print('Writing {}'.format(wcsfile))
        fits.writeto(wcsfile, img, header=wcshdr, overwrite=True)        
    else:
        wcshdr = fits.getheader(wcsfile)

    return wcsfile, wcshdr

def undistort_image(imgfile, wcsfile, pixscale=0.8, display=False):
    """Little wrapper script to project (i.e., undistort) a distorted 
    image onto a tangent plane.
    
    """
    from astropy.io import fits
    from astropy.wcs import WCS
    from reproject import reproject_interp
    
    img, hdr = fits.getdata(imgfile, header=True)
    wcshdr = fits.getheader(wcsfile)

    # Create a header describing a tangent plane.
    for key in ('WCSAXES', 'CTYPE1', 'CTYPE2', 'EQUINOX', 'LONPOLE', 'LATPOLE', 
                'CRVAL1', 'CRVAL2', 'CRPIX1', 'CRPIX2', 'CUNIT1', 'CUNIT2',
                'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2', 'IMAGEW', 'IMAGEH'):
        hdr[key] = wcshdr[key]

    # Tangent projection
    #wwcshdr = WCS(wcshdr)
    #whdr = WCS(hdr)    
    #whdr.all_world2pix(wwcshdr.all_pix2world(np.array([(1, 1), (1, dim[0]), (dim[1], dim[0]), (dim[1], 0)]), 1), 1)
    hdr['CTYPE1'] = 'RA---TAN'
    hdr['CTYPE2'] = 'DEC--TAN'
    # update the pixel scale [arcsec/pix] to the desired constant with 
    # no rotation
    hdr['CD1_2'] = 0.0
    hdr['CD2_1'] = 0.0
    hdr['CD1_1'] = -pixscale / 3600.0
    hdr['CD2_2'] =  pixscale / 3600.0

    # now project
    img_proj, footprint = reproject_interp((img, wcshdr), hdr)

    if display:
        #display_image(img, projection=WCS(wcshdr))
        display_image(img_proj, projection=WCS(hdr))

    return img_proj, footprint

def find_stars(image, imgfile, fwhm=3.0, nsigma=3, sigma=None, 
               verbose=True, overwrite=False):
    from astropy.table import Table
    
    starsfile = os.path.join(reduxdir, 'stars-{}'.format(os.path.basename(imgfile)))
    if not os.path.isfile(starsfile) or overwrite:
        from photutils import DAOStarFinder
        if sigma is None:
            sigma = np.std(image)

        daofind = DAOStarFinder(fwhm=fwhm, threshold=nsigma * sigma)
        srcs = daofind(image)
        # reverse-sort by flux 
        srcs.sort('flux')
        srcs.reverse()
        if verbose:
            print('Found {} sources'.format(len(srcs)))

        print('Writing {} stars to {}'.format(len(srcs), starsfile))
        srcs.write(starsfile, overwrite=True)
    else:
        srcs = Table.read(starsfile)
        print('Read {} stars from {}'.format(len(srcs), starsfile))
    return srcs

def get_panstarrs_catalog(imgwcs, imgfile, radius=0.2, rfaint=17, region=False):
    from astroquery.mast import Catalogs
    ra0, dec0 = imgwcs.wcs.crval
    print('Querying Pan-STARRS catalog with radius={:.3f} deg and central coordinates RA,Dec={:.5f},{:.5f}'.format(
        radius, ra0, dec0))
    if region:
        allcat = Catalogs.query_region('{} {}'.format(ra0, dec0), radius=radius,
                                       catalog='PANSTARRS', data_release='dr2', 
                                       table='mean')#, rMeanPSFMag=[12, 22])
    else:
        allcat = Catalogs.query_criteria(coordinates='{} {}'.format(ra0, dec0), radius=radius,
                                         catalog='PANSTARRS', data_release='dr2', 
                                         table='mean',
                                         columns=['objID', 'raMean', 'decMean',
                                                  'gMeanPSFMag', 'rMeanPSFMag', 'iMeanPSFMag', 'zMeanPSFMag'],
                                         gMeanPSFMag=[('lte', 18), ('gte', 12)],
                                         rMeanPSFMag=[('lte', 18), ('gte', 12)],
                                         iMeanPSFMag=[('lte', 18), ('gte', 12)],
                                         zMeanPSFMag=[('lte', 18), ('gte', 12)],
                                         sort_by=[("asc", "rMeanPSFMag")])
    
    #rmag = allcat['rMeanPSFMag']
    #good = np.isfinite(rmag) * rmag < rfaint
    #cat = allcat[good]
    #print('Keeping {}/{} Pan-STARRS sources.'.format(len(cat), len(allcat)))
    return allcat


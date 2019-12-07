"""
breyo.io
========

Code for deriving photometric zeropoints.

"""
import numpy as np

def get_sky_background(img, verbose=False):
    from astropy.stats import SigmaClip
    from photutils import Background2D, MedianBackground

    sigma_clip = SigmaClip(sigma=3.)
    bkg = Background2D(img, (100, 100), filter_size=(5, 5), sigma_clip=sigma_clip, 
                       bkg_estimator=MedianBackground())
    if verbose:
        print('Sky background median = {:.3f}, rms = {:.3f}'.format(
            bkg.background_median, bkg.background_rms_median))
        
    return bkg

def find_stars(img, fwhm=3.0, nsigma=5, verbose=False):
    from photutils import DAOStarFinder

    bkg = get_sky_background(img, verbose=verbose)
    sigma = bkg.background_rms_median
    
    daofind = DAOStarFinder(fwhm=fwhm, threshold=nsigma * sigma)
    srcs = daofind(img - bkg.background)
    #import pdb ; pdb.set_trace()
    #srcs.rename_column('objid', 'objno')
    
    ## reverse-sort by flux 
    #srcs.sort('flux')
    #srcs.reverse()
    if verbose:
        print('Found {} stars'.format(len(srcs)))
        
    return srcs

def get_panstarrs_catalog(imgwcs, radius=0.2, verbose=False):
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from astroquery.mast import Catalogs

    ra0, dec0 = imgwcs.wcs.crval
    coords = SkyCoord(ra0, dec0, unit=u.deg, frame='icrs')
    
    if verbose:
        print('Querying Pan-STARRS {:.3f} deg around RA, Dec={:.5f}, {:.5f} '.format(radius, ra0, dec0))
        
    cat = Catalogs.query_criteria(coordinates=coords, radius=radius,
                                  catalog='PANSTARRS', data_release='dr2', 
                                  table='mean',
                                  columns=['objID', 'raMean', 'decMean',
                                           'gMeanPSFMag', 'rMeanPSFMag', 'iMeanPSFMag', 'zMeanPSFMag'],
                                           gMeanPSFMag=[('lte', 18), ('gte', 8)],
                                           rMeanPSFMag=[('lte', 18), ('gte', 8)],
                                           iMeanPSFMag=[('lte', 18), ('gte', 8)],
                                           zMeanPSFMag=[('lte', 18), ('gte', 8)])
                                           #sort_by=[("asc", "rMeanPSFMag")]

    # http://legacysurvey.org/dr8/description/#photometry
    gi = cat['gMeanPSFMag'] - cat['iMeanPSFMag']
    keep = np.where( (gi > 0.4) * (gi < 2.7) )[0]
    cat = cat[keep]
    
    return cat

def reproject_one_image(hdu, pixscale=0.8):
    """Little wrapper script to project (i.e., undistort) a distorted 
    image onto a tangent plane.
    
    """
    from astropy.wcs import WCS
    from reproject import reproject_interp

    img, hdr = hdu.data, hdu.header

    dim = img.shape
    imgwcs = WCS(hdr)

    for ii, hdr1 in enumerate(hdr.cards):
        try:
            if 'Original key: "END"' in hdr1[1]:
                cut = ii
                break
        except:
            pass
    newhdr = hdr[:cut]

    # Create a header describing a tangent plane.
    for key in ('WCSAXES', 'CTYPE1', 'CTYPE2', 'EQUINOX', 'LONPOLE', 'LATPOLE', 
                'CRVAL1', 'CRVAL2', 'CRPIX1', 'CRPIX2', 'CUNIT1', 'CUNIT2',
                'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2', 'IMAGEW', 'IMAGEH'):
        hdr[key] = wcshdr[key]

    # Tangent projection
    #whdr = WCS(hdr)    
    #corn = whdr.all_world2pix(whdr.all_pix2world(np.array([(1, 1), (1, dim[0]), (dim[1], dim[0]), (dim[1], 0)]), 1), 1))


    hdr['NAXIS1'] = 2000
    hdr['NAXIS2'] = 2000
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

    return footprint, img_proj, hdr

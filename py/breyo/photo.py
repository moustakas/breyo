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
    gi = allcat['gMeanPSFMag'] - allcat['iMeanPSFMag']
    keep = np.where( (gi > 0.4) * (gi < 2.7) )[0]
    cat = cat[keep]
    
    return cat

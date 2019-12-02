"""
breyo.reduce
============

Reduction code.

"""
import numpy as np

def display_image(img, minclip=5, maxclip=95, label=None, cmap='Greys_r', 
                  srcs=None, projection=None, calibrated=False):
    """Simple wrapper to display an image.
    
    """
    from astropy.visualization import AsinhStretch as Stretch
    from astropy.visualization import ZScaleInterval as Interval
    from astropy.visualization.mpl_normalize import ImageNormalize

    #from astropy.visualization import simple_norm
    #norm = simple_norm(img, min_percent=minclip, max_percent=maxclip)

    interval = Interval(contrast=0.5)
    vmin, vmax = interval.get_limits(img)
    norm = ImageNormalize(interval=interval, stretch=Stretch(a=0.9))

    fig, ax = plt.subplots(figsize=(8, 8), subplot_kw={'projection': projection})
    im = ax.imshow(img, origin='lower', norm=norm, cmap=cmap,
                   vmin=vmin, vmax=vmax)
    if projection:
        ax.coords.grid(color='red')
        ax.coords['ra'].set_axislabel('Right Ascension')
        ax.coords['dec'].set_axislabel('Declination')      
    else:
        ax.set_xlabel('Column Number (pixels)')
        ax.set_ylabel('Row Number (pixels)')

    # Mark the locations of stars.
    if srcs:
        from photutils import CircularAperture
        pos = np.transpose((srcs['xcentroid'], srcs['ycentroid']))
        aps = CircularAperture(pos, r=12.)
        aps.plot(color='red', lw=1.5, alpha=0.6, axes=ax)
      
    # Make room for the colorbar
    fig.subplots_adjust(right=0.8)
    cax = fig.add_axes([0.85, 0.28, 0.05, 0.45])
    c = plt.colorbar(im, cax=cax)
    if label:
        c.set_label(label)
    else:
        if calibrated:
            c.set_label(r'Intensity ($e^{-}/s$)')
        else:
            c.set_label('Intensity (ADU)')

"""## Create calibration files if they don't already exist.

Build the master bias frames and flat-fields.
"""

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

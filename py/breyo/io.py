"""
breyo.io
========

I/O code.

"""
import os
import astropy.units as u

def get_rawdir(night):
    return os.path.join(os.getenv('BREYO_DATA_DIR'), 'raw', night)
    
def get_reduxdir(night=None):
    reduxdir = os.path.join(os.getenv('BREYO_DATA_DIR'), 'reduced')
    if night:
        reduxdir = os.path.join(reduxdir, night)
    if not os.path.isdir(reduxdir):
        os.makedirs(reduxdir, exist_ok=True)
    calibdir = os.path.join(reduxdir, 'calib')
    if not os.path.isdir(calibdir):
        os.makedirs(calibdir, exist_ok=True)
    return reduxdir, calibdir

def read_config(verbose=False):
    import yaml

    # not sure why this doesn't work
    #import pkg_resources
    #configfile = pkg_resources.resource_filename('breyo', os.path.join('data', 'SBIG-STL-11000M.yaml'))

    configfile = os.path.join(os.getenv('BREYO_DIR'), 'py', 'breyo', 'data', 'SBIG-STL-11000M.yaml')

    if verbose:
        print('Reading configuration file {}'.format(configfile))
    with open(configfile, 'r') as ff:
        params = yaml.safe_load(ff)

    # attach units (see config file for comments)
    params['ccd']['pixsize'] *= u.micron
    params['ccd']['gain'] *= u.electron / u.adu
    params['ccd']['readnoise'] *= u.electron 

    return params

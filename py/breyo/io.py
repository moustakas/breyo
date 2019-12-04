"""
breyo.io
========

I/O code.

"""
import os
import astropy.units as u

def read_config():
    import yaml

    # not sure why this doesn't work
    #import pkg_resources
    #configfile = pkg_resources.resource_filename('breyo', os.path.join('data', 'SBIG-STL-11000M.yaml'))

    configfile = os.path.join(os.getenv('BREYO_DIR'), 'py', 'breyo', 'data', 'SBIG-STL-11000M.yaml')
    
    print('Reading configuration file {}'.format(configfile))
    with open(configfile, 'r') as ff:
        params = yaml.safe_load(ff)

    # attach units (see config file for comments)
    params['ccd']['pixsize'] *= u.micron
    params['ccd']['gain'] *= u.electron / u.adu
    params['ccd']['readnoise'] *= u.electron 

    return params

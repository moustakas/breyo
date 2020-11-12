# Breyo Observatory Data Reduction & Analysis Pipeline

The repository collects code and documentation used to process imaging and
spectroscopy from the [Siena College](http://siena.edu) [Breyo
Observatory](https://www.siena.edu/departments/physics-and-astronomy/breyo-observatory)
0.7-meter optical telescope.

## Installing

```
conda create --name breyo python=3.6
conda activate breyo
conda install -c astropy astropy ccdproc photutils 
conda install -c conda-forge astrometry
conda install jupyterlab ipython
```

In case the conda installation of astrometry.net doesn't work, try: 

```
brew install wcslib
brew install cairo
brew install netpbm
export NETPBM_LIB="-L/usr/local/lib -lnetpbm"
export NETPBM_INC=-I/usr/local/include/netpbm
export CAIRO_LIB="-L/usr/local/lib -lcairo"
export CAIRO_INC=-I/usr/local/include/cairo
cd /usr/local/share
git clone https://github.com/dstndstn/astrometry.net
cd astrometry.net
ln -s /usr/local/include/netpbm netpbm
make
make extra
make install INSTALL_DIR=/usr/local/astrometry
```

Setting up environment variables: update your .profile file to include
```
export PYTHONPATH=/usr/local/anaconda3/bin/:/home/rfinn/github/breyo/bin/:/home/rfinn/github/breyo/py/
export BREYO_DATA_DIR=/mnt/qnap_home/rfinn/telescope-reduction/
export BREYO_DIR=/home/rfinn/github/breyo/
```

After astromety.net is installed, get the index files.  From your BREYO_DATA_DIR
```
wget --recursive --no-parent http://data.astrometry.net/5000/
```
When the download finishes, rename two directories:
```
mv data.astrometry.net astrometry.net
cd astrometry.net
mv 5000 index-5000
```

## Setting Up Data

The program will assume that the raw data is in BREYO_DATA_DIR/raw, and that each night is stored in a directory.  For example:
```
/mnt/qnap_home/rfinn/telescope_reduction/raw/2020-11-09
```

The program will create BREYO_DATA_DIR/reduced, added a new directory to store the reduced data.  For example:
```
/mnt/qnap_home/rfinn/telescope_reduction/reduced/2020-11-09
```


## Authors

* [**Rose Finn**](https://github.com/rfinn)
* [**John Moustakas**](https://github.com/moustakas)

## Acknowledgments


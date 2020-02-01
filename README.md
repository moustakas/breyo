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


## Authors

* [**Rose Finn**](https://github.com/rfinn)
* [**John Moustakas**](https://github.com/moustakas)

## Acknowledgments


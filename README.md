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

You also need to have the Schlegel dust map installed.  https://github.com/kbarbary/sfdmap/

```
pip install sfdmap

wget https://github.com/kbarbary/sfddata/archive/master.tar.gz
tar xzf master.tar.gz
cd sfddata-master
mkdir maps
cp *.fits maps/.

```


Setting up environment variables: update your .profile file to include
```
export PYTHONPATH=/usr/local/anaconda3/bin/:/home/rfinn/github/breyo/bin/:/home/rfinn/github/breyo/py/
export BREYO_DATA_DIR=/mnt/qnap_home/rfinn/telescope-reduction/
export DUST_DIR=/mnt/qnap_home/rfinn/telescope-reduction/sfddata-master/
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

Create the config file

```
cd astrometry.net
cat index-5000/cfg
#add_path ${BREYO_DATA_DIR}/astrometry.net/index-5000
add_path /Users/ioannis/research/data/breyo/astrometry.net/index-5000
#autoindex
indexset 5000
inparallel

```

## Setting Up Data

The program will assume that the raw data is in BREYO_DATA_DIR/raw, and that each night is stored in a directory.  For example:
```
/mnt/qnap_home/rfinn/telescope_reduction/raw/2020-11-09
```

When running on a virtual machine, the data are stored in /mnt/telescope/RawData/, so I did the following:

```
cd /mnt/qnap_home/rfinn/telescope-reduction
mkdir raw
cd raw
ln -s /mnt/telescope/RawData/20??-??-?? .

```
You will need to rerun this to link new nights of data.

[not sure why I can't just do, from the BREYO_DATA_DIR/

```
ln -s /mnt/telescope/RawData raw
```
I think the paths get confused somehow - will check again at a later date.]


The program will create BREYO_DATA_DIR/reduced, which is a new directory to store the reduced data.  For example:
```
/mnt/qnap_home/rfinn/telescope_reduction/reduced/2020-11-09
```


## Running the reduction pipeline



From BREYO_DATA_DIR, the program is assuming that the raw data is in a subfolder called raw, and the reduced data will go in a subfolder called reduced.

For example the raw data is in raw/2020-11-09 and the reduced data will be in reduced/2020-11-09.

To run the reduction pipeline through astrometric correction:

```
reduce-breyo 2020-11-08 --preproc 
```

### If directory contains files that shouldn't be processed...

Move to reduced data, and move all bad files into a junk directory.  These files should be highlighted red in the observing log. (Eventually, we will do this while observing, so this step will become unnecessary.)


```
cd reduced/2020-11-08
mkdir junk
mv p-M32*.fits junk/.
mv p-danae*{1..3}r.fits junk/.
mv *RGB* junk/.
```

You might also need to remove some of the sky flats if the counts were too high/low.
```
cd calib
mkdir junk
```
You can also fix any filenames that were entered incorrectly while observing.  Again, this should be done at the telescope.

### Check image coordinates

At least once, the coordinates were not updated correctly in the image headers.
```
gethead OBJCTRA OBJCTDEC p*.fits
```
If there is a problem, run py/fixcoords.py. This can happen if you forget to connect MaxIm DL to the telescope, in which case MaxIm DL will not have any of the TCS information.

```
python ~/github/breyo/py/breyo/fixcoords.py --ra "20 20 53.24" --dec "+59 26 55.6"  --filestring p-TrES
```

### Other Gotchas

The bias frames should be named 'bias*.fits'.  If they are not, you need to rename them before running the pipeline.  For example:

```
rename 's/Bias/bias/' Bias*.fits

```

### Proceed with pipeline reduction with Kepler Data

move back to BREYO_DATA_DIR and proceed with pipeline.
```
cd ../../
```

If reducing data from the Kepler K4040 CMOS detector:
```
reduce-breyo 2020-11-08 --kepler --masterdarks --masterflats --reduceall --crzap --astrometry --refstars

```

If reducing data from the Kepler K4040 CMOS detector AND feeling lucky, try:
```
reduce-breyo 2020-11-08 --kepler --masterdarks --masterflats --reduceall --crzap --astrometry --refstars --photcalib --reproject

```

### If combining data from more than one night... 


```
reduce-breyo 2021-11-03 --preproc 
reduce-breyo 2021-11-04 --preproc
```

Then merge the files from the two nights (make sure names are not the same)
```
cd reduced/
cd 2021-11-03
cp ../2021-11-04/calib/dark*.fits calib/.
```


### If reducing data from the SBIG STL-11000M CCD

```
reduce-breyo 2020-11-08 --masterbias --masterflats --reduceall --crzap --astrometry --refstars

```

## Authors

* [**Rose Finn**](https://github.com/rfinn)
* [**John Moustakas**](https://github.com/moustakas)

## Acknowledgments

Thanks to the *stellar* students who helped write and debug this code including:
* Sandy Spicer
* Ryan Cirrincione
* Daniel Allspach

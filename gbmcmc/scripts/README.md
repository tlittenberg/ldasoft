# GBMCMC Post-processing scripts

# Table of contents
1. [Introduction](#intro)
2. [build_pandas_catalog](#build_pandas_catalog)
3. [GBcatalogDemo](#GBcatalogDemo)

<a name="intro"></a>
# Introduction
This directory contains useful scripts for post-processing GBMCMC data.


<a name="build_pandas_catalog"></a>
## build_pandas_catalog.py
The script `build_pandas_catalog.py` can be used to generate PANDAS dataframes containing GBMCMC outputs such as source catalogs, frequency residuals, MCMC chains, etc.  These data can then be manipulated using python with the Jupyter Notebook [`GBcatalogDemo.ipynb`](#GBcatalogDemo)] providing an example. 

### Description
The script primarily operates on the `catalog_N` output directories of one or more GBMCMC runs.  The identified sources in each of these catalogs are combined into a single PANDAS DataFrame with the key `detections` containing point-estimates of the source parameters.  Command-line options are used to provide info for the `metadata` DataFame, which is stored together with `detections` in the specified HDF5 output file. As an option, the script will also process the `data` directories of GBMCMC runs and extract frequency-series information such as FFT of data, FFT of model, FFT of noise model, etc. This information is placed in a DataFrame called `spectrum` in the same HDF5 file.  MCMC chain files for individual sources are also stored in auxilliary HDF5 files, with each HDF5 file limited to the chains from 100 GBMCMC runs for file size reasons. The files containing MCMC chains for each source are linked in the `detections` DataFrame for easy retrieval. 

The script can operate on either a set of local directories from one or more GBMCMC runs or on a set of directories stored on an AWS S3 file service. In both cases, the desired directories for catalog building are specified with a `pattern` which is matched using wildcard characters. In the cae for AWS S3, some catalog metadata are assumed to be encoded in the filenames following a convention used in the analysis of the Radler LDC data.

### Usage
Basic usage:

`build_pandas_catalog.py [-h] [-v] [-a] [-b BUCKET] [-n NAME] [-p PARENT] [-d DURATION] [-s SPECTRALPATTERN] [-N SAMPLES] pattern catalogFile`

#### pattern
this is the search pattern for local catalog directories or ZIP files in AWS S3. Wildcards such as * can be used to select multiple directories. e.g. `runs/mycat*` or `runs/newcats*.zip`.

#### catalogFile
this is the output file for the catalog. e.g. `mycatalog.h5`

#### DURATION
this is the duration of the data in seconds, e.g. `-d 12345678`

#### NAME
A name for the catalog. This can be used to track multiple catalogs, for example different durations.  It is independet of the catalog file name. e.g. `my_12mo_catalog`

#### PARENT
This is the name of the parent catalog, if present. This is used to track the evolution of sources between multiple catalogs corresponding to different observation durations of the same systems

#### SPECTRALPATTERN
Similar to pattern, this is used to search for the `data` directories of GBMCMC outputs in order to produce frequency-series data. If not set, frequency-series data is not processed

#### SAMPLES
The number of samples in the GBMCMC run, used to compute the overlap between adjacent segments when computing frequency-series data. NOTE: for AWS data, this will be recovered from the file name if not specified at the command line

#### AWS
Flag `-a` is used to set searches to AWS S3.  Also then assumes files are in ZIP format. 

#### BUCKET
For AWS searches, specifies the S3 bucket containing the data

### Examples

#### Process 6mo data in local directories
`python build_pandas_catalog.py -v -n local_test -d 15728640 -N 128 -s /home/Data/GBMCMC/seg /home/Data/GBMCMC/seg /home/Data/GBMCMC/local_test.h5`

#### Process 12mo data from AWS
`python build_pandas_catalog.py -a -b lisaprocdata -n aws_12mo_v0 -p aws_6mo_v2 -d 31457280 -v -s 31457280/31457280*.zip 31457280_catalogs/31457280*cat.zip /home/Data/cat31457280_v0.h5`

<a name="GBcatalogDemo"></a>
## GBcatalogDemo.ipynb
This is a Jupyter Notebook that provides an example of how to access GBMCMC data using Jupyter+PANDAS and make plots, population cuts, etc. It uses GBMCMC outputs processed by [`build_pandas_catalog.py`](#build_pandas_catalog).   The Notebook is documented internally.

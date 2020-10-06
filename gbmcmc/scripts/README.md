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
Basic usage
`build_pandas_catalog.py [-h] [-v] [-a] [-b BUCKET] [-n NAME] [-p PARENT] [-d DURATION] [-s SPECTRALPATTERN] [-N SAMPLES] pattern catalogFile`

<a name="GBcatalogDemo"></a>
## GBcatalogDemo.ipynb
This is a Jupyter Notebook that provides an example of how to access GBMCMC data using Jupyter+PANDAS and make plots, population cuts, etc. It uses GBMCMC outputs processed by [`build_pandas_catalog.py`](#build_pandas_catalog).   The Notebook is documented internally.

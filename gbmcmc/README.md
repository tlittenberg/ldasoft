# GBMCMC Manual

# Table of contents
1. [Introduction](#intro)
    1. [Dependencies](#dependencies)
    2. [Useful modules](#post_proc_dependencies)
    3. [Installation instructions](#installation)
2. [Example Use Cases](#examples)
    1. [High SNR injection](#highsnr)
    2. [LDC data](#ldc)
    3. [Verification binaries](#verification_binaries)
    4. [Other examples](#other)
3. [Description of Output](#output)
    1. [Description of source parameters](#parameters)
    2. [Contents of `chains` directory](#chains)
    3. [Contents of `data` directory](#data)
    4. [Contents of `run` directory](#run)
4. [Post Processing](#post)
    1. [Contents of `catalog_N` directory](#catalog_output)

<a name="intro"></a>
# Introduction

<a name="dependencies"></a>
## C Dependencies
```bash
gsl
gslcblas
openmp
mpi
hdf5
```
Note: For OpenMP on OSX, use `sudo port install libomp +top_level`

<a name="post_proc_dependencies"></a>
## Python packages for post-production
```bash
numpy
glob
pandas
astropy
matplotlib
chainConsumer
h5py
tables
pydot
<<<<<<< HEAD
lisacattools
```
<a name="installation"></a>
## Installation
Build and install binaries in `${HOME}/ldasoft/master/bin/` with `cmake`.
```
<a name="installation"></a>
## Installation
Build and install binaries in `${HOME}/ldasoft/master/bin/` with `cmake`
>>>>>>> 3848e2f8e40eaadda41c5966f615a685041c9b4e
```bash
./install.sh ${HOME}/ldasoft/master
```

<a name="examples"></a>
# Example use cases for GBMCMC
See run options with `gb_mcmc --help`.  Detailed description of command line options found [here](doc/running.md)

<a name="highsnr"></a>
## Analyze single high SNR injection
```bash
gb_mcmc \
  --inj <path to>/ldasoft/gbmcmc/etc/sources/precision/PrecisionSource_0.txt
```
<a name="ldc"></a>
## Analyze LDC data
See [Running on LDC Documentation](`README_LDC.md`)

<a name="verification_binaries"></a>
## Analyze verification binary with EM priors
Files containing verification binary parameters and priors are located in `ldasoft/gbmcmc/etc/sources/verification`.
Current best estimates of known binaries and their parameters, with references to source material, are found [here](https://docs.google.com/spreadsheets/d/1PfwgaPjOpcEz_8RIcf87doyhErnT0cYTYJ9b6fC0yfA/edit)

### Examples for AM CVn
Run with sky-location fixed
```bash
gb_mcmc \
  --inj <path to>/ldasoft/gbmcmc/etc/sources/verification/AMCVn.dat \
  --known-source \
  --no-rj \
  --cheat \
  --samples 128
```

To run with sky-location and EM-constrained priors (currently just U[cosi]) add
```bash
  --em-prior <path to>/ldasoft/gbmcmc/etc/sources/verification/AMCVn_prior.dat 
```

To run with orbital period P and dP/dt fixed add
```bash
  --fix-freq
```

<a name="other"></a>
## Other use cases

### Use UCBs as calibration sources
```bash
gb_mcmc \
  --inj <path to>/ldasoft/gbmcmc/etc/sources/calibration/CalibrationSource_0.txt \
  --duration 604948 \
  --segments 2 \
  --fit-gap \
  --samples 128 \
  --no-rj
```

<a name="output"></a>
# GBMCMC output format

<a name="parameters"></a>
## Parameterization

Each galactic binary is fit with `NP` parameters, where `NP` is typically 8, but can be increased to 9 using the `--f-double-dot`.  The parameters are:

* `f` GW frequency at mission start (Hz)
* `fdot` first time derivative of frequency (Hz/s)
* `A` GW amplitude
* `[long,cos_colat]` sky location parameters in solar system barycenter ecliptic coordinates.  The latitude parameter is cosine(co-latitude), where colatitude is PI/2 - latitude. We use co-latitude so that the cosine is not double-valued, and U[-1,1] maps to a uniform distribution on the sky.
* `cos_inc` cosine of inclination angle [-1,1]
* `psi` polarization angle [0,PI]
* `phi` phase of GW at mission start [0,2PI]
* **OPTIONAL** `f-double-dot` Second frequency derivative (Hz/s^2)

Converting from the internal sky location parameters to galactic coordinates or RA and dec is easy using `astropy`:

```python
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord

#AM CVn as example
gbmcmc_phi = 2.9737          
gbmcmc_costheta = 0.6080   

# source = SkyCoord(lon, lat, frame)
source = SkyCoord(gbmcmc_phi*u.rad,(np.pi/2 - np.arccos(gbmcmc_costheta))*u.rad,frame='barycentrictrueecliptic')

#convert to galactic coordiantes
print("\n Galactic coordinates:")
print(source.galactic)

#convert to RA & DEC coordiantes
print("\n RA & Dec (deg)")
print(source.icrs)

#convert to RA & DEC in h:m:s format
print("\n RA & Dec (hhmmss)")
print(source.to_string(style='hmsdms'))

#check the built in location of AM CVn
print("\n Built in location for AM CVn")
print(SkyCoord.from_name("AM CVn"))
```

<a name="chains"></a>
## chains directory

### Model 
`model_chain.dat.0`: Model state information for the chain sampling the target distribution (*temperature* T=0). Each row of the file is a chain sample. The columns are

    step | N_live | logL | logL_norm | t_0
    
where `step` is the iteration of the sampler, `N_live` is the number of templates used in the model, `logL` is the unnormalized log likelihood, i.e -0.5*Chi^2, `logL_norm` is the log of the likelihood normalization, and `t_0` is the inferred start time (in seconds) of the data segement.

Samples with `step < 0` are considered part of the burn in.  The value of `step` is occasionally reset when the sampler undergoes a sufficiently large change in the maximum log likelihood. 

The total log likelihood at each step is `logL + logL_norm`.  

If the analysis includes `M` time segments, there will be `M` `t_[m=0,-1]`  columns corresponding to the start time of each segment.
The `t_m` values are constant unless the analysis is actively fitting for segment start times with the `--fit-gap` flag.

If the run includes the `--verbose` flag, the file also contains `NP x N_live` additional columns cycling through the waveform parameters.

    f | fdot | A | cos_colat | long | cos_inc | psi | phi

Thus, when parsing the file, the `N_live` column should tell the parser how many additional parameter columns are appended to the end of the row. 

### Signal parameters
`parameter_chain.dat.0`: Full set of chain samples including burn-in. 
To properly parse the file into individual chain samples, use in concert with the
`model_chain.dat.M` files to extract the number of templates used at each step `N_live` .  For each chain iteration in the `parameter_chain.dat.0` file there are `N_live` lines with each containing the parameters for a single source.  

When running with the `--verbose` flag, there will be multiple `parameter_chain.dat.M` files with the index `M` for referring to the *temperature* of the parallel-tempered chain that is producing these outputs. `M=0` is the target distribution; `M>0` is softened versions of this
distribution. If you don't know what this means, you want `M=0`.  

NOTE: These files contain *all* samples including
burnin.  Unless you care deeply about checkpointing or run diagnostics, you do
not need these files. Use `dimension_chain.dat.N` files, as described below.

### Posterior samples for specfic models
`dimension_chain.dat.N`:  The post-burnin posterior samples of
the binary parameters for the subset of models that have exactly `N` binaries.
The columns are

    f | fdot | A | cos_colat | long | cos_inc | psi | phi

If there are more than one source, the sources cycle in the output (so line 1 is
the first source, first posterior sample; line 2 is the second source, first
posterior sample, etc).

### Noise parameters
`noise_chain.dat.0`: Full set of chain samples for noise model parameters. The columns are

    step | logL | logL_norm | etaA | eta E
    
where `etaY` is the PSD parameter for channel `Y`, which is a constant multiplyer to the predicted noise PSD: `Sn_model(f) = eta * Sn_theory(f)`.

For 4-channel data there are only 4 columns and `X` replaces `A`.  

When the `--quiet` flag is in use only the post-burnin samples are recorded.

### Parallel tempering likelihoods
`log_likelihood_chain.dat`: The full log likelihood (`logL + logLnorm`) for each of the `NC` parallel tempered chains, at each chain step. Columns are

    step | logL_0 | logL_1 | ... | logL_NC-1

Useful for monitoring the effectiveness of the parallel tempering by checking that adjacent chains are exchanging likelihood values.

The `logL_0` column is the likelihood for the target distribution. 
**File is not filled when `--quiet` flag is in use.**

### Chain temperatures
`temperature_chain.dat`: Inverse temperature for each of the `NC` parallel tempered chains, at each chain step. Columns are

    step | 1/T_0 | 1/T_1 | ... | 1/T_NC-1

Useful for monitoring the adaptive temperature spacing.
**File is not filled when `--quiet` flag is in use.**

<a name="data"></a>
## data directory
`data_0_0.dat`: Fourier series of input data. For 6-link data columns are 

    f [Hz] | Re{A(f)} | Im{A(f)} | Re{E(f)} | Im{E(f)} 
    
For 4-link data there are two columns with the X channel replacing A.

The `0_0` in the filename indicate that the data is for the first (0th) time and frequency segment. If the analysis involves multiple time or frequency segments, different `data_i_j.dat` files will be filled.

`waveform_injection_0_0.dat`: Fourier series of the injected signals. For 6-link data columns are 

    f [Hz] | Re{h_A(f)} | Im{h_A(f)} | Re{h_E(f)} | Im{h_E(f)} 
    
For 4-link data there are two columns with the X channel replacing A.

The same `i_j.dat` naming convention as the `data_i_j.dat` files is used here.
In the absence of an injection for the analysis, the file contains a copy of the data.

`power_data_0_0.dat`: Power spectrum of input data. For 6-link data columns are 

    f [Hz] | A^2(f) | E^2(f) 
    
For 4-link data there are two columns with the X channel replacing A.

The same `i_j.dat` naming convention as the `data_i_j.dat` files is used here.

`power_injection_0_0.dat`: Power spectrum of injected signals. For 6-link data columns are 

    f [Hz] | h_A^2(f) | h_E^2(f) 
    
For 4-link data there are two columns with the X channel replacing A.

The same `i_j.dat` naming convention as the `data_i_j.dat` files is used here.
If no injection is being done, the file contains the power spectrum of the data and is identical to `power_data_0_0.dat`.


`frequency_proposal.dat`: Smoothed and normalized power spectrum of data used for frequency proposal. Columns are 

    frequency bin number | probability density
    
**THIS FILE COULD BE REMOVED FROM OUTPUT OR PUT IN `--verbose` OUTPUT**

`power_noise_t0_f0.dat`: Quantiles of the posterior distribution for the reconstructed noise model. For 6-link data the columns are

    f | SnA_50 | SnA_25 | SnA_75 | SnA_05 | SnA_95 | SnE_50 | SnE_25 | SnE_75 | SnE_05 | SnE_95

Where `SnY` is the noise power spectral density for channel `Y`, `X_NN` is the NN% quantile of `X`, e.g. `[SnA_25--SnA_75]` encompasses the 50% credible interval for `SnA`. For 4-link data there are 6 columns and X replaces A.

If multiple time/frequency segments are being analyzed there will be `power_noise_t[i]_f[j].dat` files for each time-frequency segment.

`power_reconstruction_t0_f0.dat`: Quantiles of the posterior distribution for the reconstructed signal model. For 6-link data the columns are

    f | hA^2_50 | hA^2_25 | hA^2_75 | hA^2_05 | hA^2_95 | hE^2_50 | hE^2_25 | hE^2_75 | hE^2_05 | hE^2_95

Where `hY^2` is the signal power spectral density for channel `Y`, `X_NN` is the NN% quantile of `X`, e.g. `[hA^2_25--hA^2_75]` encompasses the 50% credible interval for `hY`. For 4-link data there are 6 columns and X replaces A.

If multiple time/frequency segments are being analyzed there will be `power_reconstruction_t[i]_f[j].dat` files for each time-frequency segment.

`power_residual_t0_f0.dat`: Quantiles of the posterior distribution for the data residual. For 6-link data the columns are

    f | rA^2_50 | rA^2_25 | rA^2_75 | rA^2_05 | rA^2_95 | rE^2_50 | rE^2_25 | rE^2_75 | rE^2_05 | rE^2_95

Where `rY^2` is the residual power spectral density for channel `Y`, `X_NN` is the NN% quantile of `X`, e.g. `[rA^2_25--rA^2_75]` encompasses the 50% credible interval for `rA^2`. For 4-link data there are 6 columns and X replaces A.

If multiple time/frequency segments are being analyzed there will be `power_residual_t[i]_f[j].dat` files for each time-frequency segment.

`variance_residual_t0_f0.dat`: Variance of the signal model amplitude (and therefore residual) computed by summing the variance of the real and imaginary amplitudes of the joint signal model, i.e. Var(h) = Var(Re h) + Var(Im h). For 6-link data the columns are

    f | Var(hA) | Var(hE)
    
For 4-link data there are only two columns and X replaces A. 

The intent here is that this variance can be added to the noise PSD for a follow-on analysis of the residuals to take into account the uncertainty in the signal model used in the subtraction. **This needs to be tested**

`waveform_draw_0.dat`: Fair draw of the data, signal model, and residual printed periodically throughout MCMC analysis to check on health of the fit. For 6-link data the columns are 

    f | dA^2 | dE^2 | hA^2 | hE^2 | rA^2 | rE^2
    
where `dY`,`hY`, and `rY` is the data, waveform, and residual of channel `Y`. 
For 4-link data there are four columns and X replaces A.

<a name="run"></a>
## main run directory
`evidence.dat`: Posterior for number of templates used to fit the data. Columns are

    M | NN
    
Where `M` is the model "dimension", i.e. the number of templates in the model, and `NN` signifies that the sampler spent NN% of iterations in that model.

`gb_mcmc.log`: Run information including git version number, command line, summary of data settings and run flags, and  total run time.

`run.sh`: Bash script which echos command line used for analysis.

`avg_log_likelihood.dat`:  Average log likelihoood for each parallel tempered chain. These data would serve as the integrand for a thermodynamic integration step to compute the overall evidence for the model (marginalized over the number of galactic binary signals in the data). Columns are

    1/T | <logL+logL_norm>

The first row is the target distribution with `T=1`

`example_gb_catalog.sh`: Example `bash` script for post-processing with `gb_catalog`. Takes as argument the size of the model you want to post-process (e.g., maximum from `evidence.dat`).
<a name="post"></a>
# Post processing GBMCMC

Post production of Markov chain output from `gb_mcmc` is done with `gb_catalog`.  
The post processing is done on a single model from the RJMCMC analysis, identified by the number of sources `N` used in the fit, which the user specifies.
Results from `gb_catalog` are found in the a subdirectory of the main run directory called `catalog_N/`.  The catalog step can be run several times for different choices of `N`, creating different directories for each run.

For run options with `gb_catalog --help`.  Detailed description of command line options (will eventually be) found [here](doc/running.md).

An example `gb_catalog` run script is created by `gb_mcmc` in the run directory.  See `example_gb_catalog.sh`.

<a name="catalog_output"></a>
## GB Catalog output format
Within `catalog_N/` are post production files for each candidate source found by the analysis as well as aggregate files describing the whole catalog.
The source files are uniquely named by the gravitational wave frequency of the binary as `LDCFFFFFFFFFF`. 
`LDC` is the prefix for the name to make clear that production-level analysis have been using LISA Data Challenge simulations.
The `FFFFFFFFFF` is a 10-digit representation of the **median** gravitational wave frequency, in Hz.  
For example, a binary at `f_0 = 3 mHz` would be named  `LDC0030000000`.

**NOTE:** There are not necessarily `N` sources found in the `catalog_N` directory.  Additional filtering is done on the chain samples which may remove sources in the "padding" regions of the analysis window, or marginal detection candidates that did not have sufficient evidence to be included in the catalog.

Contents of each indidvidual file are as follows:

`catalog_N/LDCFFFFFFFFFF_chain.dat`: The Markov chain samples for parameters consistent with `LDCFFFFFFFFFF`. From these samples the marginalized posteriors for the source are computed. If correlations between sources are of interest the user has to go back to the full chains. Columns are the same as the raw chain files from `gb_mcmc` with **three** additional columns appended

    f | fdot | A | cos_colat | long | cos_inc | psi | phi | SNR | Waveform Match | Waveform Distance
    
Where `SNR` is computed relative to the noise level input using `--noise-file` or the analytic noise functions if not specified, `Waveform Match` is the match between the current sample and the reference waveform for the catalog _Entry_, used for filtering the raw chains, and `Waveform Distance` is the waveform distance between the current and reference waveforms--perhaps a future filtering criteria.

`catalog_N/LDCFFFFFFFFFF_params.dat`: Point estimate of the candidate source. Parameters correspond to the chain sample containing the **median** of the 1D marginalized distriubtion for the GW frequency

    f | fdot | A | cos_colat | long | cos_inc | psi | phi

`catalog_N/LDCFFFFFFFFFF_power_reconstruction.dat`: Reconstructed waveform represented by quantiles of the posterior distribution of the signal's power spectrum. For 6-link data the columns are


    f | hA^2_50 | hA^2_25 | hA^2_75 | hA^2_05 | hA^2_95 | hE^2_50 | hE^2_25 | hE^2_75 | hE^2_05 | hE^2_95
    
Where `hY^2` is the signal power spectral density for channel `Y`, `X_NN` is the NN% quantile of `X`, e.g. `[hA^2_25--hA^2_75]` encompasses the 50% credible interval for `hY`. For 4-link data there are 6 columns and X replaces A.

`catalog_N/LDCFFFFFFFFFF_waveform.dat`: Point estimate of the reconstructed waveform corresponding  to the chain sample containing the **median** of the 1D marginalized distriubtion for the GW frequency. Columns for 6-link data are

    f [Hz] | Re{h_A(f)} | Im{h_A(f)} | Re{h_E(f)} | Im{h_E(f)} 
    
For 4-link data there are two columns with the X channel replacing A.

`catalog_N/entries.dat`: Summary of sources found in catalog.  Columns are

    LDCFFFFFFFFFF | SNR | Evidence
    
Where `SNR` is computed relative to the noise level input using `--noise-file` or the analytic noise functions if not specified.  The evidence, with values from [0,1], is the probability that the candidate source is a member of the catalog, computed by the fraction of chain samples that include a set of parameters consistent with `LDCFFFFFFFFFF`.

`catalog_N/history.dat`: File associating current candidate sources to previous catalog input using `--catalog` and `--Tcatalog` arguments for `gb_catalog`.  The columns are

    LDCMMMMMMMMMM from input catalog | LDCNNNNNNNNNN from current catalog
    
Only sources in the current catalog with clearly identified _parent_ sources from the previous analysis are included in the file.


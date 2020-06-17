# GBMCMC Manual
## Dependencies
```bash
gcc
gsl
gslcblas
```

## Good things to have for post-production
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
```
## Installation
Build and install binaries in `${HOME}/ldasoft/master/bin/` 
```bash
cd <path to>/ldasoft/gbmcmc/src
make
make install
```
Edit `ldasoft/gbmcmc/src/Makefile` to change install destination.

# Example use cases for GBMCMC
See run options with `gb_mcmc --help`

## Analyze on-the-fly data & source simulation
```bash
gb_mcmc \
  --inj <path to>/ldasoft/gbmcmc/etc/sources/precision/PrecisionSource_0.txt
```

## Analyze LDC data
See `README_LDC.md`

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


# GBMCMC output format

## Parameterization

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

## `chains` directory

### `dimension_chain.dat.N`
The `dimension_chain.dat.N` files contain the post-burnin posterior samples of
the binary parameters for the subset of models that have exactly `N` binaries.
The columns are

    f | fdot | A | cos_colat | long | cos_inc | psi | phi

If there are more than one source, the sources cycle in the output (so line 1 is
the first source, first posterior sample; line 2 is the second source, first
posterior sample, etc).

### `model_chain.dat.0`

### `parameter_chain.dat.0`
The `parameter_chain.dat.M` files should be used in concert with the
`model_chain.dat.M` files.  `model_chain.dat.M` records how many binaries were
in the model at each MCMC step; there are a corresponding number of lines in
`parameter_chain.dat.M` for that step's binary parameters (i.e.
`parameter_chain.dat.M` contains a mashed-together compilation of the
`dimension_chain.dat.N` files).  The index `M` for these files refers to the
*temperature* of the parallel-tempered chain that is producing these outputs (`M
= 0` is the target distribution; `M > 0` is softened versions of this
distribution; if you don't know what this means, you want `M=0`).  Unlike
`dimension_chain.dat.N` files, these files contain *all* samples including
burnin.  Unless you care deeply about checkpointing or run diagnostics, you do
not need these files. Use `dimension_chain.dat.N` files, as described above.

### `log_likelihood_chain.dat`
### `noise_chain.dat.0`
### `temperature_chain.dat`

## `data` directory
### `data_0_0.dat`
Fourier series of input data. For 6-link data columns are 

    f [Hz] | Re{A(f)} | Im{A(f)} | Re{E(f)} | Im{E(f)} 
    
For 4-link data there are two columns with the X channel replacing A.

The `0_0` in the filename indicate that the data is for the first (0th) time and frequency segment. If the analysis involves multiple time or frequency segments, different `data_i_j.dat` files will be filled.

### `waveform_injection_0_0.dat`
Fourier series of the injected signals. For 6-link data columns are 

    f [Hz] | Re{h_A(f)} | Im{h_A(f)} | Re{h_E(f)} | Im{h_E(f)} 
    
For 4-link data there are two columns with the X channel replacing A.

The same `i_j.dat` naming convention as the `data_i_j.dat` files is used here.
In the absence of an injection for the analysis, the file contains a copy of the data.

### `power_data_0_0.dat`
Power spectrum of input data. For 6-link data columns are 

    f [Hz] | A^2(f) | E^2(f) 
    
For 4-link data there are two columns with the X channel replacing A.

The same `i_j.dat` naming convention as the `data_i_j.dat` files is used here.

### `power_injection_0_0.dat`
Power spectrum of injected signals. For 6-link data columns are 

    f [Hz] | h_A^2(f) | h_E^2(f) 
    
For 4-link data there are two columns with the X channel replacing A.

The same `i_j.dat` naming convention as the `data_i_j.dat` files is used here.
If no injection is being done, the file contains the power spectrum of the data and is identical to `power_data_0_0.dat`.


### `frequency_proposal.dat`
Smoothed and normalized power spectrum of data used for frequency proposal. Columns are 

    frequency bin number | probability density
    
**THIS FILE COULD BE REMOVED FROM OUTPUT OR PUT IN `--verbose` OUTPUT**

### `power_noise_t0_f0.dat`
Quantiles of the posterior distribution for the reconstructed noise model. For 6-link data the columns are

    f | SnA_50 | SnA_25 | SnA_75 | SnA_05 | SnA_95 | SnE_50 | SnE_25 | SnE_75 | SnE_05 | SnE_95

Where `SnY` is the noise power spectral density for channel `Y`, `X_NN` is the NN% quantile of `X`, e.g. `[SnA_25--SnA_75]` encompasses the 50% credible interval for `SnA`. For 4-link data there are 6 columns and X replaces A.

If multiple time/frequency segments are being analyzed there will be `power_noise_t[i]_f[j].dat` files for each time-frequency segment.

### `power_reconstruction_t0_f0.dat`
Quantiles of the posterior distribution for the reconstructed signal model. For 6-link data the columns are

    f | hA^2_50 | hA^2_25 | hA^2_75 | hA^2_05 | hA^2_95 | hE^2_50 | hE^2_25 | hE^2_75 | hE^2_05 | hE^2_95

Where `hY^2` is the signal power spectral density for channel `Y`, `X_NN` is the NN% quantile of `X`, e.g. `[hA^2_25--hA^2_75]` encompasses the 50% credible interval for `hY`. For 4-link data there are 6 columns and X replaces A.

If multiple time/frequency segments are being analyzed there will be `power_reconstruction_t[i]_f[j].dat` files for each time-frequency segment.

### `power_residual_t0_f0.dat`
Quantiles of the posterior distribution for the data residual. For 6-link data the columns are

    f | rA^2_50 | rA^2_25 | rA^2_75 | rA^2_05 | rA^2_95 | rE^2_50 | rE^2_25 | rE^2_75 | rE^2_05 | rE^2_95

Where `rY^2` is the residual power spectral density for channel `Y`, `X_NN` is the NN% quantile of `X`, e.g. `[rA^2_25--rA^2_75]` encompasses the 50% credible interval for `rA^2`. For 4-link data there are 6 columns and X replaces A.

If multiple time/frequency segments are being analyzed there will be `power_residual_t[i]_f[j].dat` files for each time-frequency segment.

### `variance_residual_t0_f0.dat`
Variance of the signal model amplitude (and therefore residual) computed by summing the variance of the real and imaginary amplitudes of the joint signal model, i.e. Var(h) = Var(Re h) + Var(Im h). For 6-link data the columns are

    f | Var(hA) | Var(hE)
    
For 4-link data there are only two columns and X replaces A. 

The intent here is that this variance can be added to the noise PSD for a follow-on analysis of the residuals to take into account the uncertainty in the signal model used in the subtraction. **This needs to be tested**

### `waveform_draw_0.dat`
Fair draw of the data, signal model, and residual printed periodically throughout MCMC analysis to check on health of the fit. For 6-link data the columns are 

    f | dA^2 | dE^2 | hA^2 | hE^2 | rA^2 | rE^2
    
where `dY`,`hY`, and `rY` is the data, waveform, and residual of channel `Y`. 
For 4-link data there are four columns and X replaces A.

## main run directory
### `evidence.dat`
Posterior for number of templates used to fit the data. Columns are

    M | NN
    
Where `M` is the model "dimension", i.e. the number of templates in the model, and `NN` signifies that the sampler spent NN% of iterations in that model.

### `gb_mcmc.log`
Run information including git version number, command line, summary of data settings and run flags, and  total run time.

### `run.sh`
Bash script which echos command line used for analysis.

### `avg_log_likelihood.dat`

# Post processing GBMCMC

To be continued...

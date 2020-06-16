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

## `chains` directory

### `dimension_chain.dat.N`
The `dimension_chain.dat.N` files contain the post-burnin posterior samples of
the binary parameters for the subset of models that have exactly `N` binaries.
The columns are

    f fdot A cos_colat long cos_inc psi phi

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
### `frequency_proposal.dat`
### `power_data_0_0.dat`
### `power_injection_0_0.dat`
### `power_noise_t0_f0.dat`
### `power_reconstruction_t0_f0.dat`
### `power_residual_t0_f0.dat`
### `variance_residual_t0_f0.dat`
### `waveform_draw_0.dat`
### `waveform_injection_0_0.dat`

## main run directory
### `evidence.dat`
### `gb_mcmc.log`
### `run.sh`
### `avg_log_likelihood.dat`

# Post processing GBMCMC

To be continued...

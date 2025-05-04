\mainpage

Prototype data analysis software for LISA analysis

# Doxygen code documentation
https://tlittenberg.github.io/ldasoft/html/index.html

# Acknowledgment

```bibtex
@misc{ldasoft,
	author = "Tyson B. Littenberg, Neil J. Cornish, Kristen Lackeos, Travis Robson",
	title = "LDASoft",
	howpublished = "free software (GPL)"
	doi = "10.5281/zenodo.2026177"
	year = "2020"
}
```

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2026177.svg)](https://doi.org/10.5281/zenodo.2026177)

# C Dependencies 
```bash
gsl
gslcblas
openmp
mpi
hdf5

```
The `GlobalFit` pipeline additionally requires the `global-fit` branch of [`LISA-Massive-Black-Hole`](https://github.com/eXtremeGravityInstitute/LISA-Massive-Black-Hole/tree/global-fit) to be installed.

# Installation

Example build script:
```bash
#!/bin/bash

# set prefix for install directories
export LDASOFT_PREFIX=<EDIT: path to desired install directories>

# build codes
./install.sh ${LDASOFT_PREFIX}

# add location of binaries to PATH 
export PATH=${LDASOFT_PREFIX}/bin:${PATH}
```

Directory tree for the ldasoft installation:
```
├── bin
│   ├── gaussian_mixture_model
│   ├── global_fit
│   ├── noise_mcmc
│   ├── noise_spline_mcmc
│   ├── ucb_catalog
│   ├── ucb_chirpmass
│   ├── ucb_grid
│   ├── ucb_match
│   ├── ucb_mcmc
│   └── vgb_mcmc
├── include
│   ├── gitversion.h
│   ├── glass_constants.h
│   ├── glass_data.h
│   ├── glass_gmm.h
│   ├── glass_lisa.h
│   ├── glass_math.h
│   ├── glass_noise.h
│   ├── glass_noise_io.h
│   ├── glass_noise_model.h
│   ├── glass_noise_sampler.h
│   ├── glass_ucb.h
│   ├── glass_ucb_catalog.h
│   ├── glass_ucb_data.h
│   ├── glass_ucb_fstatistic.h
│   ├── glass_ucb_io.h
│   ├── glass_ucb_model.h
│   ├── glass_ucb_prior.h
│   ├── glass_ucb_proposal.h
│   ├── glass_ucb_sampler.h
│   ├── glass_ucb_waveform.h
│   └── glass_utils.h
└── lib
    ├── libglass_noise.a
    ├── libglass_ucb.a
    └── libglass_utils.a

```


# Issue tracker
https://github.com/tlittenberg/ldasoft/issues

# Other resources
 + [Example LISA analysis flow chart](https://www.draw.io/#Htlittenberg%2Fldasoft%2Fmaster%2FLISADataFlow.drawio)
 + [Global Fit block diagram](https://app.diagrams.net/#Htlittenberg%2Fldasoft%2Fmaster%2FGlobalFit.drawio)

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
│   ├── Bright_Remove
│   ├── Confusion_Fit
│   ├── Fisher_Galaxy
│   ├── Full_Residual
│   ├── Galaxy
│   ├── OrbitFile
│   ├── Setup
│   ├── gaussian_mixture_model
│   ├── gb_catalog
│   ├── gb_match
│   ├── gb_mcmc
│   ├── gb_mcmc_brans_dicke
│   ├── gb_mcmc_chirpmass
│   ├── gb_residual
│   └── global_fit
├── include
│   ├── Constants.h
│   ├── GalacticBinary.h
│   ├── GalacticBinaryCatalog.h
│   ├── GalacticBinaryData.h
│   ├── GalacticBinaryFStatistic.h
│   ├── GalacticBinaryIO.h
│   ├── GalacticBinaryMCMC.h
│   ├── GalacticBinaryMath.h
│   ├── GalacticBinaryModel.h
│   ├── GalacticBinaryPrior.h
│   ├── GalacticBinaryProposal.h
│   ├── GalacticBinaryWaveform.h
│   ├── LISA.h
│   └── gitversion.h
└── lib
    ├── libgbmcmc.a
    ├── liblisa.a
    └── libtools.a
```


# Issue tracker
https://github.com/tlittenberg/ldasoft/issues

# Other resources
 + [Example LISA analysis flow chart](https://www.draw.io/#Htlittenberg%2Fldasoft%2Fmaster%2FLISADataFlow.drawio)
 + [Global Fit block diagram](https://app.diagrams.net/#Htlittenberg%2Fldasoft%2Fmaster%2FGlobalFit.drawio)

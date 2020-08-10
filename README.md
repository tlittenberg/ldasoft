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

# Installation
```bash
#!/bin/sh

# make install directories
mkdir -p ${HOME}/ldasoft/master/ #hard-coded installation path in Makefiles (for now)

# build codes
./install.sh

# add location of binaries to PATH (include this in login scripts)
export PATH=${HOME}/ldasoft/master/bin/:$PATH
```

## cmake
```bash
#!/bin/sh

mkdir build
cd build/
# choose where you want to install it
cmake .. -DCMAKE_INSTALL_PREFIX=${CONDA_PREFIX}
make
make install
```


# Issue tracker
https://github.com/tlittenberg/ldasoft/issues

# Other resources
 + [Example LISA analysis flow chart](https://www.draw.io/#Htlittenberg%2Fldasoft%2Fmaster%2FLISADataFlow.drawio)

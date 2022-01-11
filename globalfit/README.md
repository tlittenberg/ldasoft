# GlobalFit Manual

<a name="intro"></a>
# Introduction

The `global_fit` application is a wrapper that uses a blocked Gibbs algorithm cycling through the `gb_mcmc`, `noise_mcmc`, and `vb_mcmc` modules to sample the joint posterior. 
Communication between modules is handled with `MPI`. 
A schematic block diagram showing how the processes are divided, and how the communication between processes is managed, is found [here](https://app.diagrams.net/#Htlittenberg%2Fldasoft%2Fmaster%2FGlobalFit.drawio).

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

<a name="installation"></a>
## Installation
Build and install as part of the full `ldasoft` package.  Example for install location `${HOME}/ldasoft/master/bin/` with `cmake`.
```bash
./install.sh ${HOME}/ldasoft/master
```

<a name="examples"></a>
# Example use cases for GlobalFit
See run options with `global_fit --help`.  

The CLI is identical to `gb_mcmc` (and is actually just calling the `gb_mcmc` parsing functions which are wrapped around `getopt`). 
A detailed description of tge command line options is found [here](../gbmcmc/doc/running.md)

The `global_fit` app requires at least three `MPI` processes (one for each supported module in the model, `noise`, `ucb`, and `vgb`).

**NOTE:** There is currently only support for running on data with the same format as the LISA Data Challenge's [*Sangria* dataset](https://lisa-ldc.lal.in2p3.fr/challenge2).

To analyze two adjacent frequency bands of a 3-month segment from the LDC *Sangria* training data, using a common spline noise model across the segments, starting at a frequency of 6 mHz, use
```bash
mpirun -np 4 global_fit \
  --h5-data /path/to/LDC2_sangria_training_v1.h5 --sangria \
  --fmin 0.006 \
  --duration 7864320 \
  --samples 128 --padding 16
```
In this example the verification binary module `vb_mcmc` is not activated, so process 1 will be idle. Process 0 will contain the noise model, and processes 2 and 3 are each responsible for the UCB analysis of a narrow-band segment.
The `padding` argument is the size (in frequency bins) of the overlap between the two UCB segments to prevent artifacts in the analysis from edge-effects.

The application sets up the following tree in the run directory for storing the output data
```
.
├── mbh
│   ├── src0000
│   ├── src0001
├── noise
│   ├── chains
│   ├── checkpoint
│   └── data
├── ucb
│   ├── seg_0000
│   │   ├── chains
│   │   ├── checkpoint
│   │   └── data
│   └── seg_0001
│       ├── chains
│       ├── checkpoint
│       └── data
└── vgb
    └── seg_0000
        ├── chains
        ├── checkpoint
        └── data
```

Because `global_fit` is just a wrapper around the separately developed (and documented) sampelers, the output data is of the same format as described in the documentation for `gb_mcmc`, `noise_mcmc`, and `vb_mcmc`.
Similarly, the post-processing is all handled on a segment-by-segment basis using the same tools as for the stand-alone applications of the different samplers.

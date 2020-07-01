# Detailed Breakdown of GBMCMC Command Line Arguments

## REQUIRED

## OPTIONAL

### Run time flags 
`[-h | --help]`: Print help message and exit

`[-v | --verbose]`: Print exhaustive details to file during run, including hot chains. Increases frequency of `stdout` output and flushes file buffers after each chain sample. **Use sparingly and only during debugging**. 

`[-q | --quiet]`: Print restricted details to file during run. Sufficient output is printed to fully assess results and to resume interrupted runs. Limited `stdout` and diagnostic information. **Use when running large batches**.

`[-d | --debug]`: Use coarse grids for creating data-driven proposals (to speed up initialization time) and automatically envokes `--verbose` setting.

### LISA configuration flags

`[--orbit=FILENAME]`: Specify orbit ephemerides file. When not specified analytic Eccentric-Inclined orbits with average armlength of 2.5 Gm are used. 

`[--links=INT; default=6]`: Specify number of optical links in constellation [4=single interferometer channel `X`,6=two interferometer channels `A` and `E`].   

`[--frac-freq]`: Change units of data from phase (displacement) to fractional frequency data (velocity). **Use when analyzing LDC data**   

### Data settings and flags
`[--data=FILENAME]`: Data file containing frequency-domain strain data. For 6-link data the columns are:

    f | ReA(f) | ImA(f) | ReE(f) | ImE(f)

For 4-link data there are only 3 columns with `X` replacing `A`. When not specified data will be simulated internally.

`[--samples=INT; default=2048]`: Number of frequency bins in analysis segment. Consider expected bandwidth of signals when setting.    

`[--padding=INT; default=0]`: Number of frequency bins added to either side of data segment to prevent edge effects when analyzing full bandwidth data. Full data segment is `samples+2*padding`.
       
`[--segments=INT; default=1]`: Number of adjacent time-frequency segments to be analyzed simultaneously. Frequency band will be the same for each segment. Time difference between segments is set by `--gap-time`. The galactic binaries are coherent between the time segments but the noise model is independent in the different segments.         

`[--start-time=FLOAT; default=0.0 (s)]`: Start time of first data segment (used if particular orbital initial conditions are required).

`[--gap-time=FLOAT; default=0.0 (s)]`: Duration of data gaps between adjacent `segments`.

`[--fmin=DOUBLE; default=0.0 (Hz)]`: Minimum frequency of data `segment` before `padding` is applied. When used with `--inj` flag `fmin` is set to be `segment/2` frequency bins before the central frequency of the simulated source.

`[--duration=FLOAT; default=62914560.0 (s)]`: Observation time for data segment. Default value corresponds to 2 years. **Changes to `duration` should also include changes to `samples`**.

`[--sim-noise]`: Add a simulated noise realization to data segment. **Use when not specifying `--data` file unless zero-noise run is desired**.

`[--conf-noise]`: Include analytic model for confusion noise in noise power sepctral density.
 
`[--noiseseed=INT; default=151226]`: Seed for random number generator that produces noise realization.                  
  
`[--inj=FILENAME]`: File containing parameters of single source to be injected into data. Columns are
 
    initial frequency (Hz) | frequency derivative (Hz s^-1) | ecliptic co-latitude (rad) | ecliptic longitude (rad) | amplitude | inclination (rad) | polarization angle (rad) | initial phase (rad)
    
Up to `10` injection files can be passed at once using multiple instances of `--inj`.

`[--injseed=INT; default=151012`: Seed for random number generator used for injection parameters that are drawn internally when ingored from `--inj` file because of `--known-source` flag (e.g. when injecting verification binaries).       

`[--catalog=FILENAME]`: Parameter file with each row corresponding to an individual source, and columns the same as for `--inj` file.  Sources in the `catalog` file that are outside of the analysis region, but have bandwidth that extends into the region, are forward modeled and subtracted from the input data to prevent edge effects.                

### Markov chain settings 
`[--steps=INT; default=100000]`: Number of "burn in" AND post-burn in" MCMC iterations.  Burn-in counter is reset when large likelihood jumps are recorded in the chain, so typical runs will produce more than `2 * steps` iterations in total.       

`[--chainseed=INT; default=150914]` Seed for random number generator used in MCMC, including initializing chain state, proposal distributions, and acceptance probabilities.     

`[--chains=INT; default=12]`: Number of chains used for parallel tempering MCMC.

`[--no-burnin]`: Skip burn-in steps and assume that every chain sample is a fair draw from the posterior. **Use when only interested in parameter estimation studies of individual sources**.                  

### Signal model settings

`[--sources=INT; default=10]`: Maximum number of sources allowed in model.        

`[--cheat]`: When used with `--inj` flag, starts chain at injection parameters.  Use with `--no-burnin` if only interested in parameter estimation. 
       
`[--f-double-dot]`: Include second time derivative of frequency in signal model.  The `d^2f/dt^2` parameter will appear in an additional column appended to end of usual output.      

`[--prior]`: Replaces likelihood function with a constant to generate sample from drawn from the prior. **Use when testing detailed balance of sampler**.                   

`[--no-rj]`: Disable trans-dimensional MCMC moves resulting in a fixed dimension signal model. **FIXME: ONLY SUPPORTS SINGLE SOURCE MODEL**                
       
`[--calibration]`: Marginalize over phase and amplitudfe calibration errors. **FIXME: UNTESTED** 

`[--fit-gap]`: Use signal model coherence across segments to fit for gap duration.  Used when `segments>1` to solve for inter-segment gaps. Start time of observations is fixed. 

### Priors & Proposals
       --fix-sky     : pin sky params to injection         
       --galaxy-prior: use galaxy model for sky prior      
       --snr-prior   : use SNR-based amplitude prior       
       --em-prior    : update prior ranges from other obs  
       --known-source: injection is VB (draw orientation)  
       --detached    : detached binary(i.e., use Mc prior) 
       --update      : use chain as proposal [filename]    
       --update-cov  : use cov mtrx proposal [filename]    

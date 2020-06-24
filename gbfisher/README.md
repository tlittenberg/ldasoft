# FisherGalaxy Manual

## Input file format 
Space-separated ASCII file, one source per row. Columns are <a name="params"></a>

    f (Hz) |  df/dt (s^-2) | co-latitude (rad) | longitude (rad) | amplitude | inclination (rad) | polarization angle (rad) | ref phase (rad)`
    
The sky location parameters are in solar system barycenter ecliptic coordinates.  
We use co-latitude (pi/2 - latitude) so that its cosine is distributed from \[0,1\]
Converting to these coordinates from galactic coordinates or RA and Dec can be done with `astropy`:

```python
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord

##############################
# Using AM CVn as an example #
##############################

# Convert from galactic to ecliptic coordinate
amcvn_galactic = SkyCoord(140.23427349*u.deg, 78.93818425*u.deg, frame='galactic')
amcvn_ecliptic = amcvn_galactic.transform_to('barycentrictrueecliptic')

# Convert from RA & Dec to ecliptic coordinate
amcvn_RAandDec = SkyCoord(188.72759585*u.deg, 37.62892244*u.deg, frame='icrs')
amcvn_ecliptic = amcvn_RAandDec.transform_to('barycentrictrueecliptic')

# Extract sky angles and convert to colatitude
colatitude = np.pi/2 - amcvn_ecliptic.lat.radian
longitude = amcvn_ecliptic.lon.radian
```

## Create orbit file with spacecraft ephimeredes
Compute spacecraft coordinates during the observation time. This is a separate step to easily test different orbit models if numerically produced ephemerides are available. The process here uses the "Eccentric Inclined" analytic orbits similar to the LISA Data Challenges.

`>OrbitFile`

Output files:

`EccentricInclined.txt`: Cartesian coordinates x(t),y(t),z(t) in solar system barycentric coordinates of each spacecraft. Columns are

    t (s) | SC1_x (m) | SC1_y (m) | SC1_z (m) | SC2_x (m) | SC2_y (m) | SC2_z (m) | SC3_x (m) | SC3_y (m) | SC3_z (m)

## Simulate data with detector response to full galaxy and Gaussian instrument noise
Take a list of simulated sources and a set of spacecraft coordinates to simulate TDI data streams in the frequency domain for the given observing time. The process also filters the input galaxy for sources that are possibly detectable using a (generous) cut on signal to noise ratio SNR.

`>Galaxy <source file> <orbit file> <observation time (s)>` 

Output files:

`>Galaxy_XAE.dat`: Simulated Fourier domain TDI data with instrument noise and LISA response. Columns are

    f (Hz) | Re{X(f)} | Im{X(f)} | Re{A(f)} | Im{A(f)} | Re{E(f)} | Im{E(f)} 
    
where the X channel is for a 4-link (single interferometer) data stream while A and E are the (approximately) orthogonal TDI channels constructed from the 6-link data

`Bright.dat`: List of plausibly detectable sources from input population, using an approximate analytic expression for the signal to noise ratio, compared only to instrument noise, and thresholding on `SNR>~2`. Columns are the same as the [input galaxy](#params).


## Take data file and produce fit to confusion noise
Use a running median to determine the frequency-dependent noise level of the input data.

`>Confusion_Fit <data file> <orbit file>`

Output files:

`Confusion_XAE_0.dat`: Resulting estimate of the overall and confusion noise levels. Columns are <a name="conf"></a>

    f (Hz) | SnX(f) | SnX_conf(f) | SnA(f) | SnA_conf(f)
    
 where `SnY` is the total (instrument plus confusion) noise power spectral density of channel `Y`, and `SnY_conf` is the confusion noise component. The `E` channel noises are identical to their `A` channel counterparts. 

### Intermediate data products for checking performance
`Xfit.dat`

`Xf.dat`

`Afit.dat`

## Compute residual data and detectable catalog
Steps through list of input sources, computes their matched filter SNR compared to the current estimate of the confusion noise level, and selects those with an SNR>7 as being detectable. Detectable binaries are removed from the data stream and the median noise estimate is recomputed.

`>Bright_Remove <data file> <noise file> <source file> <orbit file>`

Output files:

`Galaxy_XAE_R1.dat`: Residual frequency-domain data after detectable binaries have been removed. Columns are 

    f (Hz) | Re{X_res(f)} | Im{X_res(f)} | Re{A_res(f)} | Im{A_res(f)} | Re{E_res(f)} | Im{E_res(f)} 
    
where the X_res channel is for a 4-link (single interferometer) residual data stream while A_res and E_res are the (approximately) orthogonal TDI channels constructed from the 6-link data residuals.

`Confusion_XAE_1.dat`: Final estimate of confusion noise after `SNR>7` binaries have been removed from the data. Columns are the same as [the inital estimate](#conf).

`Confusion_XAE_DS.dat`: Same as `Confusion_XAE_1.dat` but downsampled in frequency by a factor of `100` to make the file a more manageable size.

`BrightX.dat`: List of detectable binares by a 4-link (single interferometer) detector. Columns are the same as [input galaxy](#params).

`BrightAE.dat` List of detectable binares by a 6-link (three interferometer) detector. Columns are the same as [input galaxy](#params).

### Intermediate data products for checking performance (overwriting output from `ConfusionFit`)
`Xfit.dat`

`Xf.dat`

`Afit.dat`


## Estimate errors
From list of detectable binaries, make final detection assessment using input confusion noise estimate and, for detectable binaries, estimate statistical error of parameter estimation with the Fisher approximation to the parameter covariance matrix.
The process assesses tallies statistics of the catalog including overall number of detectable sources and the number which are well localized and/or with constrained orbital dynamics.

`>Fisher_Galaxy <source file> <noise file> <4-link results file> <6-link results file> <orbit file>`

Output files:

`<4-link results file>`,`<6-link results file>`: List of detectable binaries and estimated statistical uncertainties for [4/6]-link ([1/3] interferometer) mission. Columns are

```bash
 1. f (Hz) 
 2. co-latitude (rad) 
 3. longitude (rad) 
 4. amplitude 
 5. inclination (rad) 
 6. polarization angle (rad) 
 7. ref phase (rad)`
 8. df/dt (s^-2)
 9. d^2f/dt^2 (s^-3)
 10. sigma_f (fractional)
 11. sigma_co-latitude (radians)
 12. sigma_longitude (radians)
 13. sigma_amplitude (fractional)
 14. sigma_inclination (radians)
 15. sigma_polarization (radians)
 16. sigma_phase (radians)
 17. sigma_fdot (fractional)
 18. sigma_fddot (fractional)
 19. sigma_omega (sq. deg.
 20. SNR
```

`Brightest.dat`: 100 highest-SNR signals and their errors.  Columns are the same as the above results files 

`DrawAE.dat`: **UNDER CONSTRUCTION** Fair draws from Gaussian using Fisher as covariance matrix for all detectable binaries. Colums are the same as [input galaxy](#params).

`sky_all.dat`: 3D location (in cartesian coordinates with x-y plane in the ecliptic) of all detected binaries. Columns are

    x_ecliptic (kpc) | y_ecliptic (kpc) |z_ecliptic (kpc)

`sky_3D.dat`: 3D location (in cartesian coordinates with x-y plane in the ecliptic) of all detected binaries **with distance measured to within 10% and sky location measured to within 1 sq. deg**. Columns are

    x_ecliptic (kpc) | y_ecliptic (kpc) |z_ecliptic (kpc)


## Example for input file `FisherGalaxy_LDC_Radler_galaxy_key.dat`:
```bash
#Galaxy simulation
GALAXYFILE=FisherGalaxy_LDC_Radler_galaxy_key.dat

#Observation time
T=31457280 #1 year

#Calculate LISA orbits
OrbitFile

#Simulate LISA response to galaxy
Galaxy $GALAXYFILE EccentricInclined.txt $T

#Get initial estiamte of confusion noise level
Confusion_Fit Galaxy_XAE.dat EccentricInclined.txt

#Find and regress resolvable binaries
Bright_Remove Galaxy_XAE.dat Confusion_XAE_0.dat Bright.dat EccentricInclined.txt

#Compute error estimates for resolvable binaries
Fisher_Galaxy BrightAE.dat Confusion_XAE_1.dat SigmasX.dat SigmasAE.dat DrawAE.dat EccentricInclined.txt 
```

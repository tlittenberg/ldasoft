## Input file format 
Space-separated ASCII file, one source per row. Columns are:

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
`./OrbitFile`

Output files:

`EccentricInclined.txt`: Cartesian coordinates x(t),y(t),z(t) in solar system barycentric coordinates of each spacecraft. Columns are:

    t (s) | SC1_x (m) | SC1_y (m) | SC1_z (m) | SC2_x (m) | SC2_y (m) | SC2_z (m) | SC3_x (m) | SC3_y (m) | SC3_z (m)

## Simulate data with detector response to full galaxy and Gaussian instrument noise
`./Galaxy <source file> <orbit file> <observation time (s)>`

Outputs files:

`Galaxy_XAE.dat`

`Bright.dat`


## Take data file and produce fit to confusion noise
`./Confusion_Fit <data file> <orbit file>`

Outputs files:

`Confusion_XAE_0.dat`

`Xfit.dat`

`Xf.dat`

`Afit.dat`

## Compute residual data and detectable catalog
`./Bright_Remove <data file> <noise file> <source file> <orbit file>`

Output files:

`Confusion_XAE_1.dat`

`Confusion_XAE_DS.dat`

`BrightX.dat`

`BrightAE.dat`


## Estimate errors
`./Fisher_Galaxy <source file> <noise file> <4-link results file> <6-link results file> <orbit file>`

Output files:

`SigmasAE.dat`

`SigmasX.dat`

`Brightest.dat`

`DrawAE.dat`

`sky_all.dat`

`sky_3D.dat`

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

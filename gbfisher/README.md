## Input file format 
Space-separated ASCII file, one source per row. Columns are:
`f (Hz) |  df/dt (s^-2) | co-latitude (rad) | longitude (rad) | amplitude | inclination (rad) | polarization angle (rad) | ref phase (rad)`

## Create detector.h for different mission configurations
`./Setup kappa, lambda, Larm (m), Sacc (m^2 s^-4 Hz^-1), Spos (m^2 Hz^-1), Cadence (s), Nsamples`

## Create orbit file for mission configuration EccentricInclined.txt
`./OrbitFile`

## Simulation of the full galaxy Galaxy_XAE.dat with noise
`./Galaxy <galaxy file> <orbit file>`

## Take galaxy simulation and produce fit to confusion noise
`./Confusion_Fit <input file> <orbit file>`

## Remove bright sources
`./Bright_Remove XAE.dat Noise.dat Bright.dat <orbit file>`

## Estimate errors
`./Fisher_Galaxy detections.dat confusion.dat sigmasX.dat sigmasAE.dat <orbit file>`

## Example for input file Catalog_NoID.txt:
```bash
gcc -O3 -o Setup Setup.c -lm
./Setup 0 0 2.5e9 9.e-30 8.321e-23 1 31457280
gcc -O3 -o OrbitFile OrbitFile.c arrays.c -lm
gcc -O3 -o Galaxy Galaxy.c arrays.c -lm 
gcc -O3 -o Confusion_Fit Confusion_Fit.c arrays.c -lm
gcc -O3 -o Bright_Remove Bright_Remove.c arrays.c -lm
```

```bash
./OrbitFile
./Galaxy Catalog_NoID.txt EccentricInclined.txt
./Confusion_Fit Galaxy_XAE.dat EccentricInclined.txt
./Bright_Remove Galaxy_XAE.dat Confusion_XAE_0.dat Bright.dat EccentricInclined.txt
```

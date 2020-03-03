# Example use cases for GBMCMC

## Analyze on-the-fly data & source simulation
`gb_mcmc --inj <path to>/ldasoft/gbmcmc/etc/sources/precision/PrecisionSource_0.txt`

## Analyze LDC simulations
See `src/README.md`

## Analyze verification binary with EM "priors"
Run with sky-location fixed
`gb_mcmc --inj <path to>/ldasoft/gbmcmc/etc/sources/verification/ZTF1539.dat --known-source --no-rj --duration 62914560 --cheat --samples 128`

Run with sky-location and orbital period fixed
`gb_mcmc --inj /Users/tlittenb/ldasoft/galactic_binaries/etc/sources/verification/ZTF1539_tides.dat --known-source --no-rj --duration 62914560 --cheat --samples 128 --fix-freq`

## Use UCBs as calibration sources
`gb_mcmc --inj <path to>/ldasoft/gbmcmc/etc/sources/calibration/CalibrationSource_0.txt --duration 604948 --segments 2 --fit-gap --samples 128 --no-rj`

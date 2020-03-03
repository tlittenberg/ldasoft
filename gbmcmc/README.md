# Example use cases for GBMCMC
Get help with `gb_mcmc --help`

## Analyze on-the-fly data & source simulation
```bash
gb_mcmc \
  --inj <path to>/ldasoft/gbmcmc/etc/sources/precision/PrecisionSource_0.txt
```

## Analyze LDC simulations
See `src/README.md`

## Analyze verification binary with EM "priors"
Run with sky-location fixed
```bash
gb_mcmc \
  --inj <path to>/ldasoft/gbmcmc/etc/sources/verification/ZTF1539.dat \
  --known-source \
  --no-rj \
  --cheat \
  --samples 128
```

Run with sky-location and orbital period fixed
```bash
gb_mcmc \
  --inj /Users/tlittenb/ldasoft/galactic_binaries/etc/sources/verification/ZTF1539_tides.dat \
  --known-source \
  --no-rj \
  --cheat \
  --samples 128 \
  --fix-freq
```

## Use UCBs as calibration sources
```bash
gb_mcmc \
  --inj <path to>/ldasoft/gbmcmc/etc/sources/calibration/CalibrationSource_0.txt \
  --duration 604948 \
  --segments 2 \
  --fit-gap \
  --samples 128 \
  --no-rj
```

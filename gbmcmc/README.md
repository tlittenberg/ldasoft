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

## Output format

The `dimension_chain.dat.N` files contain the post-burnin posterior samples of
the binary parameters for the subset of models that have exactly `N` binaries.
The columns are

    f fdot A cos_colat long cos_inc psi phi

If there are more than one source, the sources cycle in the output (so line 1 is
the first source, first posterior sample; line 2 is the second source, first
posterior sample, etc).

The `parameter_chain.dat.M` files should be used in concert with the
`model_chain.dat.M` files.  `model_chain.dat.M` records how many binaries were
in the model at each MCMC step; there are a corresponding number of lines in
`parameter_chain.dat.M` for that step's binary parameters (i.e.
`parameter_chain.dat.M` contains a mashed-together compilation of the
`dimension_chain.dat.N` files).  The index `M` for these files refers to the
*temperature* of the parallel-tempered chain that is producing these outputs (`M
= 0` is the target distribution; `M > 0` is softened versions of this
distribution; if you don't know what this means, you want `M=0`).  Unline
`dimension_chain.dat.N` files, these files contain *all* samples including
burnin.  Unless you care deeply about checkpointing or run diagnostics, you do
not need these files. Use `dimension_chain.dat.N` files, as described above.

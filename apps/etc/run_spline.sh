#!/bin/bash

noise_sampler="${GLASS_PREFIX}/noise_spline_mcmc"
outdir="${GLASS_TEST}/noise_spline_mcmc"
fmin=1e-4 #Hz
fmax=1e-2
Tobs=2621440 #one month
threads=12

cmd="${noise_sampler} --sim-noise --conf-noise --fmin ${fmin} --fmax ${fmax} --duration ${Tobs} --thread ${threads} --rundir ${outdir}"

echo $cmd

$cmd

#!/bin/bash

noise_sampler="glass_noise_mcmc"
fmin=1e-4 #Hz
fmax=1e-2
Tobs=2621440 #one month
threads=12

cmd="${noise_sampler} --sim-noise --conf-noise --fmin ${fmin} --fmax ${fmax} --duration ${Tobs} --thread ${threads} --rundir noise_mcmc"

echo $cmd

$cmd

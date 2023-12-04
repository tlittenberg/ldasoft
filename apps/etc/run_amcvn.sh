#!/bin/bash

ucb_sampler="${GLASS_PREFIX}/ucb_mcmc"
injection_file="${GLASS_REPO}/ucb/etc/sources/verification/amcvn.dat"
outdir="${GLASS_TEST}/amcvn"

threads=12

cmd="${ucb_sampler} --known-source --rundir ${outdir} --no-rj --cheat --threads ${threads} --inj ${injection_file}"

echo $cmd

$cmd

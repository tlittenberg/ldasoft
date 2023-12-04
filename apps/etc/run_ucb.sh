#!/bin/bash

ucb_sampler="${GLASS_PREFIX}/ucb_mcmc"
injection_file="${GLASS_REPO}/ucb/etc/sources/precision/PrecisionSource_1.txt"
outdir="${GLASS_TEST}/ucb"
threads=12

cmd="${ucb_sampler} --rundir ${outdir} --cheat --threads ${threads} --inj ${injection_file}"

echo $cmd

$cmd

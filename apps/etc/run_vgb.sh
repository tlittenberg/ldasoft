#!/bin/bash

vgb_sampler="${GLASS_PREFIX}/vgb_mcmc"
vgb_list="${GLASS_REPO}/apps/etc/verification_binaries.dat"
outdir="${GLASS_TEST}/vgb_mcmc"
bandwidth=128
threads=12

cmd="${vgb_sampler} --known-sources ${vgb_list} --samples ${bandwidth} --threads ${threads} --rundir ${outdir}"

echo $cmd

$cmd

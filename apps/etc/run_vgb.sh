#!/bin/bash

vgb_sampler="glass_vgb_mcmc"
vgb_list="~/ldasoft/apps/etc/verification_binaries.dat"
bandwidth=128
threads=12

cmd="${vgb_sampler} --known-sources ${vgb_list} --samples ${bandwidth} --threads ${threads} --rundir vgb_mcmc"

echo $cmd

$cmd

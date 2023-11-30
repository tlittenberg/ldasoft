#!/bin/bash

ucb_sampler="glass_ucb_mcmc"
injection_file="~/ldasoft/ucb/etc/sources/verification/amcvn.dat"
threads=12

cmd="${ucb_sampler} --known-source --rundir amcvn --no-rj --cheat --threads ${threads} --inj ${injection_file}"

echo $cmd

$cmd

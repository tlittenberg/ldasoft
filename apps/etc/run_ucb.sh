#!/bin/bash

ucb_sampler="glass_ucb_mcmc"
injection_file="~/ldasoft/ucb/etc/sources/precision/PrecisionSource_1.txt"
threads=12

cmd="${ucb_sampler} --rundir ucb --cheat --threads ${threads} --inj ${injection_file}"

echo $cmd

$cmd

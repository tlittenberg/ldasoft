#!/bin/bash

nohup bash run_ucb.sh > ucb.out 
nohup bash run_amcvn.sh > amcvn.out 
nohup bash run_noise.sh > noise.out
nohup bash run_spline.sh > spline.out 
nohup bash run_vgb.sh > vbmcmc.out
#nohup bash run_globalfit.sh > globalfit.out 

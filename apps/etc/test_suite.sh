#!/bin/bash

export GLASS_PREFIX="/Users/tyson/opt/bin"
export GLASS_REPO="/Users/tyson/ldasoft"
export GLASS_TEST="/Users/tyson/ldasoft/test"

mkdir -p ${GLASS_TEST}

bash run_ucb.sh | tee ${GLASS_TEST}/ucb.out 
#bash run_amcvn.sh | tee ${GLASS_TEST}/amcvn.out 
#bash run_noise.sh | tee ${GLASS_TEST}/noise.out
#bash run_spline.sh | tee ${GLASS_TEST}/spline.out 
#bash run_vgb.sh | tee ${GLASS_TEST}/vgb.out
#bash run_globalfit.sh | tee ${GLASS_TEST}/globalfit.out 

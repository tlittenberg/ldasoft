#!/bin/bash

export GLASS_PREFIX="/Users/tlittenb/ldasoft/reorg/bin"
export GLASS_REPO="/Users/tlittenb/ldasoft"
export GLASS_TEST="/Users/tlittenb/ldasoft/test"


#bash run_ucb.sh | tee ucb.out 
#bash run_amcvn.sh | tee amcvn.out 
#bash run_noise.sh | tee noise.out
#bash run_spline.sh | tee spline.out 
#bash run_vgb.sh | tee vgb.out
bash run_globalfit.sh | tee globalfit.out 

#!/usr/bin/env bash
LANG=en_US
export LANG

read n_pop N K dir<<<  "$1 $2 $3 $4" 

J0=1.0
J1=0.0
sig1=1.0

sed -ie "s/J0 .*/J0 -${J0}/" globals.h 
sed -ie "s/MEAN_XI .*/MEAN_XI ${J1}/" globals.h 
sed -ie "s/VAR_XI .*/VAR_XI ${sig1}/" globals.h 

sed -ie "s/IF_TRIALS .*/IF_TRIALS 1/" globals.h 

sed -ie "s/FIXED_XI .*/FIXED_XI 1/" globals.h 

for trial in $(seq 0 1 9); do 

    sed -ie "s/TRIAL_ID .*/TRIAL_ID ${trial}/" globals.h 

    echo "#########################################################################"
    echo "trial ${trial}"
    echo "#########################################################################"

    echo "generating connectivity matrix ..."
    (cd ../../../cuda/connectivity/; make -B &>/dev/null)
    (cd ../../../cuda/connectivity/; ./a.out &>/dev/null)
    
    echo "compiling ..."
    g++ main.cpp -std=c++11 -Ofast -s -o loop_trials.out
    
    screen -dmS N${N}_K${K}_J0${i}_MEAN_XI_${J1}_VAR_XI_${sig1}_TRIAL_${trial} ./loop_trials.out $n_pop $N $K $dir
    echo "running J0_${i}_MEAN_XI_${J1}_VAR_XI_${sig1}_trial_${trial} "
done

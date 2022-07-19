#!/usr/bin/env bash 

read n_pop N K dir <<< "$1 $2 $3 $4" 

make 
screen -dmS a_${n_pop}_pop_${dir}_N_${N}_K_${K}_off ./a.out $n_pop $N $K ${dir}_off
screen -dmS b_${n_pop}_pop_${dir}_N_${N}_K_${K}_on ./a.out $n_pop $N $K ${dir}_on 

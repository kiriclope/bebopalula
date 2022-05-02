#!/usr/bin/env bash 

temp_files=$(./bash_utils.sh)
temp_globals=$(echo $temp_files | awk 'NR==1{print $1}') 
temp_main=$(echo $temp_files | awk 'NR==1{print $2}') 
temp_out=$(echo $temp_files | awk 'NR==1{print $3}') 

IF_LIF=1 
IF_BIN=0 

IF_STP=1 
IF_GEN_CON=0 

RANK=1
IF_GAUSS=1 
IF_SPEC=0 
FIX_MAP_SEED=1 

IF_LOW_RANK=0 
FIX_KSI_SEED=1  

sed -ie "s/ DURATION .*/ DURATION (double) 60E3 /" "$temp_globals" ; 
sed -ie "s/ TIME_INI .*/ TIME_INI (double) 0E3 /" "$temp_globals" ; 
sed -ie "s/ TIME_STEADY .*/ TIME_STEADY (double) 30E3 /" "$temp_globals" ; 
sed -ie "s/ TIME_WINDOW .*/ TIME_WINDOW (double) .250E3 /" "$temp_globals" ; 
sed -ie "s/ TIME_REC .*/ TIME_REC (double) 60E3 /" "$temp_globals" ; 
sed -ie "s/ TIME_REC_SPIKES .*/ TIME_REC_SPIKES (double) 0E3 /" "$temp_globals" ; 

sed -ie "s/ IF_LIF .*/ IF_LIF ${IF_LIF} /" "$temp_globals" ; 
sed -ie "s/ IF_BIN .*/ IF_BIN ${IF_BIN} /" "$temp_globals" ; 

sed -ie "s/ IF_STP .*/ IF_STP ${IF_STP} /" "$temp_globals" ; 

sed -ie "s/ IF_TRIALS .*/ IF_TRIALS 1 /" "$temp_globals" ; 
sed -ie "s/ IF_INI_COND .*/ IF_INI_COND 0 /" "$temp_globals" ; 

sed -ie "s/ IF_GEN_CON .*/ IF_GEN_CON ${IF_GEN_CON} /" "$temp_globals" ; 
sed -ie "s/ IF_SAVE_CON_VEC .*/ IF_SAVE_CON_VEC 0 /" "$temp_globals" ; 
sed -ie "s/ IF_SAVE_SPARSE_REP .*/ IF_SAVE_SPARSE_REP 0 /" "$temp_globals" ; 

sed -ie "s/ RANK .*/ RANK ${RANK} /" "$temp_globals" ; 

sed -ie "s/ IF_SPEC .*/ IF_SPEC ${IF_SPEC} /" "$temp_globals" ; 
sed -ie "s/ FIX_MAP_SEED .*/ FIX_MAP_SEED ${FIX_MAP_SEED} /" "$temp_globals" ; 

sed -ie "s/ IF_LOW_RANK .*/ IF_LOW_RANK ${IF_LOW_RANK} /" "$temp_globals" ; 
sed -ie "s/ FIX_KSI_SEED .*/ FIX_KSI_SEED ${FIX_KSI_SEED} /" "$temp_globals" ; 

sed -ie "s/ IF_HYSTERESIS .*/ IF_HYSTERESIS 0 /" "$temp_globals" ; 
sed -ie "s/ IF_DPA .*/ IF_DPA 0 /" "$temp_globals" ; 
sed -ie "s/ IF_DUAL .*/ IF_DUAL 0 /" "$temp_globals" ; 
sed -ie "s/ IF_CHRISTOS .*/ IF_CHRISTOS 0 /" "$temp_globals" ; 

sed -ie "s/ IF_STEP .*/ IF_STEP 1 /" "$temp_globals" ; 

read n_pop N K dir n_trials <<<  "$1 $2 $3 $4 $5" 

for trial in $(seq 1 1 $n_trials); do 
    
    echo "#########################################################################" 
    ./mem_usage.sh 
    ./cpu_usage.sh 
    echo "#########################################################################" 

    dum=$(echo "print(${trial} / ${n_trials})" | python3)
    sed -ie "s/ PHI_EXT (double) .*/ PHI_EXT (double) $dum /" "$temp_globals" ; 
    
    echo "simulation parameters:" 
    echo "n_pop ${n_pop} n_neurons ${N}0000 K ${K} ${dir} trial ${trial}" 
    echo "#########################################################################" 
    
    sed -ie "s/ TRIAL_ID .*/ TRIAL_ID ${trial} /" "$temp_globals" 
    
    echo "g++ $temp_main -Ofast -s -std=c++11 -o $temp_out -lgsl -lblas" 
    g++ ${temp_main} -Ofast -s -std=c++11 -o ${temp_out} -lgsl -lblas 
    # ./${temp_out} $n_pop $N $K $dir 
    screen -dmS ${n_pop}_pop_${dir}_N_${N}_K_${K}_trial_${trial} ./${temp_out} $n_pop $N $K $dir 
    
    # sleep 5s 
done

rm $temp_files

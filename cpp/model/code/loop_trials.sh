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
IF_GAUSS=0
IF_RING=1

IF_SPEC=0
FIX_MAP_SEED=1 

IF_LOW_RANK=0 
FIX_KSI_SEED=1 

sed -ie "s/ DURATION .*/ DURATION (double) 10E3 /" "$temp_globals" ; 
sed -ie "s/ TIME_INI .*/ TIME_INI (double) 0E3 /" "$temp_globals" ; 
sed -ie "s/ TIME_STEADY .*/ TIME_STEADY (double) 10E3 /" "$temp_globals" ; 
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

sed -ie "s/ IF_CHRISTOS .*/ IF_CHRISTOS 1 /" "$temp_globals" ; 
sed -ie "s/ IF_STEP .*/ IF_STEP 0 /" "$temp_globals" ; 

sed -ie "s/ IF_CON_DIR .*/ IF_CON_DIR 1 /" "$temp_globals" ; 

read n_pop N K dir n_trials cue <<<  "$1 $2 $3 $4 $5 $6" 

sed -ie "s/ CUE .*/ CUE (double) ${cue} /" "$temp_globals" ;

# generating connectivity with cuda 
( cd ../../../cuda/connectivity2/ ; 
  sed -ie "s/ n_pop .*/ n_pop ${n_pop} /" globals.h ;
  sed -ie "s/ N_NEURONS .*/ N_NEURONS (unsigned long) ${N}0000 /" globals.h ;
  sed -ie "s/ K .*/ K (double) ${K} /" globals.h ; 
  sed -ie "s/ IF_RING .*/ IF_RING ${IF_RING} /" globals.h ;
  sed -ie "s/ IF_SAVE_SPARSE_REP .*/ IF_SAVE_SPARSE_REP 1 /" globals.h ; 
  sed -ie "s/ IF_CON_DIR .*/ IF_CON_DIR 1 /" globals.h ; 
  sed -ie "s/ CUE .*/ CUE (double) ${cue} /" globals.h ; 
  make > /dev/null 2>&1 ; 
) 

for trial in $(seq 1 1 $n_trials); do
    
    echo "#########################################################################" 
    ./mem_usage.sh 
    # ./cpu_usage.sh 
    echo "#########################################################################" 
    
    # ## generating connectivity with cpp
    # sed -ie "s/ SEED_CON .*/ SEED_CON (double) ${trial} /" "$temp_globals" ; 
    # sed -ie "s/ IF_LIF .*/ IF_LIF 0 /" "$temp_globals" ; 
    # sed -ie "s/ IF_GEN_CON .*/ IF_GEN_CON 1 /" "$temp_globals" ; 
    # sed -ie "s/ IF_SAVE_SPARSE_REP .*/ IF_SAVE_SPARSE_REP 1 /" "$temp_globals" ; 
    
    # g++ -L/home/leon/bebopalula/cpp/libs/gsl/lib -I/home/leon/bebopalula/cpp/libs/gsl/include -std=c++11 ${temp_main} -Ofast -s -o matrix.out -lgsl -lgslcblas 
    # ./matrix.out $n_pop $N $K ${dir}_off 
    # # srun --priority="TOP" ./matrix.out $n_pop $N $K ${dir}_off
    
    # sed -ie "s/ IF_LIF .*/ IF_LIF 1 /" "$temp_globals" ; 
    # sed -ie "s/ IF_GEN_CON .*/ IF_GEN_CON 0 /" "$temp_globals" ; 
    # sed -ie "s/ IF_SAVE_SPARSE_REP .*/ IF_SAVE_SPARSE_REP 0 /" "$temp_globals" ; 
    
    # generating connectivity with cuda 
    echo "generating connectivity trial ${trial} cue ${cue}:" 
    echo "#########################################################################" 
    ( cd ../../../cuda/connectivity2/ ; 
      # ./a.out ; 
      ./a.out > /dev/null ; 
    ) 
    
    echo "#########################################################################" 
    echo "simulation parameters:" 
    echo "n_pop ${n_pop} n_neurons ${N}0000 K ${K} ${dir} trial ${trial}" 
    echo "#########################################################################" 
    
    sed -ie "s/ TRIAL_ID .*/ TRIAL_ID ${trial} /" "$temp_globals" 
    sed -ie "s/ PHI_CUE (double) .*/ PHI_CUE (double) .375 /" "$temp_globals" ; 
    sed -ie "s/ PHI_ERASE (double) .*/ PHI_ERASE (double) 1.0-.375 /" "$temp_globals" ; 
    
    echo "#########################################################################" 
    echo "compiling close condition"
    echo "#########################################################################" 
    
    g++ -L/home/leon/bebopalula/cpp/libs/gsl/lib -I/home/leon/bebopalula/cpp/libs/gsl/include -std=c++11 ${temp_main} -Ofast -s -o ${temp_out}_${trial}.out -lgsl -lgslcblas 

    # ./${temp_out}_${trial}.out $n_pop $N $K ${dir}_off 
    screen -dmS ${n_pop}_pop_${dir}_N_${N}_K_${K}_trial_${trial}_off_close ./${temp_out}_${trial}.out $n_pop $N $K ${dir}_off 
    screen -dmS ${n_pop}_pop_${dir}_N_${N}_K_${K}_trial_${trial}_on_close ./${temp_out}_${trial}.out $n_pop $N $K ${dir}_on 
    # screen -dmS ${n_pop}_pop_${dir}_N_${N}_K_${K}_trial_${trial}_off_close srun --priority="TOP" ./${temp_out}_${trial}.out $n_pop $N $K ${dir}_off 
    # screen -dmS ${n_pop}_pop_${dir}_N_${N}_K_${K}_trial_${trial}_on_close srun --priority="TOP" ./${temp_out}_${trial}.out $n_pop $N $K ${dir}_on 
    
    echo "#########################################################################" 
    ./mem_usage.sh     
    echo "#########################################################################" 
    echo "compiling far condition"
    echo "#########################################################################"
    
    sed -ie "s/ PHI_CUE (double) .*/ PHI_CUE (double) .25 /" "$temp_globals" ; 
    sed -ie "s/ PHI_ERASE (double) .*/ PHI_ERASE (double) .75 /" "$temp_globals" ;
    
    g++ -L/home/leon/bebopalula/cpp/libs/gsl/lib -I/home/leon/bebopalula/cpp/libs/gsl/include -std=c++11 ${temp_main} -Ofast -s -o ${temp_out}_${trial}.out -lgsl -lgslcblas 
    
    screen -dmS ${n_pop}_pop_${dir}_N_${N}_K_${K}_trial_${trial}_off_far ./${temp_out}_${trial}.out $n_pop $N $K ${dir}_off    
    screen -dmS ${n_pop}_pop_${dir}_N_${N}_K_${K}_trial_${trial}_on_far ./${temp_out}_${trial}.out $n_pop $N $K ${dir}_on 
    # screen -dmS ${n_pop}_pop_${dir}_N_${N}_K_${K}_trial_${trial}_off_far srun --priority="TOP" ./${temp_out}_${trial}.out $n_pop $N $K ${dir}_off 
    # screen -dmS ${n_pop}_pop_${dir}_N_${N}_K_${K}_trial_${trial}_on_far srun --priority="TOP" ./${temp_out}_${trial}.out $n_pop $N $K ${dir}_on 
    
done

rm $temp_files

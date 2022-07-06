#!/usr/bin/env bash 

read n_pop N K kappa <<<  "$1 $2 $3 $4" 

# generating connectivity with cuda 
( cd ../../../cuda/connectivity/ ; 
  sed -ie "s/ n_pop .*/ n_pop ${n_pop} /" globals.h ;
  sed -ie "s/ N_NEURONS .*/ N_NEURONS (unsigned long) ${N}0000 /" globals.h ;
  sed -ie "s/ K .*/ K (double) ${K} /" globals.h ; 
  sed -ie "s/ IF_RING .*/ IF_RING 1 /" globals.h ;
  sed -ie "s/ KAPPA (double) .*/ KAPPA (double) ${kappa} /" globals.h ;
  sed -ie "s/ IF_SAVE_SPARSE_REP .*/ IF_SAVE_SPARSE_REP 1 /" globals.h ; 
  sed -ie "s/ IF_CON_DIR .*/ IF_CON_DIR 0 /" globals.h ; 
  make > /dev/null 2>&1 ; 
  ./a.out ; 
) 

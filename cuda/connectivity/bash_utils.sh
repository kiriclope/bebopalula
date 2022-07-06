#!/bin/bash

read rand <<< "$1" 

temp_cuda_main=$(printf 'temp_cuda_main_%d.cu' $rand) 
temp_cuda_globals=$(printf 'temp_cuda_globals_%d.h' $rand) 
temp_cuda_out=$(printf 'temp_cuda_%d_a' $rand) 

cp main.cu $temp_cuda_main 
cp globals.h $temp_cuda_globals 

sed -i '5s/ .*/ "'"$temp_cuda_globals"'" /' "$temp_cuda_main" ; 

head $temp_cuda_main

echo $temp_cuda_globals $temp_cuda_main $temp_cuda_out

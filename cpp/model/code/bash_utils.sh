#!/bin/bash

read rand <<< "$1" 

temp_main=$(printf 'temp_%d_main.cpp' $rand) 
temp_globals=$(printf 'temp_%d_globals.h' $rand) 
temp_out=$(printf 'temp_%d_a' $rand) 

cp main.cpp $temp_main 
cp globals.h $temp_globals 

sed -i 's/ "globals.h" .*/ "'"$temp_globals"'" /' "$temp_main" ; 

# echo $temp_globals $temp_main $temp_out


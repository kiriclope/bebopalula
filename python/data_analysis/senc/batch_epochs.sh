#!/usr/bin/env bash

declare -a days=('last' 'first')
declare -a objects=('score' 'cos')

# for mouse in 2 3 1; do
# for day in ${days[@]}; do
# echo $mouse $day
# python3 sel_epochs.py $mouse all $day 'correct' cos sample 'LD' > cos.txt
#     done
# done

python3 sel_epochs.py 2 all 'all' 'correct' cos sample 'ED'

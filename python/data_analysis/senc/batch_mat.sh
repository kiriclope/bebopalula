#!/usr/bin/env bash

declare -a days=('all' 'first' 'last')
declare -a objects=('score' 'cos')
declare -a tasks=('DPA' 'DualGo' 'DualNoGo')

for mouse in 1 2 3; do
    for day in ${days[@]}; do
	for task in ${tasks[@]}; do
	    python3 cross_temp_mat.py $mouse $task $day 'correct'
	done
    done
done

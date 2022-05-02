#!/usr/bin/env bash

declare -a days=('last' 'first')
declare -a tasks=('DPA' 'DualGo' 'DualNoGo')
declare -a trials=('correct' 'incorrect')

for mouse in 3 ; do
    for day in ${days[@]}; do
	for task in ${tasks[@]}; do
	    echo $mouse $day $task
	    python3 overlap_2D.py $mouse $task $day 'correct'
	done
    done
done

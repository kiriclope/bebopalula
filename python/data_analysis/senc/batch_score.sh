#!/usr/bin/env bash

declare -a days=('last' 'first')
declare -a objects=('score' 'cos')

for mouse in 2 3 1; do
    for day in ${days[@]}; do
	echo $mouse $day
	python3 sel_epochs.py $mouse all $day 'correct' score sample 'ED' > batch_score_ED.txt 
	python3 sel_epochs.py $mouse all $day 'correct' score sample 'LD' > batch_score_LD.txt 
	python3 sel_epochs.py $mouse all $day 'correct' cos sample 'ED' > batch_cos_ED.txt 
    done
done

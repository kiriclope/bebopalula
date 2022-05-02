#!/usr/bin/env bash

declare -a days=('last' 'first')

# for mouse in 1 2 3; do
#     for day in ${days[@]}; do
# 	echo $mouse $day
# 	# python3 sel_time.py $mouse all $day 'correct' cos sample 'ED' > batch_time.txt 
# 	python3 sel_time.py $mouse all $day 'correct' score sample 'ED' > batch_time.txt 
# 	python3 sel_time.py $mouse all $day 'correct' score sample 'LD' > batch_time.txt 
# 	# python3 sel_time.py $mouse all $day 'correct' non_zero sample 'ED' 
# 	# python3 sel_time.py $mouse all $day 'correct' proj sample 'ED' 
#     done 
# done 

# for mouse in 2 3; do
#     python3 sel_time.py $mouse all all 'correct' score sample 'ED' > batch_time.txt
# done

# for mouse in 1 2 3; do
#     for day in ${days[@]}; do
# 	python3 sel_time.py $mouse all $day 'correct' score sample 'LD' > batch_time.txt 
#     done
# done

for day in ${days[@]}; do
    python3 sel_time.py 2 all $day 'correct' cos sample 'ED'  
    # python3 sel_time.py 2 all $day 'correct' score sample 'LD' > batch_time_2.txt 
done

#!/usr/bin/env bash

python sel_epochs.py 1 all 'first' 'correct' cos sample 
python sel_epochs.py 1 all 'last' 'correct' cos sample 
python sel_epochs.py 1 all 'all' 'correct' cos sample 

python sel_time.py 1 all 'first' 'correct' cos sample 
python sel_time.py 1 all 'last' 'correct' cos sample 
python sel_time.py 1 all 'all' 'correct' cos sample 



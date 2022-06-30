#!/bin/bash

mem_usage=`free -m | awk 'NR==2{print $3/$2*100 }'` 
mem_usage=$( printf "%.0f" $mem_usage )
echo -e "mem_usage" $mem_usage "%" 

if [ $mem_usage -gt 60 ]; then
    echo " MEM_USAGE > 60.0, sleeping for a while ..."    
fi

while [ $mem_usage -gt 60 ]; 
do
    sleep 2s ; 
    mem_usage=`free -m | awk 'NR==2{print $3/$2*100 }'` ;
    mem_usage=$( printf "%.0f" $mem_usage )
    
    if [ $mem_usage -gt 95 ]; then
	echo " MEM_USAGE > 95.0, killing everything ..."
	pkill screen
    fi
    
done

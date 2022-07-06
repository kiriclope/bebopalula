#!/bin/bash

mem_usage=`free -m | awk 'NR==2{print $3/$2*100 }'` 
mem_usage=$( printf "%.0f" $mem_usage )
echo "mem_usage" $mem_usage "%"

if [ $mem_usage -gt 65 ]; then
    echo " MEM_USAGE > 65.0, sleeping for a while ..."    
fi

while [ $mem_usage -gt 65 ]; 
do
    sleep 2s ; 
    mem_usage=`free -m | awk 'NR==2{print $3/$2*100 }'` ;
    mem_usage=$( printf "%.0f" $mem_usage )
    
    if [ $mem_usage -gt 99 ]; then
	echo " MEM_USAGE > 99.0, killing everything ..."
	pkill screen
    fi
    
done

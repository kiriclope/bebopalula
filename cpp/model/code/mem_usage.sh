#!/bin/bash

mem_usage=`free -m | awk 'NR==2{print $3/$2*100 }'` 
echo -e "mem_usage" $mem_usage "%\c" 

if [ `echo "${mem_usage} > 65.0" | bc` -eq 1 ] ; then
    echo " MEM_USAGE > 65.0, sleeping for a while ..." 
fi

while [ `echo "${mem_usage} > 65.0" | bc` -eq 1 ] ; 
do
    sleep 20s ; 
    mem_usage=`free -m | awk 'NR==2{print $3/$2*100 }'` ;
done

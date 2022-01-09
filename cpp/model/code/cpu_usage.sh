#!/usr/bin/env bash 

calc(){ awk "BEGIN { print "$*" *100 }"; } 

cpu_use () {
    n_cpu=$(nproc --all) ;
    n_threads=$(top -bn1 | awk 'NR==2{print $4}') ;
    cpu_usage=$( calc $n_threads/$n_cpu ) ;
    # echo "n_cpu" $n_cpu ", n_threads" $n_threads ", usage" $cpu_usage "%" ;
    echo $cpu_usage ;
}

cpu_usage=$(cpu_use)
echo " cpu_usage" $cpu_usage "%"

if [ `echo "${cpu_usage} > 90.0" | bc` -eq 1 ] ; then
    echo " CPU_USAGE > 75.0, sleeping for a while ..." 
fi

while [ `echo "${cpu_usage} > 90.0" | bc` -eq 1 ] ;
do
    sleep 5s ;
    cpu_usage=$(cpu_use) 
done

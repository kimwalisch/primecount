#!/bin/bash

if [ $# -ne 3 ]
then
    echo "Usage example:"
    echo "$ ./kill_resume.sh 1e21 16 600"
    echo "Compute pi(1e21) using up to 16 threads and kill & resume"
    echo "the computation after at most 600 seconds."
    exit 1
fi

command -v ./primecount >/dev/null 2>/dev/null
if [ $? -ne 0 ]
then
    echo "Error: no primecount binary in current directory."
    exit 1
fi

max_threads=$2
max_seconds=$3

./primecount $1 -s &

while true
do 
    threads=$(( ( RANDOM % $max_threads )  + 1 ))
    seconds=$(( ( RANDOM % $max_seconds )  + 1 ))
    echo "sleep $seconds"
    sleep $seconds

    primecount_pid=$(ps aux | grep '\./primecount' | grep -v grep | tr -s ' ' | cut -f2 -d' ')

    if [ -z "$primecount_pid" ]
    then
        break
    fi

    kill $primecount_pid
    sleep 1
    ./primecount --resume --threads=$threads &
done

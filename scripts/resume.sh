#!/bin/bash

# Usage: ./resume.sh [x] [seconds]
# Kills primecount every n seconds and then resumes the
# computation

x=1e18
seconds=3

# Execute in base directory
if [ "$(basename $(pwd))" = "scripts" ]
then
    cd ..
fi

command -v ./primecount >/dev/null 2>/dev/null
if [ $? -ne 0 ]
then
    echo "Error: no primecount binary in current directory."
    exit 1
fi

if [ $# -ge 1 ]
then
    x=$1
    echo "x = $x"
fi

if [ $# -ge 2 ]
then
    seconds=$2
    echo "seconds = $seconds"
fi

rm -f primecount.backup
./primecount $x --S2_hard &
run=0

while [ $run -eq 0 ]
do
    sleep $seconds;
    ps_primecount=$(ps -A -o pid,args | grep primecount | grep -v grep)
    run=$?

    if [ $run -eq 0 ]
    then
        # Get primecount PID
        pid=$(echo "$ps_primecount" | awk '{print $1}')

        # Send SIGKILL (9) to simulate crash
        kill -9 $pid

        # Print current status
        echo ""
        echo "Percent: $(cat primecount.backup | grep percent | cut -d':' -f2 | cut -c2-6)"

        # Resume from primecount.backup
        ./primecount --resume &
    fi
done

echo "finished!"

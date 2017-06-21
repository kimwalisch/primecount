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
./primecount $x --S2_hard -s &
sleep $seconds;

while [ "$(grep -o -c threads primecount.backup)" -ne 0 ]
do
    sleep $seconds;
    # Send SIGKILL (9) to simulate crash
    kill -9 $(ps | grep primecount | grep -v grep | cut -f1 -d' ')
    ./primecount --resume &
done

echo "finished!"

#!/bin/bash

# Usage: ./worktodo.sh
# Iterates over the numbers in worktodo.txt (one number per line with
# optional flags) and processes them using primecount.

command -v ./primecount >/dev/null 2>/dev/null
if [ $? -ne 0 ]
then
    echo "Error: no primecount binary in current directory."
    exit 1
fi

# check if primecount supports --log option
log_flag=$(./primecount --help | grep Log > /dev/null && echo --log)

while read first_line < worktodo.txt
do
    # Skip empty lines
    if [ "$first_line" != "" ]
    then
        ./primecount $first_line $log_flag

        # Check if primecount exited successfully
        if [ $? -ne 0 ]
        then
            echo ""
            echo "Error in primecount_from_file.sh:"
            echo "The following command failed: ./primecount $first_line"
            exit 1
        fi
    fi

    # delete first line from worktodo.txt
    tail -n +2 worktodo.txt > .tmp_worktodo.txt
    cp -f .tmp_worktodo.txt worktodo.txt
done

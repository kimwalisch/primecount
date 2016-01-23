#!/bin/bash

while read first_line < worktodo.txt
do
    # Skip empty lines
    if [ "$first_line" != "" ]; then

        # Appends result into results.txt
        ./primecount $first_line --log

        # Check if primecount exited successfully
        if [ $? -ne 0 ]; then
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

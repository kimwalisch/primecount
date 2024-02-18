#!/bin/bash

# This scripts compares the results of the Li(x) implementation
# when run with and without float128 enabled. The float128 result
# is always correct whereas the long double result is frequently
# off by 1. The goal of this script is to find a threshold below
# which all long double results are correct.

if [ "$(basename $(pwd))" != "primecount" ];
then
    echo "ERROR: Run this script from the primecount root directory!"
    exit 1
fi

rm -rf build*
mkdir build; cd build; cmake ..; make -j
cd ..
mkdir build128; cd build128; cmake .. -DWITH_FLOAT128=ON; make -j

for p in 19 18 17 16 15 14
do
    for x in 8 6 4 2 1
    do
        for i in {1..30000}
        do
            n="${x}e$p-10703*$i"
            res128=$(./primecount "$n" --Li);
            res=$(../build/./primecount "$n" --Li);
            if [ "$res" -ne "$res128" ]
            then
                echo "$n: $res128 (float128) != $res (long double)";
                break;
            fi
        done
    done
done

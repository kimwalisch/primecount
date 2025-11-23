#!/bin/bash

# Find the fastest alpha_y tuning factors for primecount.
# Usage:
#   ./find_fastest_alpha_y_gourdon.sh [--start=n] [--stop=n] [-t=threads]
# Description:
#   This script calculates pi(10^n) for start <= n <= stop using different
#   alpha tuning factors and prints out the fastest alphas. The fastest
#   alpha factors will vary slightly for different CPUs. Once a list of
#   fast alpha factors has been generated you can follow the
#   instructions in doc/alpha-factor-tuning.pdf and generate a new
#   alpha tuning function for use in primecount's source code.

# Execute in base directory
if [[ "$(basename $(pwd))" = "scripts" ]]
then
    cd ..
fi

command -v ./primecount >/dev/null 2>/dev/null
if [[ $? -ne 0 ]]
then
    echo "Error: no primecount binary in current directory."
    exit 1
fi

command -v bc >/dev/null 2>/dev/null
if [[ $? -ne 0 ]]
then
    echo "Error: GNU bc is not installed."
    exit 1
fi

start=1
stop=25
seconds=0
repeat=3
threads=$(./primecount 1e18 --Sigma --alpha-z=1 -s | grep threads | cut -d'=' -f2 | cut -d' ' -f2)

for i in "$@"
do
    case $i in
        -t=*)
        threads="${i#*=}"
        shift
        ;;
        --start=*)
        start="${i#*=}"
        shift
        ;;
        --stop=*)
        stop="${i#*=}"
        shift
        ;;
        *)
        echo "Find the fastest alpha_y tuning factors for primecount."
        echo "Usage:"
        echo "  ./find_fastest_alpha_y_gourdon.sh [--start=n] [--stop=n] [-t=threads]"
        echo "Description:"
        echo "  This script calculates pi(10^n) for start <= n <= stop using different"
        echo "  alpha tuning factors and prints out the fastest alphas. The fastest"
        echo "  alpha factors will vary slightly for different CPUs. Once a list of"
        echo "  fast alpha factors has been generated you can follow the"
        echo "  instructions in doc/alpha-factor-tuning.pdf and generate a new"
        echo "  alpha tuning function for use in primecount's source code."
        exit 1
        ;;
    esac
done

# Returns 1 if $1 == $2, else 0
function is_equal
{
    echo $1'=='$2 | bc -l
}

# Returns 1 if $1 < $2, else 0
function is_smaller
{
    echo $1'<'$2 | bc -l
}

# Returns 1 if $1 <= $2, else 0
function is_smaller_equal
{
    echo $1'<='$2 | bc -l
}

# Returns 1 if $1 > $2, else 0
function is_greater
{
    echo $1'>'$2 | bc -l
}

# Floating point calculator
# $1: String containing an arithmetic expression
function calc
{
    echo "scale=3; $1" | bc -l
}

# Given 2 numbers, returns the greatest number 
function maximum
{
    if [[ $(is_greater $1 $2) -eq 1 ]]
    then
        echo $1
    else
        echo $2
    fi
}

# Given 2 numbers, returns the greatest number 
function minimum
{
    if [[ $(is_smaller $1 $2) -eq 1 ]]
    then
        echo $1
    else
        echo $2
    fi
}

# $1: primecount args
function get_primecount_alpha_y
{
    alpha_y=$(./primecount $1 --Sigma --alpha-z=1 -s | grep alpha_y | cut -d'=' -f2 | cut -d' ' -f2)
    echo $alpha_y
}

# $1: primecount args
function get_primecount_seconds
{
    sleep 0.2
    seconds=$(./primecount $1 --time | grep Seconds | cut -d':' -f2 | cut -d' ' -f2)
    echo $seconds
}

# Calculate pi(10^i) for start <= i <= stop
for ((i = start; i <= stop; i++))
do
    alpha_y=$(get_primecount_alpha_y "1e$i")
    fastest_seconds=10^30
    fastest_alpha_y=$alpha_y
    found_fastest=false
    too_fast=false
    iter_count=0
    div_increment=2

    echo ""
    echo "PrimePi(10^$i)"
    echo "==============================================================="

    for div in 2 4 8 16 32 64;
    do
        copy_fastest_alpha_y=$fastest_alpha_y

        for ((j = 0; j < repeat; j++))
        do
            # Benchmark runs too quickly for this small input
            if [[ "$too_fast" = "true" ]]
            then
                break
            fi

            pivot=$(calc "$copy_fastest_alpha_y / $div")
            max_alpha_y=$(calc "$copy_fastest_alpha_y + $pivot")
            new_alpha_y=$(calc "$copy_fastest_alpha_y - $pivot")
            new_alpha_y=$(maximum 1.000 $new_alpha_y)
            increment=$(calc "($max_alpha_y - $new_alpha_y) / $div_increment")
            increment=$(maximum 0.01 $increment)

            while [[ $(is_smaller_equal $new_alpha_y $max_alpha_y) -eq 1 ]]
            do
                seconds=$(get_primecount_seconds "1e$i -t$threads --alpha-z=1 --alpha-y=$new_alpha_y")
                echo "1e$i --threads=$threads --alpha-z=1 --alpha-y=$new_alpha_y, seconds: $seconds"
                iter_count=$(($iter_count + 1))

                if [[ $iter_count -eq 2 ]]
                then
                    old_seconds=$seconds
                fi

                # Benchmark runs too quickly for this small input
                if [[ $(is_equal $seconds 0) -eq 1 ]]
                then
                    too_fast=true
                    break
                fi

                if [[ $(is_smaller $seconds $fastest_seconds) -eq 1 ]]
                then
                    found_fastest=true
                    fastest_alpha_y=$new_alpha_y
                    fastest_seconds=$seconds
                fi

                # Reduce the number of long running primecount benchmarks
                if [[ $(is_greater $repeat 1) -eq 1 ]]
                then
                    if [[ $(is_greater $seconds 20) -eq 1 ]]
                    then
                        repeat=1
                    elif [[ $(is_greater $seconds 3) -eq 1 ]]
                    then
                        repeat=2
                    fi
                fi

                new_alpha_y=$(calc "$new_alpha_y + $increment")
            done
        done

        div_increment=1
    done

    if [[ "$found_fastest" != "true" ]] || \
       [[ "$too_fast" = "true" ]]
    then
        fastest_alpha_y="undef"
        fastest_seconds="undef"
        old_seconds="undef"
    fi

    echo "Result: fastest_alpha_y=$fastest_alpha_y, seconds=$fastest_seconds, old_alpha_y=$alpha_y, old_seconds=$old_seconds"
done

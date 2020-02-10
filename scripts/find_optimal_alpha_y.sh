#!/bin/bash

# Find the optimal alpha_y tuning factors for primecount (Gourdon's algorithm).
# Optimal means that these alpha_y factors use the fewest number of
# instructions. However usually a silghtly smaller alpha_y factor will
# run faster due to CPU cache effects.
#
# Usage:
#   sudo ./find_optimal_alpha_y.sh [--start=n] [--stop=n] [-t=threads]
# Description:
#   This script calculates pi(10^n) for start <= n <= stop using different
#   alpha_y tuning factors and prints out the optimal alpha_y.
#   Once a list of optimal alpha_y factors has been generated you can
#   follow the instructions in doc/alpha-factor-tuning.pdf and generate
#   a new alpha tuning function for use in primecount's source code.

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

command -v bc >/dev/null 2>/dev/null
if [ $? -ne 0 ]
then
    echo "Error: GNU bc is not installed."
    exit 1
fi

start=1
stop=24
repeat=1
threads=$(./primecount 1e18 --Sigma --alpha-z=1 -s | grep threads | cut -f3 -d' ')

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
        echo "Find the optimal alpha_y tuning factors for primecount (Gourdon's algorithm)."
        echo "Optimal means that these alpha_y factors use the fewest number of"
        echo "instructions. However usually a silghtly smaller alpha_y factor will"
        echo "run faster due to CPU cache effects."
        echo ""
        echo "Usage:"
        echo "  sudo ./find_optimal_alpha_y.sh [--start=n] [--stop=n] [-t=threads]"
        echo "Description:"
        echo "  This script calculates pi(10^n) for start <= n <= stop using different"
        echo "  alpha_y tuning factors and prints out the optimal alpha_y."
        echo "  Once a list of optimal alpha_y factors has been generated you can"
        echo "  follow the instructions in doc/alpha-factor-tuning.pdf and generate"
        echo "  a new alpha tuning function for use in primecount's source code."
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
    if [ $(is_greater $1 $2) -eq 1 ]
    then
        echo $1
    else
        echo $2
    fi
}

# Given 2 numbers, returns the greatest number 
function minimum
{
    if [ $(is_smaller $1 $2) -eq 1 ]
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
function get_primecount_instructions
{
    instructions=$(perf stat ./primecount $1 2>&1 | grep instructions | sed "s/^[ \t]*//" | cut -f1 -d' ' | sed "s/\.//g")
    echo $instructions
}

# Calculate pi(10^i) for start <= i <= stop
for ((i = start; i <= stop; i++))
do
    alpha_y=$(get_primecount_alpha_y "1e$i")
    fewest_instructions=100000000000000000
    optimal_alpha_y="1.000"

    # This loops tries to get close (< 25%) to the optimal alpha_y
    for ((j = 0; j < repeat; j++))
    do
        max_alpha_y=$(calc "$alpha_y * 3")
        new_alpha_y=$(calc "$alpha_y / 2")
        new_alpha_y=$(maximum 1.000 $new_alpha_y)
        increment=$(calc "($max_alpha_y - $new_alpha_y) / 6")
        increment=$(maximum 0.1 $increment)

        while [ $(is_smaller_equal $new_alpha_y $max_alpha_y) -eq 1 ]
        do
            instructions=$(get_primecount_instructions "1e$i -t$threads --alpha-y=$new_alpha_y --alpha-z=1")

            if [ $(is_smaller $instructions $fewest_instructions) -eq 1 ]
            then
                optimal_alpha_y=$new_alpha_y
                fewest_instructions=$instructions
            fi

            # Benchmark runs too quickly for this small input
            if [ $(is_equal $instructions 0) -eq 1 ]
            then
                break
            fi

            new_alpha_y=$(calc "$new_alpha_y + $increment")
        done
    done

    if [ $(is_greater $fewest_instructions 0) -eq 1 ]
    then
        # This loops tries to get very close (< 1.6%) to the optimal alpha_y
        for ((j = 0; j < repeat; j++))
        do
            pivot=$(calc "$optimal_alpha_y / 2")
            max_alpha_y=$(calc "$optimal_alpha_y + $pivot")
            new_alpha_y=$(calc "$optimal_alpha_y - $pivot")
            new_alpha_y=$(maximum 1.000 $new_alpha_y)
            increment=$(calc "($max_alpha_y - $new_alpha_y) / 10")
            increment=$(maximum 0.1 $increment)

            while [ $(is_smaller_equal $new_alpha_y $max_alpha_y) -eq 1 ]
            do
                instructions=$(get_primecount_instructions "1e$i -t$threads --alpha-y=$new_alpha_y --alpha-z=1")

                if [ $(is_smaller $instructions $fewest_instructions) -eq 1 ]
                then
                    optimal_alpha_y=$new_alpha_y
                    fewest_instructions=$instructions
                fi

                new_alpha_y=$(calc "$new_alpha_y + $increment")
            done
        done
    fi

    # Print optimal alpha_y found for pi(10^$i)
    printf '%-11s %-12s %-18s %-22s %-20s\n' "pi(10^$i)" "threads=$threads" "old_alpha_y=$alpha_y" "optimal_alpha_y=$optimal_alpha_y" "instructions=$fewest_instructions"
done

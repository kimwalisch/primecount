#!/bin/bash

# Find the optimal alpha tuning factors for primecount.
# Usage:
#   ./find_optimal_alpha_dr.sh [--start=n] [--stop=n] [-t=threads]
# Description:
#   This script calculates pi(10^n) for start <= n <= stop using different
#   alpha tuning factors and prints out the optimal alphas. The optimal
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
instructions=0
repeat=2
threads=$(./primecount 1e18 --Sigma -s | grep threads | cut -d'=' -f2 | cut -d' ' -f2)

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
        echo "Find the optimal alpha tuning factors for primecount."
        echo "Usage:"
        echo "  ./find_optimal_alpha_dr.sh [--start=n] [--stop=n] [-t=threads]"
        echo "Description:"
        echo "  This script calculates pi(10^n) for start <= n <= stop using different"
        echo "  alpha tuning factors and prints out the optimal alphas. The optimal"
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
function get_primecount_alpha
{
    alpha=$(./primecount $1 --S2-trivial -s | grep alpha | cut -d'=' -f2 | cut -d' ' -f2)
    echo $alpha
}

# $1: primecount args
function get_primecount_instructions
{
    instructions=$(perf stat ./primecount $1 -d 2>&1 | grep instructions | sed "s/^[ \t]*//" | cut -f1 -d' ' | sed 's/[^0-9]*//g')
    echo $instructions
}

# Calculate pi(10^i) for start <= i <= stop
for ((i = start; i <= stop; i++))
do
    alpha=$(get_primecount_alpha "1e$i")
    optimal_instructions=10^30
    optimal_alpha=$alpha
    found_optimal=false
    too_fast=false
    iter_count=0
    div_increment=2

    echo ""
    echo "PrimePi(10^$i)"
    echo "==============================================================="

    for div in 2 4 8 16 32;
    do
        copy_optimal_alpha=$optimal_alpha

        for ((j = 0; j < repeat; j++))
        do
            # Benchmark runs too quickly for this small input
            if [[ "$too_fast" = "true" ]]
            then
                break
            fi

            pivot=$(calc "$copy_optimal_alpha / $div")
            max_alpha=$(calc "$copy_optimal_alpha + $pivot")
            new_alpha=$(calc "$copy_optimal_alpha - $pivot")
            new_alpha=$(maximum 1.000 $new_alpha)
            increment=$(calc "($max_alpha - $new_alpha) / $div_increment")
            increment=$(maximum 0.1 $increment)

            if [[ $(is_smaller $increment 0.1) -eq 1 ]]
            then
                break
            fi

            while [[ $(is_smaller_equal $new_alpha $max_alpha) -eq 1 ]]
            do
                instructions=$(get_primecount_instructions "1e$i -t$threads --alpha=$new_alpha")
                echo "1e$i --threads=$threads --alpha=$new_alpha, instructions: $instructions"
                iter_count=$(($iter_count + 1))

                if [[ $iter_count -eq 2 ]]
                then
                    old_instructions=$instructions
                fi

                if [[ $(is_smaller $instructions $optimal_instructions) -eq 1 ]]
                then
                    found_optimal=true
                    optimal_alpha=$new_alpha
                    optimal_instructions=$instructions
                fi

                # Reduce the number of long running primecount benchmarks
                if [[ $(is_greater $i 15) -eq 1 ]]
                then
                    repeat=1
                fi

                new_alpha=$(calc "$new_alpha + $increment")
            done
        done

        div_increment=1
    done

    if [[ "$found_optimal" != "true" ]] || \
       [[ "$too_fast" = "true" ]]
    then
        optimal_alpha="undef"
        optimal_instructions="undef"
        old_instructions="undef"
    fi

    echo "Result: optimal_alpha=$optimal_alpha, instructions=$optimal_instructions, old_alpha=$alpha, old_instructions=$old_instructions"
done

#!/bin/bash

# Usage: ./generate_html_benchmarks.sh
# Generates an html table with pi(x) benchmark results
# for use in the README.md file.
# HTML is written to a timestamped benchmark file in the current directory.
# Benchmark progress is printed to stdout.

export LC_ALL=C

# Progress and status messages are written to stdout, while HTML is written
# to a temporary file and renamed after all benchmarks finish successfully.

if [[ ! -x ./primecount ]]
then
    echo "Error: no primecount binary in current directory." >&2
    exit 1
fi

html_file=$(mktemp .benchmark.tmp.XXXXXX) || exit 1
trap 'rm -f "$html_file"' EXIT
exec 3>"$html_file" || exit 1

function html
{
    if ! printf '%s\n' "$1" >&3
    then
        echo "Error: failed to write benchmark HTML." >&2
        exit 1
    fi
}

function format_integer
{
    local value=$1
    awk -v value="$value" 'BEGIN {
        len = length(value)
        for (i = 1; i <= len; i++) {
            printf "%s", substr(value, i, 1)
            if (i < len && (len - i) % 3 == 0)
                printf ","
        }
    }'
}

function round_seconds
{
    local secs=$1
    awk -v secs="$secs" 'BEGIN {
        split(secs, parts, ".")
        integer = parts[1] + 0
        fraction = parts[2] "000"
        hundredths = substr(fraction, 1, 2) + 0
        if (substr(fraction, 3, 1) >= 5)
            hundredths++
        if (hundredths >= 100) {
            integer++
            hundredths -= 100
        }
        printf "%d.%02d", integer, hundredths
    }'
}

function pretty_print_secs
{
    local rounded integer fraction
    rounded=$(round_seconds "$1")
    integer=${rounded%.*}
    fraction=${rounded#*.}
    echo "$(format_integer "$integer").$fraction"
}

function get_num_runs
{
    local option=$1
    local secs=$2

    if [[ $option == "--gourdon" ]]
    then
        echo 5
    elif awk -v secs="$secs" 'BEGIN { exit !(secs < 3) }'
    then
        echo 5
    elif awk -v secs="$secs" 'BEGIN { exit !(secs <= 10) }'
    then
        echo 3
    elif awk -v secs="$secs" 'BEGIN { exit !(secs <= 100) }'
    then
        echo 2
    else
        echo 1
    fi
}

function run_benchmark
{
    local x=$1
    local option=$2
    local output

    sleep 1
    if ! output=$(./primecount "$x" "$option" --time)
    then
        echo "Error: ./primecount $x $option --time failed." >&2
        return 1
    fi

    run_prime_count=$(echo "$output" | awk '/^[0-9]+$/ { print; exit }')
    run_seconds=$(echo "$output" | awk '$1 == "Seconds:" { print $2; exit }')

    if [[ ! $run_prime_count =~ ^[0-9]+$ || ! $run_seconds =~ ^[0-9]+([.][0-9]+)?$ ]]
    then
        echo "Error: could not parse output of ./primecount $x $option --time." >&2
        return 1
    fi
}

function primecount_benchmark
{
    local x=$1
    local option=$2
    local algorithm=$3
    local run

    run_benchmark "$x" "$option" || return 1
    benchmark_prime_count=$run_prime_count
    benchmark_seconds=$run_seconds
    benchmark_runs=$(get_num_runs "$option" "$run_seconds")

    for ((run = 2; run <= benchmark_runs; run++))
    do
        run_benchmark "$x" "$option" || return 1

        if [[ $run_prime_count != "$benchmark_prime_count" ]]
        then
            echo "Error: inconsistent prime counts for $x $option." >&2
            return 1
        fi

        if awk -v a="$run_seconds" -v b="$benchmark_seconds" 'BEGIN { exit !(a < b) }'
        then
            benchmark_seconds=$run_seconds
        fi
    done

    printf "%-24s %-8s %4d %12ss\n" "$algorithm" "$x" "$benchmark_runs" "$benchmark_seconds"
}

function check_prime_count
{
    if [[ $benchmark_prime_count != "$prime_count" ]]
    then
        echo "Error: inconsistent prime counts for $x." >&2
        exit 1
    fi
}

printf "%-24s %-8s %4s %13s\n" "Algorithm" "x" "Runs" "Fastest"
printf "%-24s %-8s %4s %13s\n" "------------------------" "--------" "----" "-------------"

html "<table>"
html '  <tr align="center">'
html "    <td><b>x</b></td>"
html "    <td><b>Prime Count</b></td>"
html "    <td><b>Legendre</b></td>"
html "    <td><b>Meissel</b></td>"
html "    <td><b>Lagarias<br/>Miller<br/>Odlyzko</b></td>"
html "    <td><b>Deleglise<br/>Rivat</b></td>"
html "    <td><b>Gourdon</b></td>"
html "  </tr>"

for i in {10..18}
do
    x=1e$i
    ((i != 10)) && printf '\n'

    primecount_benchmark "$x" --legendre "Legendre" || exit 1
    prime_count=$benchmark_prime_count
    legendre=$benchmark_seconds

    primecount_benchmark "$x" --meissel "Meissel" || exit 1
    check_prime_count
    meissel=$benchmark_seconds

    primecount_benchmark "$x" --lmo "Lagarias-Miller-Odlyzko" || exit 1
    check_prime_count
    lmo=$benchmark_seconds

    primecount_benchmark "$x" --deleglise-rivat "Deleglise-Rivat" || exit 1
    check_prime_count
    deleglise_rivat=$benchmark_seconds

    primecount_benchmark "$x" --gourdon "Gourdon" || exit 1
    check_prime_count
    gourdon=$benchmark_seconds

    html '  <tr align="right">'
    html "    <td>10<sup>$i</sup></td>"
    html "    <td>$(format_integer "$prime_count")</td>"
    html "    <td>$(pretty_print_secs "$legendre")s</td>"
    html "    <td>$(pretty_print_secs "$meissel")s</td>"
    html "    <td>$(pretty_print_secs "$lmo")s</td>"
    html "    <td>$(pretty_print_secs "$deleglise_rivat")s</td>"
    html "    <td>$(pretty_print_secs "$gourdon")s</td>"
    html "  </tr>"
done

for i in {19..22}
do
    x=1e$i
    printf '\n'

    primecount_benchmark "$x" --deleglise-rivat "Deleglise-Rivat" || exit 1
    prime_count=$benchmark_prime_count
    deleglise_rivat=$benchmark_seconds

    primecount_benchmark "$x" --gourdon "Gourdon" || exit 1
    check_prime_count
    gourdon=$benchmark_seconds

    html '  <tr align="right">'
    html "    <td>10<sup>$i</sup></td>"
    html "    <td>$(format_integer "$prime_count")</td>"
    html "    <td>NaN</td>"
    html "    <td>NaN</td>"
    html "    <td>NaN</td>"
    html "    <td>$(pretty_print_secs "$deleglise_rivat")s</td>"
    html "    <td>$(pretty_print_secs "$gourdon")s</td>"
    html "  </tr>"
done

html "</table>"
exec 3>&-

while true
do
    benchmark_file="benchmark-$(date '+%Y-%m-%d_%H-%M-%S').txt"
    [[ ! -e $benchmark_file ]] && break
    sleep 1
done

if ! mv "$html_file" "$benchmark_file"
then
    echo "Error: failed to create $benchmark_file." >&2
    exit 1
fi
echo "Benchmark table written to $PWD/$benchmark_file"

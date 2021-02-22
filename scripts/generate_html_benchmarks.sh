#!/bin/bash

# Usage: ./generate_html_benchmarks.sh
# Generates an html table with pi(x) benchmark results
# for use in the README.md file.

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

# Returns 1 if $1 < $2, else 0
function is_smaller
{
    echo $1'<'$2 | bc -l
}

function round()
{
    echo $(printf %.$2f $(echo "scale=$2;(((10^$2)*$1)+0.5)/(10^$2)" | bc))
}

function get_prime_count()
{
    echo $(./primecount $1 | awk '{ len=length($0); res=""; for (i=0;i<=len;i++) { res=substr($0,len-i+1,1) res; if (i > 0 && i < len && i % 3 == 0) { res = "," res } }; print res }')
}

function primecount_benchmark()
{
    sleep 1;
    secs1=$(./primecount $1 $2 --time | grep Seconds | cut -f2 -d' ')
    sleep 1;
    secs2=$(./primecount $1 $2 --time | grep Seconds | cut -f2 -d' ')

    if [[ $(is_smaller $secs1 $secs2) -eq 1 ]]
    then
        echo $secs1
    else
        echo $secs2
    fi
}

echo "<table>"
echo '  <tr align="center">'
echo "    <td><b>x</b></td>"
echo "    <td><b>Prime Count</b></td>"
echo "    <td><b>Legendre</b></td>"
echo "    <td><b>Meissel</b></td>"
echo "    <td><b>Lagarias<br/>Miller<br/>Odlyzko</b></td>"
echo "    <td><b>Deleglise<br/>Rivat</b></td>"
echo "    <td><b>Gourdon</b></td>"
echo "  </tr>"

for i in {10..18};
do
    echo '  <tr align="right">'
    echo "    <td>10<sup>$i</sup></td>"
    echo "    <td>$(get_prime_count 1e$i)</td>"
    echo "    <td>$(round $(primecount_benchmark 1e$i --legendre) 2)s</td>"
    echo "    <td>$(round $(primecount_benchmark 1e$i --meissel) 2)s</td>"
    echo "    <td>$(round $(primecount_benchmark 1e$i --lmo) 2)s</td>"
    echo "    <td>$(round $(primecount_benchmark 1e$i --deleglise-rivat) 2)s</td>"
    echo "    <td>$(round $(primecount_benchmark 1e$i --gourdon) 2)s</td>"
    echo "  </tr>"
done

for i in {19..22};
do
    echo '  <tr align="right">'
    echo "    <td>10<sup>$i</sup></td>"
    echo "    <td>$(get_prime_count 1e$i)</td>"
    echo "    <td>NaN</td>"
    echo "    <td>NaN</td>"
    echo "    <td>NaN</td>"
    echo "    <td>$(round $(primecount_benchmark 1e$i --deleglise-rivat) 2)s</td>"
    echo "    <td>$(round $(primecount_benchmark 1e$i --gourdon) 2)s</td>"
    echo "  </tr>"
done

echo "</table>"

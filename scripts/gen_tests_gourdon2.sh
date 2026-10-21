#!/bin/bash

for formula in "B" "Sigma"
do
    echo "=== $formula ======================================================="
    echo ""

    for i in {1..13};
    do
        res=$(./primecount 1e$i --$formula)
        x=$(./primecount 1e$i --$formula -s | grep '^x =' | cut -f3 -d' ')
        y=$(./primecount 1e$i --$formula -s | grep '^y =' | cut -f3 -d' ')

        verify=$(./primecount 1e$i -g -s | grep "$formula = $res\$")
        if [[ -z "$verify" ]] || [[ "$(./primecount 1e$i -g)" != "$(./primecount 1e$i -m)" ]]
        then
            echo ""
            echo "Error!"
            exit 1
        fi

        echo "{ ${x}, ${y}, ${res} },"

        res=$(./primecount 1e$i --$formula --alpha=10000000)
        x=$(./primecount 1e$i --$formula --alpha=10000000 -s | grep '^x =' | cut -f3 -d' ')
        y=$(./primecount 1e$i --$formula --alpha=10000000 -s | grep '^y =' | cut -f3 -d' ')

        verify=$(./primecount 1e$i --alpha=10000000 -g -s | grep "$formula = $res\$")
        if [[ -z "$verify" ]] || [[ "$(./primecount 1e$i -g --alpha=10000000)" != "$(./primecount 1e$i -m)" ]]
        then
            echo ""
            echo "Error!"
            exit 1
        fi

        echo "{ ${x}, ${y}, ${res} },"

        res=$(./primecount 1e$i --$formula --alpha=1)
        x=$(./primecount 1e$i --$formula --alpha=1 -s | grep '^x =' | cut -f3 -d' ')
        y=$(./primecount 1e$i --$formula --alpha=1 -s | grep '^y =' | cut -f3 -d' ')

        verify=$(./primecount 1e$i --alpha=1 -g -s | grep "$formula = $res\$")
        if [[ -z "$verify" ]] || [[ "$(./primecount 1e$i -g --alpha=1)" != "$(./primecount 1e$i -m)" ]]
        then
            echo ""
            echo "Error!"
            exit 1
        fi

        echo "{ ${x}, ${y}, ${res} },"
    done

    # A few larger computations using default alpha
    for i in 14 15;
    do
        res=$(./primecount 1e$i --$formula)
        x=$(./primecount 1e$i --$formula -s | grep '^x =' | cut -f3 -d' ')
        y=$(./primecount 1e$i --$formula -s | grep '^y =' | cut -f3 -d' ')

        verify=$(./primecount 1e$i -g -s | grep "= $res\$")
        if [[ -z "$verify" ]] || [[ "$(./primecount 1e$i -g)" != "$(./primecount 1e$i -m)" ]]
        then
            echo ""
            echo "Error!"
            exit 1
        fi

        echo "{ ${x}, ${y}, ${res} },"
    done

    echo ""
done

echo "Suceess!"

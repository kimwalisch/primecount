#!/bin/bash

for formula in "S2-trivial" "S2-easy" "S2-hard"
do
    echo "=== $formula ======================================================="
    echo ""

    for i in {1..11};
    do
        res=$(./primecount 1e$i --$formula)
        x=$(./primecount 1e$i --$formula -s | grep '^x =' | cut -f3 -d' ')
        y=$(./primecount 1e$i --$formula -s | grep '^y =' | cut -f3 -d' ')
        z=$(./primecount 1e$i --$formula -s | grep '^z =' | cut -f3 -d' ')
        c=$(./primecount 1e$i --$formula -s | grep '^c =' | cut -f3 -d' ')

        verify=$(./primecount 1e$i -d -s | grep "${formula//-/_} = $res\$")
        if [[ -z "$verify" ]] || [[ "$(./primecount 1e$i -d)" != "$(./primecount 1e$i -m)" ]]
        then
            echo ""
            echo "Error!"
            exit 1
        fi

        echo "{ ${x}LL, ${y}, ${z}, ${c}, ${res}LL },"

        res=$(./primecount 1e$i --$formula --alpha=1)
        x=$(./primecount 1e$i --$formula --alpha=1 -s | grep '^x =' | cut -f3 -d' ')
        y=$(./primecount 1e$i --$formula --alpha=1 -s | grep '^y =' | cut -f3 -d' ')
        z=$(./primecount 1e$i --$formula --alpha=1 -s | grep '^z =' | cut -f3 -d' ')
        c=$(./primecount 1e$i --$formula --alpha=1 -s | grep '^c =' | cut -f3 -d' ')

        verify=$(./primecount 1e$i --alpha=1 -d -s | grep "${formula//-/_} = $res\$")
        if [[ -z "$verify" ]] || [[ "$(./primecount 1e$i -d --alpha=1)" != "$(./primecount 1e$i -m)" ]]
        then
            echo ""
            echo "Error!"
            exit 1
        fi

        echo "{ ${x}LL, ${y}, ${z}, ${c}, ${res}LL },"

        res=$(./primecount 1e$i --$formula --alpha=10000000)
        x=$(./primecount 1e$i --$formula --alpha=10000000 -s | grep '^x =' | cut -f3 -d' ')
        y=$(./primecount 1e$i --$formula --alpha=10000000 -s | grep '^y =' | cut -f3 -d' ')
        z=$(./primecount 1e$i --$formula --alpha=10000000 -s | grep '^z =' | cut -f3 -d' ')
        c=$(./primecount 1e$i --$formula --alpha=10000000 -s | grep '^c =' | cut -f3 -d' ')

        verify=$(./primecount 1e$i --alpha=10000000 -d -s | grep "${formula//-/_} = $res\$")
        if [[ -z "$verify" ]] || [[ "$(./primecount 1e$i -d --alpha=10000000)" != "$(./primecount 1e$i -m)" ]]
        then
            echo ""
            echo "Error!"
            exit 1
        fi

        echo "{ ${x}LL, ${y}, ${z}, ${c}, ${res}LL },"
    done

    # A few larger computations using default alpha
    for i in 13 14;
    do
        res=$(./primecount 1e$i --$formula)
        x=$(./primecount 1e$i --$formula -s | grep '^x =' | cut -f3 -d' ')
        y=$(./primecount 1e$i --$formula -s | grep '^y =' | cut -f3 -d' ')
        z=$(./primecount 1e$i --$formula -s | grep '^z =' | cut -f3 -d' ')
        c=$(./primecount 1e$i --$formula -s | grep '^c =' | cut -f3 -d' ')

        verify=$(./primecount 1e$i -d -s | grep "= $res\$")
        if [[ -z "$verify" ]] || [[ "$(./primecount 1e$i -d)" != "$(./primecount 1e$i -m)" ]]
        then
            echo ""
            echo "Error!"
            exit 1
        fi

        echo "{ ${x}LL, ${y}, ${z}, ${c}, ${res}LL },"
    done

    echo ""
done

echo "Suceess!"

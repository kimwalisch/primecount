#!/bin/bash

for formula in "D" "AC" "Phi0"
do
    echo "=== $formula ======================================================="
    echo ""

    # When using primecount x -g -s the result of the
    # AC algorithm is printed as "A + C = res"  
    if [[ "$formula" == "AC" ]]
    then
        formula_label="A + C"
    else
        formula_label=$formula
    fi

    for i in {1..12};
    do
        res=$(./primecount 1e$i --$formula)
        x=$(./primecount 1e$i --$formula -s | grep '^x =' | cut -f3 -d' ')
        y=$(./primecount 1e$i --$formula -s | grep '^y =' | cut -f3 -d' ')
        z=$(./primecount 1e$i --$formula -s | grep '^z =' | cut -f3 -d' ')
        k=$(./primecount 1e$i --$formula -s | grep '^k =' | cut -f3 -d' ')

        verify=$(./primecount 1e$i -g -s | grep "$formula_label = $res\$")
        if [[ -z "$verify" ]] || [[ "$(./primecount 1e$i -g)" != "$(./primecount 1e$i -m)" ]]
        then
            echo ""
            echo "Error!"
            exit 1
        fi

        echo "{ ${x}LL, ${y}, ${z}, ${k}, ${res}LL },"

        res=$(./primecount 1e$i --$formula --alpha-y=10000000)
        x=$(./primecount 1e$i --$formula --alpha-y=10000000 -s | grep '^x =' | cut -f3 -d' ')
        y=$(./primecount 1e$i --$formula --alpha-y=10000000 -s | grep '^y =' | cut -f3 -d' ')
        z=$(./primecount 1e$i --$formula --alpha-y=10000000 -s | grep '^z =' | cut -f3 -d' ')
        k=$(./primecount 1e$i --$formula --alpha-y=10000000 -s | grep '^k =' | cut -f3 -d' ')

        verify=$(./primecount 1e$i --alpha-y=10000000 -g -s | grep "$formula_label = $res\$")
        if [[ -z "$verify" ]] || [[ "$(./primecount 1e$i -g --alpha-y=10000000)" != "$(./primecount 1e$i -m)" ]]
        then
            echo ""
            echo "Error!"
            exit 1
        fi

        echo "{ ${x}LL, ${y}, ${z}, ${k}, ${res}LL },"

        res=$(./primecount 1e$i --$formula --alpha-z=10000000)
        x=$(./primecount 1e$i --$formula --alpha-z=10000000 -s | grep '^x =' | cut -f3 -d' ')
        y=$(./primecount 1e$i --$formula --alpha-z=10000000 -s | grep '^y =' | cut -f3 -d' ')
        z=$(./primecount 1e$i --$formula --alpha-z=10000000 -s | grep '^z =' | cut -f3 -d' ')
        k=$(./primecount 1e$i --$formula --alpha-z=10000000 -s | grep '^k =' | cut -f3 -d' ')

        verify=$(./primecount 1e$i --alpha-z=10000000 -g -s | grep "$formula_label = $res\$")
        if [[ -z "$verify" ]] || [[ "$(./primecount 1e$i -g --alpha-z=10000000)" != "$(./primecount 1e$i -m)" ]]
        then
            echo ""
            echo "Error!"
            exit 1
        fi

        echo "{ ${x}LL, ${y}, ${z}, ${k}, ${res}LL },"
    done

    # A few larger computations using default alpha
    for i in 13 14 15;
    do
        res=$(./primecount 1e$i --$formula)
        x=$(./primecount 1e$i --$formula -s | grep '^x =' | cut -f3 -d' ')
        y=$(./primecount 1e$i --$formula -s | grep '^y =' | cut -f3 -d' ')
        z=$(./primecount 1e$i --$formula -s | grep '^z =' | cut -f3 -d' ')
        k=$(./primecount 1e$i --$formula -s | grep '^k =' | cut -f3 -d' ')

        verify=$(./primecount 1e$i -g -s | grep "$formula_label = $res\$")
        if [[ -z "$verify" ]] || [[ "$(./primecount 1e$i -g)" != "$(./primecount 1e$i -m)" ]]
        then
            echo ""
            echo "Error!"
            exit 1
        fi

        echo "{ ${x}LL, ${y}, ${z}, ${k}, ${res}LL },"
    done

    echo ""
done

echo "Suceess!"

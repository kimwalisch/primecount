#!/bin/sh
# Script that prints the CPUs cache line size in bytes
# Usage: ./cache_line_size.sh

command -v getconf >/dev/null 2>/dev/null
if [ $? -eq 0 ]
then
    CACHE_LINE_SIZE=$(getconf LEVEL1_DCACHE_LINESIZE 2>/dev/null)
fi

if test "x$CACHE_LINE_SIZE" = "x" || \
   test "$CACHE_LINE_SIZE" = "0"
then
    CACHE_LINE_SIZE=$(cat /sys/devices/system/cpu/cpu0/cache/index0/coherency_line_size 2>/dev/null)

    if test "x$CACHE_LINE_SIZE" = "x"
    then
        command -v sysctl >/dev/null 2>/dev/null
        if [ $? -eq 0 ]
        then
            CACHE_LINE_SIZE=$(sysctl hw.cachelinesize 2>/dev/null | sed -e 's/^.* //')
        fi
    fi
fi

if test "x$CACHE_LINE_SIZE" != "x"
then
    # Check detected cache line size
    if [ $CACHE_LINE_SIZE -ge 16 2>/dev/null ] && \
       [ $CACHE_LINE_SIZE -le 8192 2>/dev/null ]
    then
        # CACHE_LINE_SIZE must be power of 2
        old=$CACHE_LINE_SIZE
        CACHE_LINE_SIZE=1
        while [ $CACHE_LINE_SIZE -lt $old ]
        do
            CACHE_LINE_SIZE=$(expr $CACHE_LINE_SIZE '*' 2)
        done

        echo "$CACHE_LINE_SIZE"
        exit 0
    fi
fi

# Failed to detect cache line size
exit 1

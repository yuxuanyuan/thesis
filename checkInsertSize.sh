#!/bin/bash

head -n 20000 $1 | awk '{ if ($9 > 0) { N+=1; S+=$9; S2+=$9*$9 }} END { M=S/N; print "n="N", mean="int (M)", stdev="int (sqrt ((S2-M*M*N)/(N-1)))}' 2>/dev/null
if [[ $? != 0 ]]; then
    echo "n=XXXX, mean=NA, stdev=NA"
fi

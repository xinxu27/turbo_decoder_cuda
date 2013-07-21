#!/bin/bash

max=25
for i in *Max*
do
    block_num=$(cut -d '-' -f1 $i)
    echo $block_num
    for (( iter=6; iter<=$max; iter=iter+1 ))
    do
        echo hello
        #cat $i | grep 'Ber' | cut -d '=' -f2 | awk 'NR%25 == "'$iter'"{print $0}' > Iters${iter}_${i}
    done
done

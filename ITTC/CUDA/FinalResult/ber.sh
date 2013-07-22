#!/bin/bash

max=25
for i in *Max*
do
    block_num=$(echo $i)
    Block=${block_num%%_*txt}
    BlockNum=$((${Block}*4))
    for (( iter=6; iter<=$max; iter=iter+1 ))
    do
        cat $i | grep 'Ber' | cut -d '=' -f2 | awk 'NR%25 == "'$iter'"{print $0}' > Iters${iter}_${BlockNum}
    done
done


#! /usr/bin/env bash

echo "Compiling..."
rm *.o
rm ./poly-eval
make
echo "done. Creating files..."

pmax=12
p=2
while [ $p -le $pmax ] do
    echo "" > n-100k_p-$p.txt
    let p=p+1
done

seq 100000 | shuf > batch_consts.txt

echo "" > batch_vars.txt
i=0
while [ $i -lt 10 ] do
    RANDOM=$$
    r1=$((${RANDOM}%98+1))
    printf -v r1 "0.%.2d" $r1
    echo "$r1" >> batch_vars.txt
    let i=i+1
done

echo "done. Running..."

i=0
p=2
while [ $p -le $pmax]
    while [ $i -lt 10 ] do
        mpirun -np $p batch_consts.txt batch_vars.txt >> n-100k_p-$p.txt
        let i=i+1
    done

echo "done. Exiting."


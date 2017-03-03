#! /usr/bin/env bash

echo "Compiling..."
rm *.o
rm ./poly-eval
make
echo "Creating files..."

pmax=24
p=2
while [ $p -le $pmax ]; do
    echo "" > n-100k_p-${p}.txt
    (( p++ ))
done

echo "100000" > batch_consts.txt
seq 100000 | shuf >> batch_consts.txt

echo "10" > batch_vars.txt
i=0
RANDOM=$$
while [ $i -lt 10 ]; do
    r1=$((${RANDOM}%98+1))
    printf -v r1 "0.%.2d" $r1
    echo "$r1" >> batch_vars.txt
    (( i++ ))
done

echo "Running..."

MPIRUN=/usr/lib64/openmpi/bin/mpirun

j=0
p=2
while [ $p -le $pmax ]; do
    while [ $j -lt 10 ]; do
        $MPIRUN -np $p --hostfile $PBS_NODEFILE ./poly-eval batch_consts.txt batch_vars.txt >> n-100k_p-$p.txt
        let j=j+1
    done
    let p=p+1
    j=0
done

echo "Exiting."


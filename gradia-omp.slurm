#!/bin/bash
#SBATCH --time=0:05:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
gcc -O2 -std=c99 -o randomgraph randomgraph.c
gcc -O2 -std=c99 -fopenmp -o "$USER-gradia" "$USER-gradia.c"
for t in 12; do
    echo OMP_NUM_THREADS=$t
    ./randomgraph 3000 0.1 999 20190824 | OMP_NUM_THREADS=$t "./$USER-gradia"
done

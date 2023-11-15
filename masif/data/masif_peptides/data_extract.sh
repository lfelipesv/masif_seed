#!/bin/bash
#SBATCH --nodes 1
#SBATCH --partition=serial
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 16384
#SBATCH --time 03:00:00
#SBATCH --array=1-5
#SBATCH --output=exelogs/_masif_precompute.%A_%a.out
#SBATCH --error=exelogs/_masif_precompute.%A_%a.err

N=20
i=1
while read p; do
    echo $i
    FIELD1=$(echo $p| cut -d" " -f1)
    PDBID=$(echo $FIELD1| cut -d"_" -f1)
    CHAIN1=$(echo $FIELD1| cut -d"_" -f2)
    ./01-data_copy_and_extract.sh $PDBID\_$CHAIN1 &
    if [ $(( i % N )) -eq 0 ]; then
        wait
    fi
    i=$((i+1))
    echo $i
done < lists/ab_list_synthetic_5000.txt

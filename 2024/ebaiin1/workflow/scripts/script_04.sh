#!/bin/bash

mkdir -p fastqc

module load fastqc/0.11.9

for SAMPLE in FILE1 FILE2 FILE3
do
    echo ">>> Processing $SAMPLE"
    fastqc --outdir fastqc data/${SAMPLE}.fastq.gz 2> fastqc/${SAMPLE}.log
done


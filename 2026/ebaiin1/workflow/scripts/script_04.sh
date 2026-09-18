#!/bin/bash

mkdir -p fastqc

module load fastqc/0.11.9

for SAMPLE in FILE1 FILE2 FILE3
do
    echo ">>> Processing $SAMPLE"
    fastqc --outdir fastqc data/${SAMPLE}_R1.fastq.gz 2> fastqc/${SAMPLE}_R1.log
    fastqc --outdir fastqc data/${SAMPLE}_R2.fastq.gz 2> fastqc/${SAMPLE}_R2.log
done

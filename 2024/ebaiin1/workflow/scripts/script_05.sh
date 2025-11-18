#!/bin/bash

mkdir -p fastqc

module load fastqc/0.11.9

for SAMPLE in $(ls data/ | sed 's/.fastq.gz//')
do
    echo ">>> Processing $SAMPLE"
    fastqc --outdir fastqc data/${SAMPLE}_R1.fastq.gz &> fastqc/${SAMPLE}_R1.log
    fastqc --outdir fastqc data/${SAMPLE}_R2.fastq.gz &> fastqc/${SAMPLE}_R2.log
done






#!/bin/bash

mkdir -p fastqc

module load fastqc/0.11.9

for SAMPLE in $(ls --color=none data/ | sed 's/.fastq.gz//')
do
    echo ">>> Processing $SAMPLE"
    fastqc --outdir fastqc data/${SAMPLE}.fastq.gz &> fastqc/${SAMPLE}.log
done







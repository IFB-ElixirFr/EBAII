#!/bin/bash

mkdir -p fastqc

module load fastqc/0.11.9

SAMPLE=FILE1
echo ">>> Processing $SAMPLE"
fastqc --outdir fastqc data/${SAMPLE}.fastq.gz

SAMPLE=FILE2
echo ">>> Processing $SAMPLE"
fastqc --outdir fastqc data/${SAMPLE}.fastq.gz

SAMPLE=FILE3
echo ">>> Processing $SAMPLE"
fastqc --outdir fastqc data/${SAMPLE}.fastq.gz


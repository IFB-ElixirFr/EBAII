#!/bin/bash

mkdir -p fastqc

module load fastqc/0.11.9

SAMPLE=FILE1
echo ">>> Processing $SAMPLE"
fastqc --outdir fastqc data/${SAMPLE}.fastq.gz 2> fastqc/${SAMPLE}.log

SAMPLE=FILE2
echo ">>> Processing $SAMPLE"
fastqc --outdir fastqc data/${SAMPLE}.fastq.gz 2> fastqc/${SAMPLE}.log

SAMPLE=FILE3
echo ">>> Processing $SAMPLE"
fastqc --outdir fastqc data/${SAMPLE}.fastq.gz 2> fastqc/${SAMPLE}.log



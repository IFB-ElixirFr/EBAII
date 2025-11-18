#!/bin/bash

mkdir -p fastqc

module load fastqc/0.11.9

SAMPLE=FILE1
echo ">>> Processing $SAMPLE"
fastqc --outdir fastqc data/${SAMPLE}_R1.fastq.gz
fastqc --outdir fastqc data/${SAMPLE}_R2.fastq.gz

SAMPLE=FILE2
echo ">>> Processing $SAMPLE"
fastqc --outdir fastqc data/${SAMPLE}_R1.fastq.gz
fastqc --outdir fastqc data/${SAMPLE}_R2.fastq.gz

SAMPLE=FILE3
echo ">>> Processing $SAMPLE"
fastqc --outdir fastqc data/${SAMPLE}_R1.fastq.gz
fastqc --outdir fastqc data/${SAMPLE}_R2.fastq.gz
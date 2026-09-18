#!/bin/bash

mkdir -p fastqc

module load fastqc/0.11.9

SAMPLE=FILE1
echo ">>> Processing $SAMPLE"
fastqc --outdir fastqc data/${SAMPLE}_R1.fastq.gz 2> fastqc/${SAMPLE}_R1.log
fastqc --outdir fastqc data/${SAMPLE}_R2.fastq.gz 2> fastqc/${SAMPLE}_R2.log

SAMPLE=FILE2
echo ">>> Processing $SAMPLE"
fastqc --outdir fastqc data/${SAMPLE}_R1.fastq.gz 2> fastqc/${SAMPLE}_R1.log
fastqc --outdir fastqc data/${SAMPLE}_R2.fastq.gz 2> fastqc/${SAMPLE}_R2.log

SAMPLE=FILE3
echo ">>> Processing $SAMPLE"
fastqc --outdir fastqc data/${SAMPLE}_R1.fastq.gz 2> fastqc/${SAMPLE}_R1.log
fastqc --outdir fastqc data/${SAMPLE}_R2.fastq.gz 2> fastqc/${SAMPLE}_R2.log
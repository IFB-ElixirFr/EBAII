#!/bin/bash


mkdir -p fastqc

module load fastqc/0.11.9

fastqc --outdir fastqc data/FILE1_R1.fastq.gz
fastqc --outdir fastqc data/FILE1_R2.fastq.gz
#!/bin/bash
 
mkdir -p fastqc

module load fastqc/0.11.9


fastqc --outdir fastqc data/FILE1.fastq.gz

fastqc --outdir fastqc data/FILE2.fastq.gz

fastqc --outdir fastqc data/FILE3.fastq.gz


#!/bin/bash

#SBATCH --account=2538_eb3i_n1_2025
#SBATCH --job-name=bwa_align
#SBATCH --cpus-per-task=8
#SBATCH --mem=20GB
#SBATCH --array=1-3

module load bwa/0.7.17
module load samtools/1.18
  
OUT_DIR=$PWD"/bwa_res"
mkdir -p $OUT_DIR

GENOME_INDEX=$PWD"/TAIR10"

SAMPLE=$(ls data/*R1* | sed 's/.fastq.gz//' | awk "NR==${SLURM_ARRAY_TASK_ID}")
BASENAME=$(basename ${SAMPLE/_R1/})

  
srun bwa mem -t 2 ${GENOME_INDEX}  ${SAMPLE}.fastq.gz  ${SAMPLE/R1/R2}.fastq.gz > $OUT_DIR/${BASENAME}.unsorted.sam
srun samtools view -@ 2 -q 10 -b $OUT_DIR/${BASENAME}.unsorted.sam | samtools sort -@ 2 - -o $OUT_DIR/${BASENAME}.sorted.bam 

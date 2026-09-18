#!/bin/bash

#SBATCH --account=2538_eb3i_n1_2025
#SBATCH --job-name=minimap_align
#SBATCH --cpus-per-task=8
#SBATCH --mem=20GB
#SBATCH --array=1-3

module load minimap2/2.28
module load samtools/1.18
  
OUT_DIR=$PWD"/minimap_res"
mkdir -p $OUT_DIR

GENOME="/shared/bank/arabidopsis_thaliana/TAIR10.1/fasta/GCF_000001735.4_TAIR10.1_genomic.fna"

SAMPLE=$(ls data/*R1* | sed 's/.fastq.gz//' | awk "NR==${SLURM_ARRAY_TASK_ID}")
BASENAME=$(basename ${SAMPLE/_R1/})

  
srun minimap2 -t 2 -a ${GENOME}  ${SAMPLE}.fastq.gz -o $OUT_DIR/${BASENAME}.unsorted.sam
srun samtools view -@ 2 -q 10 -b $OUT_DIR/${BASENAME}.unsorted.sam | samtools sort -@ 2 - -o $OUT_DIR/${BASENAME}.sorted.bam 

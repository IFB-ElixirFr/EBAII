#!/bin/bash

#SBATCH --account=2538_eb3i_n1_2025
#SBATCH --job-name=bowtie_align
#SBATCH --cpus-per-task=8
#SBATCH --mem=20GB
#SBATCH --array=1-3

module load bowtie2/2.5.1
module load samtools/1.18
  
OUT_DIR="bowtie_res"
mkdir -p $OUT_DIR

GENOME_INDEX="/shared/bank/arabidopsis_thaliana/TAIR10.1/bowtie2/GCF_000001735.4_TAIR10.1_genomic"


SAMPLE=$(ls data/*R1* | sed 's/.fastq.gz//' | awk "NR==${SLURM_ARRAY_TASK_ID}")
BASENAME=$(basename ${SAMPLE/_R1/.bowtie_TAIR10.1})

  
srun bowtie2 -p 2 --mm -x ${GENOME_INDEX} -1 ${SAMPLE}.fastq.gz -2 ${SAMPLE/R1/R2}.fastq.gz -S bowtie_res/${BASENAME}.sam
srun samtools view -@ 2 -q 10 -b bowtie_res/${BASENAME}.sam | samtools sort -@ 2 - -o bowtie_res/${BASENAME}.bam 
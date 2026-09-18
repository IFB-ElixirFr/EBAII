#!/bin/bash

#SBATCH --account=2538_eb3i_n1_2025
#SBATCH --job-name=star_align
#SBATCH --cpus-per-task=8
#SBATCH --mem=20GB
#SBATCH --array=1-3

module load star/2.7.11a

OUT_DIR="/path/to/my/project/TP_rnaseq/workflow/star_res"
mkdir -p $OUT_DIR

DATA_DIR="/shared/projects/2422_ebaii_n1/atelier_rnaseq/04-Workflow/data/"
STAR_INDEX="/shared/bank/arabidopsis_thaliana/TAIR10.1/star-2.7.9a/"
GTF="/shared/bank/arabidopsis_thaliana/TAIR10.1/gtf/GCF_000001735.4_TAIR10.1_genomic.gtf"

R1IN=$(ls $DATA_DIR/*_R1.fastq.gz | awk "NR==${SLURM_ARRAY_TASK_ID}")
R2IN=${R1IN/_R1/_R2}
BASENAME=${R1IN/_R1.fastq.gz/.STAR_TAIR10.1_}

srun STAR --runThreadN ${SLURM_CPUS_PER_TASK} --genomeDir ${STAR_INDEX} \
  --sjdbGTFfile ${GTF} --readFilesCommand zcat --readFilesIn ${R1IN} ${R2IN} \
  --outFileNamePrefix ${OUT_DIR}/${BASENAME} --outSAMtype BAM SortedByCoordinate \
  --outSAMunmapped Within KeepPairs
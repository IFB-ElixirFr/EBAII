#!/bin/bash
#SBATCH --partition=fast
#SBATCH --job-name=my_fastqc
#SBATCH --account=2538_eb3i_n1_2025 # Modifier en fonction du projet
#SBATCH --cpus-per-task=1       # Modifier en fonction des besoins
#SBATCH --mem=4GB               # Idem
#SBATCH --array=1-3             # Modifier en fonction du nb de tâches à lancer en parallèle

mkdir -p fastqc
module load fastqc/0.11.9
mkdir -p trimmomatic
module load trimmomatic/0.39

# le Nième fichier de ma liste
SAMPLE=$(ls data/ | sed 's/R1.fastq.gz//' | \
head -n ${SLURM_ARRAY_TASK_ID} | tail -n 1)

srun --job-name FASTQC-$SAMPLE fastqc --outdir fastqc data/${SAMPLE}.fastq.gz
srun --job-name TRIM-$SAMPLE trimmomatic PE -threads 4 -phred33 \
                            data/${SAMPLE}.fastq.gz  data/${SAMPLE/R1/R2}.fastq.gz \
                            trimmomatic/${SAMPLE}.fastq.gz trimmomatic/${SAMPLE/R1/R2}.fastq.gz \
                            SLIDINGWINDOW:4:20 MINLEN:20
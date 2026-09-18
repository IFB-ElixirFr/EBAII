#!/bin/bash
#SBATCH --partition=fast
#SBATCH --job-name=my_fastqc
#SBATCH --account=2538_eb3i_n1_2025  # Modifier en fonction du projet
#SBATCH --cpus-per-task=1        # Modifier en fonction des besoins
#SBATCH --mem=4GB                # Idem
#SBATCH --array=1-6              # Modifier en fonction du nb de tâches à lancer en parallèle

mkdir -p fastqc
module load fastqc/0.11.9

# le Nième fichier de ma liste
SAMPLE=$(ls data/ | sed 's/.fastq.gz//' | \
head -n ${SLURM_ARRAY_TASK_ID} | tail -n 1)

srun --job-name $SAMPLE fastqc --outdir fastqc data/${SAMPLE}.fastq.gz
#!/bin/bash

#SBATCH --partition=fast
#SBATCH --job-name=my_fastqc
#SBATCH --account=2538_eb3i_n1_2025  # Modifier en fonction du projet
#SBATCH --cpus-per-task=1        # Modifier en fonction des besoins
#SBATCH --mem=4GB                # Idem
  
mkdir -p fastqc
module load fastqc/0.11.9

for SAMPLE in $(ls data/ | sed 's/.fastq.gz//')
do
	echo ">>> Processing $SAMPLE"
	srun --job-name $SAMPLE fastqc --outdir fastqc data/${SAMPLE}_R1.fastq.gz
	srun --job-name $SAMPLE fastqc --outdir fastqc data/${SAMPLE}_R2.fastq.gz
done
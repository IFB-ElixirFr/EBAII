#!/usr/bin/sh

## Charge tools you need for analysis
module load star/2.6
module load samtools/1.9
module load subread/1.6.1
module load multiqc/1.7

#indexation of genome for STAR !!!!!!! DON'T RUN THIS COMMAND !!!!
mkdir TAIR10
srun -c 4 STAR --runMode genomeGenerate --genomeFastaFiles Arabidopsis_thaliana.TAIR10.dna.toplevel.fa --genomeDir TAIR10 --runThreadN 4
srun samtools faidx Arabidopsis_thaliana.TAIR10.dna.toplevel.fa


# Copy fastq files in your home: Choose one WT and one KO
mkdir data	
cp -r /shared/projects/ebaii2020/atelier_rnaseq/01-Bioinfo/data/WT1* data/
cp -r /shared/projects/ebaii2020/atelier_rnaseq/01-Bioinfo/data/KO1* data/

# create a directory to store the future results
mkdir results/


#example of command with one sample:

srun -c 4 STAR --genomeDir /shared/projects/ebaii2020/atelier_rnaseq/01-Bioinfo/TAIR10 --runThreadN 4 --readFilesCommand zcat --outFileNamePrefix results/WT1_ --readFilesIn data/WT1_1.fastq.gz data/WT1_2.fastq.gz --outSAMtype BAM SortedByCoordinate --alignIntronMax 1000 --alignMatesGapMax 10000 --sjdbGTFfile /shared/projects/ebaii2020/atelier_rnaseq/01-Bioinfo/TAIR10/Arabidopsis_thaliana.TAIR10.45.gtf

srun samtools index results/WT1_Aligned.sortedByCoord.out.bam


# TODO: complete the command of featureCount. Google is your friend !!
srun -c 4 featureCounts -T 4 ...


#copy commands and change them to run all analysis (5min by sample)

#KO1
srun -c 4 STAR --genomeDir /shared/projects/ebaii2020/atelier_rnaseq/01-Bioinfo/TAIR10 --runThreadN 4 --readFilesCommand zcat --outFileNamePrefix results/KO1_ --readFilesIn data/KO1_1.fastq.gz data/KO1_2.fastq.gz --outSAMtype BAM SortedByCoordinate --alignIntronMax 1000 --alignMatesGapMax 10000 --sjdbGTFfile /shared/projects/ebaii2020/atelier_rnaseq/01-Bioinfo/TAIR10/Arabidopsis_thaliana.TAIR10.45.gtf 
srun samtools index results/KO1_Aligned.sortedByCoord.out.bam
srun -c 4 featureCounts -T 4 ...

#run mapping and counting in a few lines, doesn't matter the number of samples.
for sample in WT1 KO1 ...



#deconnect your session
logout

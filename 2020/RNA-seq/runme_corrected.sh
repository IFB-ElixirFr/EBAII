#!/usr/bin/sh

## Charge tools you need for analysis
module load star/2.6
module load samtools/1.9
module load subread/1.6.1
module load multiqc/1.7

#indexation of genome for STAR !!!!!!! DON'T RUN THIS COMMAND !!!!
mkdir TAIR10
srun -c 4 STAR --runMode genomeGenerate --genomeFastaFiles Arabidopsis_thaliana.TAIR10.dna.toplevel.fa --genomeDir TAIR10 --runThreadN 4
samtools faidx Arabidopsis_thaliana.TAIR10.dna.toplevel.fa


# Copy fastq files in your home: Choose one WT and one KO
mkdir data	
cp -r /shared/projects/ebai2019/atelier_rnaseq/data/WT* data/
cp -r /shared/projects/ebai2019/atelier_rnaseq/data/KO* data/

# create a directory to store the future results
mkdir results/


#example of command with one sample:

srun -c 4 STAR --genomeDir /shared/projects/ebai2019/atelier_rnaseq/TAIR10 --runThreadN 4 --readFilesCommand zcat --outFileNamePrefix results/WT1_ --readFilesIn data/WT1_1.fastq.gz data/WT1_2.fastq.gz --outSAMtype BAM SortedByCoordinate --alignIntronMax 1000 --alignMatesGapMax 10000 --sjdbGTFfile /shared/projects/ebai2019/atelier_rnaseq/TAIR10/Arabidopsis_thaliana.TAIR10.45.gtf

srun samtools index results/WT1_Aligned.sortedByCoord.out.bam


#copy commands and change them to run all analysis (5min by sample)
#WT1
srun -c 4 STAR --genomeDir /shared/projects/ebai2019/atelier_rnaseq/TAIR10 --runThreadN 4 --readFilesCommand zcat --outFileNamePrefix results/WT1_ --readFilesIn data/WT1_1.fastq.gz data/WT1_2.fastq.gz --outSAMtype BAM SortedByCoordinate --alignIntronMax 1000 --alignMatesGapMax 10000 --sjdbGTFfile /shared/projects/ebai2019/atelier_rnaseq/TAIR10/Arabidopsis_thaliana.TAIR10.45.gtf
srun samtools index results/WT1_Aligned.sortedByCoord.out.bam
srun -c 4 featureCounts -T 4 -t exon -g gene_id -s 1 -a /shared/projects/ebai2019/atelier_rnaseq/TAIR10/Arabidopsis_thaliana.TAIR10.45.gtf -o results/WT1_feature.out results/WT1_Aligned.sortedByCoord.out.bam 2> results/WT1_counting.logs  

#WT2
srun -c 4 STAR --genomeDir /shared/projects/ebai2019/atelier_rnaseq/TAIR10 --runThreadN 4 --readFilesCommand zcat --outFileNamePrefix results/WT2_ --readFilesIn data/WT2_1.fastq.gz data/WT2_2.fastq.gz --outSAMtype BAM SortedByCoordinate --alignIntronMax 1000 --alignMatesGapMax 10000 --sjdbGTFfile /shared/projects/ebai2019/atelier_rnaseq/TAIR10/Arabidopsis_thaliana.TAIR10.45.gtf
srun samtools index results/WT2_Aligned.sortedByCoord.out.bam
srun -c 4 featureCounts -T 4 -t exon -g gene_id -s 1 -a /shared/projects/ebai2019/atelier_rnaseq/TAIR10/Arabidopsis_thaliana.TAIR10.45.gtf -o results/WT2_feature.out results/WT2_Aligned.sortedByCoord.out.bam 2> results/WT2_counting.logs 

#WT3
srun -c 4 STAR --genomeDir /shared/projects/ebai2019/atelier_rnaseq/TAIR10 --runThreadN 4 --readFilesCommand zcat --outFileNamePrefix results/WT3_ --readFilesIn data/WT3_1.fastq.gz data/WT3_2.fastq.gz --outSAMtype BAM SortedByCoordinate --alignIntronMax 1000 --alignMatesGapMax 10000 --sjdbGTFfile /shared/projects/ebai2019/atelier_rnaseq/TAIR10/Arabidopsis_thaliana.TAIR10.45.gtf
srun samtools index results/WT3_Aligned.sortedByCoord.out.bam
srun -c 4 featureCounts -T 4 -t exon -g gene_id -s 1 -a /shared/projects/ebai2019/atelier_rnaseq/TAIR10/Arabidopsis_thaliana.TAIR10.45.gtf -o results/WT3_feature.out results/WT3_Aligned.sortedByCoord.out.bam 2> results/WT3_counting.logs 

#KO1
srun -c 4 STAR --genomeDir /shared/projects/ebai2019/atelier_rnaseq/TAIR10 --runThreadN 4 --readFilesCommand zcat --outFileNamePrefix results/KO1_ --readFilesIn data/KO1_1.fastq.gz data/KO1_2.fastq.gz --outSAMtype BAM SortedByCoordinate --alignIntronMax 1000 --alignMatesGapMax 10000 --sjdbGTFfile /shared/projects/ebai2019/atelier_rnaseq/TAIR10/Arabidopsis_thaliana.TAIR10.45.gtf
srun samtools index results/KO1_Aligned.sortedByCoord.out.bam
srun -c 4 featureCounts -T 4 -t exon -g gene_id -s 1 -a /shared/projects/ebai2019/atelier_rnaseq/TAIR10/Arabidopsis_thaliana.TAIR10.45.gtf -o results/KO1_feature.out results/KO1_Aligned.sortedByCoord.out.bam 2> results/KO1_counting.logs 

#KO2
srun -c 4 STAR --genomeDir /shared/projects/ebai2019/atelier_rnaseq/TAIR10 --runThreadN 4 --readFilesCommand zcat --outFileNamePrefix results/KO2_ --readFilesIn data/KO2_1.fastq.gz data/KO2_2.fastq.gz --outSAMtype BAM SortedByCoordinate --alignIntronMax 1000 --alignMatesGapMax 10000 --sjdbGTFfile /shared/projects/ebai2019/atelier_rnaseq/TAIR10/Arabidopsis_thaliana.TAIR10.45.gtf
srun samtools index results/KO2_Aligned.sortedByCoord.out.bam
srun -c 4 featureCounts -T 4 -t exon -g gene_id -s 1 -a /shared/projects/ebai2019/atelier_rnaseq/TAIR10/Arabidopsis_thaliana.TAIR10.45.gtf -o results/KO2_feature.out results/KO2_Aligned.sortedByCoord.out.bam 2> results/KO2_counting.logs 

#KO3
srun -c 4 STAR --genomeDir /shared/projects/ebai2019/atelier_rnaseq/TAIR10 --runThreadN 4 --readFilesCommand zcat --outFileNamePrefix results/KO3_ --readFilesIn data/KO3_1.fastq.gz data/KO3_2.fastq.gz --outSAMtype BAM SortedByCoordinate --alignIntronMax 1000 --alignMatesGapMax 10000 --sjdbGTFfile /shared/projects/ebai2019/atelier_rnaseq/TAIR10/Arabidopsis_thaliana.TAIR10.45.gtf
srun samtools index results/KO3_Aligned.sortedByCoord.out.bam
srun -c 4 featureCounts -T 4 -t exon -g gene_id -s 1 -a /shared/projects/ebai2019/atelier_rnaseq/TAIR10/Arabidopsis_thaliana.TAIR10.45.gtf -o results/KO3_feature.out results/KO3_Aligned.sortedByCoord.out.bam 2> results/KO3_counting.logs 


#run mapping and counting in a few lines, doesn't matter the number of samples.
for sample in WT1 WT2 WT3 KO1 KO2 KO3 ; do
	echo $sample	
	srun -c 4 STAR --genomeDir /shared/projects/ebai2019/atelier_rnaseq/TAIR10 --runThreadN 4 --readFilesCommand zcat --outFileNamePrefix results/${sample}_ --readFilesIn data/${sample}_1.fastq.gz data/${sample}_2.fastq.gz --outSAMtype BAM SortedByCoordinate --alignIntronMax 1000 --alignMatesGapMax 10000 --sjdbGTFfile /shared/projects/ebai2019/atelier_rnaseq/TAIR10/Arabidopsis_thaliana.TAIR10.45.gtf
	srun samtools index results/${sample}_Aligned.sortedByCoord.out.bam
	srun -c 4 featureCounts  -T 4 -t exon -g gene_id -s 1 -a /shared/projects/ebai2019/atelier_rnaseq/TAIR10/Arabidopsis_thaliana.TAIR10.45.gtf -o ${sample}_feature.out results/${sample}_Aligned.sortedByCoord.out.bam 2> results/${sample}_counting.logs
done


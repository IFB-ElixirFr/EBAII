#!/usr/bin/env bash
# Commands to run step by step, in the order shown in the course.
# Replace <your_project> with your actual project directory.
# Wait for all SLURM mapping jobs to finish before processing SAM files.

cd "/shared/projects/<your_project>"
mkdir EBA3I_2026_chipseq
cd EBA3I_2026_chipseq
cp -r /shared/projects/2621_eb3in1_2026/atelier_chipseq/training/data .

# FastQC
mkdir 01-QualityControl
cd 01-QualityControl
module load fastqc/0.11.9
fastqc --help
fastqc ../data/RtrC_ChIP_rep1_1.fastq.gz ../data/RtrC_ChIP_rep1_2.fastq.gz -o .
fastqc ../data/RtrC_ChIP_rep2_1.fastq.gz ../data/RtrC_ChIP_rep2_2.fastq.gz -o .
fastqc ../data/RtrC_input_rep1_1.fastq.gz ../data/RtrC_input_rep1_2.fastq.gz -o .
fastqc ../data/RtrC_input_rep2_1.fastq.gz ../data/RtrC_input_rep2_2.fastq.gz -o .
module unload fastqc/0.11.9
module load multiqc/1.35
multiqc . -o . -n multiqc_report.html
module unload multiqc/1.35
cd ..

# Bowtie2
module load bowtie2/2.5.1
bowtie2
mkdir 02-Mapping
cd 02-Mapping
mkdir index
cd index
bowtie2-build
bowtie2-build ../../data/Caulobacter_CP001340.fasta Caulobacter_CP001340
cd ..
mkdir bam
cd bam
sbatch -p fast -o RtrC_ChIP_rep1.mapping.out --cpus-per-task 10 --wrap="bowtie2 -p 10 --mm -3 1 -x ../index/Caulobacter_CP001340 -1 ../../data/RtrC_ChIP_rep1_1.fastq.gz -2 ../../data/RtrC_ChIP_rep1_2.fastq.gz -S RtrC_ChIP_rep1.sam"
sbatch -p fast -o RtrC_ChIP_rep2.mapping.out --cpus-per-task 10 --wrap="bowtie2 -p 10 --mm -3 1 -x ../index/Caulobacter_CP001340 -1 ../../data/RtrC_ChIP_rep2_1.fastq.gz -2 ../../data/RtrC_ChIP_rep2_2.fastq.gz -S RtrC_ChIP_rep2.sam"
sbatch -p fast -o RtrC_input_rep1.mapping.out --cpus-per-task 10 --wrap="bowtie2 -p 10 --mm -3 1 -x ../index/Caulobacter_CP001340 -1 ../../data/RtrC_input_rep1_1.fastq.gz -2 ../../data/RtrC_input_rep1_2.fastq.gz -S RtrC_input_rep1.sam"
sbatch -p fast -o RtrC_input_rep2.mapping.out --cpus-per-task 10 --wrap="bowtie2 -p 10 --mm -3 1 -x ../index/Caulobacter_CP001340 -1 ../../data/RtrC_input_rep2_1.fastq.gz -2 ../../data/RtrC_input_rep2_2.fastq.gz -S RtrC_input_rep2.sam"

# Wait for all four jobs to complete and check mapping logs before continuing.
module load samtools/1.21
samtools view -@ 2 -q 10 -b RtrC_ChIP_rep1.sam | samtools sort -@ 2 - -o RtrC_ChIP_rep1.bam
samtools index RtrC_ChIP_rep1.bam
gzip RtrC_ChIP_rep1.sam
samtools view -@ 2 -q 10 -b RtrC_ChIP_rep2.sam | samtools sort -@ 2 - -o RtrC_ChIP_rep2.bam
samtools index RtrC_ChIP_rep2.bam
gzip RtrC_ChIP_rep2.sam
samtools view -@ 2 -q 10 -b RtrC_input_rep1.sam | samtools sort -@ 2 - -o RtrC_input_rep1.bam
samtools index RtrC_input_rep1.bam
gzip RtrC_input_rep1.sam
samtools view -@ 2 -q 10 -b RtrC_input_rep2.sam | samtools sort -@ 2 - -o RtrC_input_rep2.bam
samtools index RtrC_input_rep2.bam
gzip RtrC_input_rep2.sam
module unload samtools/1.21 bowtie2/2.5.1
cd ..
module load multiqc/1.35
multiqc bam -o . -n multiqc_mapping.html
module unload multiqc/1.35
cd bam

# Picard duplicates
module load picard/2.23.5
module load samtools/1.21
picard MarkDuplicates \
-CREATE_INDEX true \
-INPUT RtrC_ChIP_rep1.bam \
-OUTPUT Marked_RtrC_ChIP_rep1.bam \
-METRICS_FILE RtrC_ChIP_rep1.dup_metrics.txt
# samtools flagstat Marked_RtrC_ChIP_rep1.bam > Marked_RtrC_ChIP_rep1.flagstat.txt
picard MarkDuplicates \
-CREATE_INDEX true \
-INPUT RtrC_ChIP_rep2.bam \
-OUTPUT Marked_RtrC_ChIP_rep2.bam \
-METRICS_FILE RtrC_ChIP_rep2.dup_metrics.txt
# samtools flagstat Marked_RtrC_ChIP_rep2.bam > Marked_RtrC_ChIP_rep2.flagstat.txt
picard MarkDuplicates \
-CREATE_INDEX true \
-INPUT RtrC_input_rep1.bam \
-OUTPUT Marked_RtrC_input_rep1.bam \
-METRICS_FILE RtrC_input_rep1.dup_metrics.txt
# samtools flagstat Marked_RtrC_input_rep1.bam > Marked_RtrC_input_rep1.flagstat.txt
picard MarkDuplicates \
-CREATE_INDEX true \
-INPUT RtrC_input_rep2.bam \
-OUTPUT Marked_RtrC_input_rep2.bam \
-METRICS_FILE RtrC_input_rep2.dup_metrics.txt
# samtools flagstat Marked_RtrC_input_rep2.bam > Marked_RtrC_input_rep2.flagstat.txt
module load multiqc/1.35
multiqc . -o . -n multiqc_duplicates.html
module unload multiqc/1.35
module unload samtools/1.21 picard/2.23.5
cd ../..

# ChIP fingerprint
mkdir 03-ChIPQualityControls
cd 03-ChIPQualityControls
module load deeptools/3.5.4
plotFingerprint \
  -p 2 \
  --numberOfSamples 10000 \
  -b ../02-Mapping/bam/Marked_RtrC_ChIP_rep1.bam \
     ../02-Mapping/bam/Marked_RtrC_ChIP_rep2.bam \
     ../02-Mapping/bam/Marked_RtrC_input_rep1.bam \
     ../02-Mapping/bam/Marked_RtrC_input_rep2.bam \
  -plot fingerprint_10000.png

plotFingerprint \
  -p 2 \
  -b ../02-Mapping/bam/Marked_RtrC_ChIP_rep1.bam \
     ../02-Mapping/bam/Marked_RtrC_ChIP_rep2.bam \
     ../02-Mapping/bam/Marked_RtrC_input_rep1.bam \
     ../02-Mapping/bam/Marked_RtrC_input_rep2.bam \
  -plot fingerprint.png

module unload deeptools/3.5.4
cd ..

# BigWig
module load deeptools/3.5.4
bamCoverage --help
mkdir 04-Visualization
cd 04-Visualization
bamCoverage \
  --bam ../02-Mapping/bam/Marked_RtrC_ChIP_rep1.bam \
  --outFileName RtrC_ChIP_rep1_nodup.bw \
  --outFileFormat bigwig \
  --normalizeUsing CPM \
  --skipNonCoveredRegions \
  --extendReads \
  --ignoreDuplicates
bamCoverage \
  --bam ../02-Mapping/bam/Marked_RtrC_ChIP_rep2.bam \
  --outFileName RtrC_ChIP_rep2_nodup.bw \
  --outFileFormat bigwig \
  --normalizeUsing CPM \
  --skipNonCoveredRegions \
  --extendReads \
  --ignoreDuplicates
bamCoverage \
  --bam ../02-Mapping/bam/Marked_RtrC_input_rep1.bam \
  --outFileName RtrC_input_rep1_nodup.bw \
  --outFileFormat bigwig \
  --normalizeUsing CPM \
  --skipNonCoveredRegions \
  --extendReads \
  --ignoreDuplicates
bamCoverage \
  --bam ../02-Mapping/bam/Marked_RtrC_input_rep2.bam \
  --outFileName RtrC_input_rep2_nodup.bw \
  --outFileFormat bigwig \
  --normalizeUsing CPM \
  --skipNonCoveredRegions \
  --extendReads \
  --ignoreDuplicates
cd ..
module unload deeptools/3.5.4

# MACS2
mkdir 05-PeakCalling
mkdir 05-PeakCalling/replicates
cd 05-PeakCalling/replicates
module load macs2/2.2.7.1
macs2 callpeak --help
macs2 callpeak \
  -t ../../02-Mapping/bam/RtrC_ChIP_rep1.bam \
  -c ../../02-Mapping/bam/RtrC_input_rep1.bam \
  --format BAMPE \
  --gsize 4042929 \
  --name 'RtrC_rep1' \
  -p 1e-2 \
  &> rep1_MACS.out
macs2 callpeak \
  -t ../../02-Mapping/bam/RtrC_ChIP_rep2.bam \
  -c ../../02-Mapping/bam/RtrC_input_rep2.bam \
  --format BAMPE \
  --gsize 4042929 \
  --name 'RtrC_rep2' \
  -p 1e-2 \
  &> rep2_MACS.out
cd ..
mkdir pool
cd pool
macs2 callpeak \
  -t ../../02-Mapping/bam/RtrC_ChIP_rep1.bam \
     ../../02-Mapping/bam/RtrC_ChIP_rep2.bam \
  -c ../../02-Mapping/bam/RtrC_input_rep1.bam \
     ../../02-Mapping/bam/RtrC_input_rep2.bam \
  --format BAMPE \
  --gsize 4042929 \
  --name 'RtrC_pool' \
  -p 1e-2 \
  &> pool_MACS.out
cd ..
module load multiqc/1.35
multiqc replicates pool -o . -n multiqc_peakcalling.html
module unload multiqc/1.35
mkdir idr
cd idr
module load idr/2.0.4.2
idr --help
idr --samples ../replicates/RtrC_rep1_peaks.narrowPeak \
              ../replicates/RtrC_rep2_peaks.narrowPeak \
    --peak-list ../pool/RtrC_pool_peaks.narrowPeak \
    --input-file-type narrowPeak \
    --output-file RtrC_idr_results.bed --plot

# center peaks and trim +- 250bp around to prepare input for motif analysis
awk 'BEGIN{OFS="\t"} $12 >= 1.30103 {print $1,$2+$10-250,$2+$10+250,"RtrC_idr_peak_"NR,$12}' RtrC_idr_results.bed > RtrC_idr_peaks.bed

module unload macs2/2.2.7.1
module unload idr/2.0.4.2
cd ../..

# Motif sequence extraction
mkdir 06-MotifAnalysis
cd 06-MotifAnalysis
module load samtools/1.21
samtools faidx ../data/Caulobacter_CP001340.fasta
module load bedtools/2.30.0
bedtools getfasta \
  -fi ../data/Caulobacter_CP001340.fasta \
  -bed ../05-PeakCalling/idr/RtrC_idr_peaks.bed \
  -fo RtrC_idr_peaks.fa
cd ..

# Peak annotation
mkdir 07-PeakAnnotation
cd 07-PeakAnnotation
cut \
  -f 1-5 \
  ../05-PeakCalling/idr/RtrC_idr_peaks.bed | \
  awk -F "\t" '{print $0"\t+"}' \
  > RtrC_idr_peaks.bed
module load homer/4.11
annotatePeaks.pl --help
annotatePeaks.pl \
  RtrC_idr_peaks.bed \
  ../data/Caulobacter_CP001340.fasta \
  -gtf ../data/Caulobacter_CP001340.gtf \
  > RtrC_idr_annotated_peaks.tsv
cd ..

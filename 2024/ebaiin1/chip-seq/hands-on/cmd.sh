USERID=slegras

## Création de l'environnement
cd /shared/projects/2538_eb3i_n1_2025/
mkdir $USERID
cd $USERID
cp -r /shared/projects/2538_eb3i_n1_2025/training_data/data .

## Controle qualité
module load fastqc/0.12.1

mkdir 01-QualityControl
cd 01-QualityControl
fastqc ../data/FNR_IP_ChIP-seq_Anaerobic_A.fastq.gz -o .
fastqc ../data/FNR_IP_ChIP-seq_Anaerobic_B.fastq.gz -o .
fastqc ../data/Anaerobic_INPUT_DNA.fastq.gz -o .
cd ..

## Mapping
module load bowtie2/2.5.1

mkdir 02-Mapping
cd 02-Mapping
mkdir index
cd index
bowtie2-build ../../data/Escherichia_coli_K12.fasta Escherichia_coli_K12
cd ..
mkdir bam
cd bam
sbatch -A 2538_eb3i_n1_2025 -p fast -o FNR_IP_ChIP-seq_Anaerobic_A.mapping.out --cpus-per-task 10 --wrap="bowtie2 -p 10 --mm -3 1 -x ../index/Escherichia_coli_K12 -U ../../data/FNR_IP_ChIP-seq_Anaerobic_A.fastq.gz -S FNR_IP_ChIP-seq_Anaerobic_A.sam"
sbatch -A 2538_eb3i_n1_2025 -p fast -o FNR_IP_ChIP-seq_Anaerobic_B.mapping.out --cpus-per-task 10 --wrap="bowtie2 -p 10 --mm -3 1 -x ../index/Escherichia_coli_K12 -U ../../data/FNR_IP_ChIP-seq_Anaerobic_B.fastq.gz -S FNR_IP_ChIP-seq_Anaerobic_B.sam"
sbatch -A 2538_eb3i_n1_2025 -p fast -o Anaerobic_INPUT_DNA.mapping.out --cpus-per-task 10 --wrap="bowtie2 -p 10 --mm -3 1 -x ../index/Escherichia_coli_K12 -U ../../data/Anaerobic_INPUT_DNA.fastq.gz -S Anaerobic_INPUT_DNA.sam"

# sbatch --cpus-per-task 10 --wrap="bowtie -p 10 ../index/Escherichia_coli_K12 ../../data/Anaerobic_INPUT_DNA.fastq.gz -v 2 -m 1 -3 1 -S 2> Anaerobic_INPUT_DNA.out > Anaerobic_INPUT_DNA.sam"
# sbatch --cpus-per-task 10 --wrap="bowtie -p 10 ../index/Escherichia_coli_K12 ../../data/FNR_IP_ChIP-seq_Anaerobic_B.fastq.gz -v 2 -m 1 -3 1 -S 2> FNR_IP_ChIP-seq_Anaerobic_B.out > FNR_IP_ChIP-seq_Anaerobic_B.sam"
# sbatch --cpus-per-task 10 --wrap="bowtie -p 10 ../index/Escherichia_coli_K12 ../../data/FNR_IP_ChIP-seq_Anaerobic_A.fastq.gz -v 2 -m 1 -3 1 -S 2> FNR_IP_ChIP-seq_Anaerobic_A.out > FNR_IP_ChIP-seq_Anaerobic_A.sam"

# Création de fichiers sam ordonnés
module load samtools/1.21
samtools view -@ 2 -q 10 -b FNR_IP_ChIP-seq_Anaerobic_A.sam | samtools sort -@ 2 - -o FNR_IP_ChIP-seq_Anaerobic_A.bam 
samtools view -@ 2 -q 10 -b FNR_IP_ChIP-seq_Anaerobic_B.sam | samtools sort -@ 2 - -o FNR_IP_ChIP-seq_Anaerobic_B.bam 
samtools view -@ 2 -q 10 -b Anaerobic_INPUT_DNA.sam | samtools sort -@ 2 - -o Anaerobic_INPUT_DNA.bam 

# On index les fichiers bam
samtools index FNR_IP_ChIP-seq_Anaerobic_B.bam
samtools index FNR_IP_ChIP-seq_Anaerobic_A.bam
samtools index Anaerobic_INPUT_DNA.bam

# On compresse les fichiers sam
gzip FNR_IP_ChIP-seq_Anaerobic_A.sam &
gzip FNR_IP_ChIP-seq_Anaerobic_B.sam &
gzip Anaerobic_INPUT_DNA.sam &

# marquage des duplicats
module load picard/2.23.5

picard MarkDuplicates \
  -CREATE_INDEX true \
  -INPUT FNR_IP_ChIP-seq_Anaerobic_A.bam \
  -OUTPUT Marked_FNR_IP_ChIP-seq_Anaerobic_A.bam \
  -METRICS_FILE Marked_FNR_IP_ChIP-seq_Anaerobic_A.metric &

picard MarkDuplicates \
  -CREATE_INDEX true \
  -INPUT FNR_IP_ChIP-seq_Anaerobic_B.bam \
  -OUTPUT Marked_FNR_IP_ChIP-seq_Anaerobic_B.bam \
  -METRICS_FILE Marked_FNR_IP_ChIP-seq_Anaerobic_B.metric &

picard MarkDuplicates \
  -CREATE_INDEX true \
  -INPUT Anaerobic_INPUT_DNA.bam \
  -OUTPUT Marked_Anaerobic_INPUT_DNA.bam \
  -METRICS_FILE Anaerobic_INPUT_DNA.metric &

cd ../..

## Control qualité
module load deeptools/3.5.4

mkdir 03-ChIPQualityControls
cd 03-ChIPQualityControls

plotFingerprint \
  -p 2 \
  --numberOfSamples 10000 \
  -b ../02-Mapping/bam/FNR_IP_ChIP-seq_Anaerobic_A.bam \
     ../02-Mapping/bam/FNR_IP_ChIP-seq_Anaerobic_B.bam \
     ../02-Mapping/bam/Anaerobic_INPUT_DNA.bam \
  -plot fingerprint_10000.png &

# plotFingerprint \
#   -p 2 \
#   -b ../02-Mapping/bam/FNR_IP_ChIP-seq_Anaerobic_A.bam \
#   ../02-Mapping/bam/FNR_IP_ChIP-seq_Anaerobic_B.bam \
#   ../02-Mapping/bam/Anaerobic_INPUT_DNA.bam \
#   -plot fingerprint.png &

cd ..

## Génération de fichiers bigwig
module load deeptools/3.5.4

mkdir 04-Visualization
cd 04-Visualization/
bamCoverage \
  --bam ../02-Mapping/bam/Marked_FNR_IP_ChIP-seq_Anaerobic_A.bam \
  --outFileName FNR_IP_ChIP-seq_Anaerobic_A_nodup.bw \
  --outFileFormat bigwig \
  --effectiveGenomeSize 4639675 \
  --normalizeUsing CPM \
  --skipNonCoveredRegions \
  --extendReads 200 \
  --ignoreDuplicates &

bamCoverage \
  --bam ../02-Mapping/bam/Marked_FNR_IP_ChIP-seq_Anaerobic_B.bam \
  --outFileName FNR_IP_ChIP-seq_Anaerobic_B_nodup.bw \
  --outFileFormat bigwig \
  --effectiveGenomeSize 4639675 \
  --normalizeUsing CPM \
  --skipNonCoveredRegions \
  --extendReads 200 \
  --ignoreDuplicates &

bamCoverage \
  --bam ../02-Mapping/bam/Marked_Anaerobic_INPUT_DNA.bam \
  --outFileName Anaerobic_INPUT_DNA_nodup.bw \
  --outFileFormat bigwig \
  --effectiveGenomeSize 4639675 \
  --normalizeUsing CPM \
  --skipNonCoveredRegions \
  --extendReads 200 \
  --ignoreDuplicates &

cd ..

## Peak calling
module load macs2/2.2.7.1

mkdir 05-PeakCalling
# Peak calling sur les réplicats
mkdir 05-PeakCalling/replicates
cd 05-PeakCalling/replicates
macs2 callpeak \
  -t ../../02-Mapping/bam/FNR_IP_ChIP-seq_Anaerobic_A.bam \
  -c ../../02-Mapping/bam/Anaerobic_INPUT_DNA.bam \
  --format BAM \
  --gsize 4639675 \
  --name 'FNR_Anaerobic_A' \
  --bw 400 \
  --fix-bimodal \
  -p 1e-2 \
  &> repA_MACS.out &

macs2 callpeak \
  -t ../../02-Mapping/bam/FNR_IP_ChIP-seq_Anaerobic_B.bam \
  -c ../../02-Mapping/bam/Anaerobic_INPUT_DNA.bam \
  --format BAM \
  --gsize 4639675 \
  --name 'FNR_Anaerobic_B' \
  --bw 400 \
  --fix-bimodal \
  -p 1e-2 \
  &> repB_MACS.out &
cd ..

# Peak calling sur le pool de réplicat
mkdir pool
cd pool
macs2 callpeak \
  -t ../../02-Mapping/bam/FNR_IP_ChIP-seq_Anaerobic_A.bam \
  ../../02-Mapping/bam/FNR_IP_ChIP-seq_Anaerobic_B.bam \
  -c ../../02-Mapping/bam/Anaerobic_INPUT_DNA.bam \
  --format BAM \
  --gsize 4639675 \
  --name 'FNR_Anaerobic_pool' \
  --bw 400 \
  --fix-bimodal \
  -p 1e-2 \
  &> pool_MACS.out &
cd ..

# Analyse IDR
module load idr/2.0.4.2
mkdir idr
cd idr
idr \
  --samples ../replicates/FNR_Anaerobic_A_peaks.narrowPeak \
  ../replicates/FNR_Anaerobic_B_peaks.narrowPeak \
  --peak-list ../pool/FNR_Anaerobic_pool_peaks.narrowPeak \
  --input-file-type narrowPeak \
  --output-file FNR_anaerobic_idr_peaks.bed \
  --plot &
cd ../..

## Préparation des fichiers pour l'analyse de motif
module load bedtools/2.30.0

mkdir 06-MotifAnalysis
cd 06-MotifAnalysis
samtools faidx ../data/Escherichia_coli_K12.fasta
bedtools getfasta -fi ../data/Escherichia_coli_K12.fasta \
  -bed ../05-PeakCalling/idr/FNR_anaerobic_idr_peaks.bed -fo FNR_anaerobic_idr_peaks.fa

cd ..

## Annotation des pics
module load homer/4.11

mkdir 07-PeakAnnotation
cd 07-PeakAnnotation

# On met le fichier de pics dans le bon format
cut -f1-5 ../05-PeakCalling/idr/FNR_anaerobic_idr_peaks.bed | \
  awk -F "\t" '{print $0"\t+"}'  > FNR_anaerobic_idr_peaks.bed

# lancement de l'annotation
annotatePeaks.pl \
  FNR_anaerobic_idr_peaks.bed \
  ../data/Escherichia_coli_K12.fasta \
  -gtf ../data/Escherichia_coli_K_12_MG1655.annotation.fixed.gtf \
  > FNR_anaerobic_idr_annotated_peaks.tsv

# # calcul de statistiques sur les résultats
tail -n 2 FNR_anaerobic_idr_final_peaks_annotation.tsv | awk '{print $11}'
vi FNR_anaerobic_idr_final_peaks_annotation.tsv
tail -n +2 FNR_anaerobic_idr_final_peaks_annotation.tsv | awk '{print $11}'
tail -n +2 FNR_anaerobic_idr_final_peaks_annotation.tsv | awk '{print $8}' | sort | uniq -c
tail -n +2 FNR_anaerobic_idr_final_peaks_annotation.tsv | awk '{if ($8=="promoter-TSS") print $11}'
tail -n +2 FNR_anaerobic_idr_final_peaks_annotation.tsv | awk '{if ($8=="promoter-TSS") print $11}' | wc -l
tail -n +2 FNR_anaerobic_idr_final_peaks_annotation.tsv | awk '{if ($8=="promoter-TSS") print $11}' > FNR_anaerobic_idr_final_peaks_annotation_officialGeneSymbols.tsv

# # On compresse de nouveau les fichiers
# gzip ../data/Escherichia_coli_K_12_MG1655.annotation.fixed.gtf
# gzip ../data/Escherichia_coli_K12.fasta

cd ..

## Création de l'environnement
mkdir M2.2-BIMS-epigenomique
cd M2.2-BIMS-epigenomique
cp -r ../EBAII2021_chipseq/data .

## Controle qualité
module add fastqc/0.11.9

mkdir 01-QualityControl
cd 01-QualityControl
fastqc ../data/SRR576933.fastq.gz -o .
fastqc ../data/SRR576934.fastq.gz -o .
fastqc ../data/SRR576938.fastq.gz -o .
cd ..

## Mapping
module add bowtie/1.2.3

mkdir 02-Mapping
cd 02-Mapping
mkdir index
cd index
bowtie-build ../../data/Escherichia_coli_K12.fasta Escherichia_coli_K12
cd ..
mkdir bam
cd bam
sbatch --cpus-per-task 10 --wrap="bowtie -p 10 ../index/Escherichia_coli_K12 ../../data/SRR576938.fastq.gz -v 2 -m 1 -3 1 -S 2> SRR576938.out > SRR576938.sam"
sbatch --cpus-per-task 10 --wrap="bowtie -p 10 ../index/Escherichia_coli_K12 ../../data/SRR576934.fastq.gz -v 2 -m 1 -3 1 -S 2> SRR576934.out > SRR576934.sam"
sbatch --cpus-per-task 10 --wrap="bowtie -p 10 ../index/Escherichia_coli_K12 ../../data/SRR576933.fastq.gz -v 2 -m 1 -3 1 -S 2> SRR576933.out > SRR576933.sam"

# Création de fichiers sam ordonnés
module add samtools/1.10
samtools sort SRR576933.sam | samtools view -b > SRR576933.bam
samtools sort SRR576934.sam | samtools view -b > SRR576934.bam
samtools sort SRR576938.sam | samtools view -b > SRR576938.bam

# On index les fichiers bam
samtools index SRR576934.bam
samtools index SRR576933.bam
samtools index SRR576938.bam

# On compresse les fichiers sam
gzip SRR576933.sam &
gzip SRR576934.sam &
gzip SRR576938.sam &

# marquage des duplicats
module add picard/2.22.0

picard MarkDuplicates CREATE_INDEX=true INPUT=SRR576933.bam OUTPUT=Marked_SRR576933.bam METRICS_FILE=metrics VALIDATION_STRINGENCY=STRICT
picard MarkDuplicates CREATE_INDEX=true INPUT=SRR576934.bam OUTPUT=Marked_SRR576934.bam METRICS_FILE=metrics VALIDATION_STRINGENCY=STRICT
picard MarkDuplicates CREATE_INDEX=true INPUT=SRR576938.bam OUTPUT=Marked_SRR576938.bam METRICS_FILE=metrics VALIDATION_STRINGENCY=STRICT

cd ../..

## Control qualité
module add deeptools/3.2.0

mkdir 03-ChIPQualityControls
cd 03-ChIPQualityControls
plotFingerprint --numberOfSamples 10000 -b ../02-Mapping/bam/SRR576933.bam ../02-Mapping/bam/SRR576934.bam ../02-Mapping/bam/SRR576938.bam -plot fingerprint_10000.png &
plotFingerprint -b ../02-Mapping/bam/SRR576933.bam ../02-Mapping/bam/SRR576934.bam ../02-Mapping/bam/SRR576938.bam -plot fingerprint.png &

cd ..

## Génération de fichiers bigwig
module add deeptools/3.2.0

mkdir 04-Visualization
cd 04-Visualization/
bamCoverage --bam ../02-Mapping/bam/Marked_SRR576933.bam --outFileName SRR576933_nodup.bw --outFileFormat bigwig --effectiveGenomeSize 4639675 --normalizeUsing RPGC --skipNonCoveredRegions --extendReads 200 --ignoreDuplicates
bamCoverage --bam ../02-Mapping/bam/Marked_SRR576934.bam --outFileName SRR576934_nodup.bw --outFileFormat bigwig --effectiveGenomeSize 4639675 --normalizeUsing RPGC --skipNonCoveredRegions --extendReads 200 --ignoreDuplicates
bamCoverage --bam ../02-Mapping/bam/Marked_SRR576938.bam --outFileName SRR576938_nodup.bw --outFileFormat bigwig --effectiveGenomeSize 4639675 --normalizeUsing RPGC --skipNonCoveredRegions --extendReads 200 --ignoreDuplicates

cd ..

## Peak calling
module add macs2/2.2.7.1

mkdir 05-PeakCalling
# Peak calling sur les réplicats
mkdir 05-PeakCalling/replicates
cd 05-PeakCalling/replicates
macs2 callpeak -t ../../02-Mapping/bam/SRR576933.bam -c ../../02-Mapping/bam/SRR576938.bam --format BAM --gsize 4639675 --name 'FNR_Anaerobic_A' --bw 400 --fix-bimodal -p 1e-2 &> repA_MACS.out
macs2 callpeak -t ../../02-Mapping/bam/SRR576934.bam -c ../../02-Mapping/bam/SRR576938.bam --format BAM --gsize 4639675 --name 'FNR_Anaerobic_B' --bw 400 --fix-bimodal -p 1e-2 &> repB_MACS.out
cd ..

# Peak calling sur le pool de réplicat
mkdir pool
cd pool
macs2 callpeak -t ../../02-Mapping/bam/SRR576933.bam ../../02-Mapping/bam/SRR576934.bam -c ../../02-Mapping/bam/SRR576938.bam --format BAM --gsize 4639675 --name 'FNR_Anaerobic_pool' --bw 400 --fix-bimodal -p 1e-2 &> pool_MACS.out
cd ..

# Analyse IDR
mkdir idr
cd idr
idr --samples ../replicates/FNR_Anaerobic_A_peaks.narrowPeak ../replicates/FNR_Anaerobic_B_peaks.narrowPeak --peak-list ../pool/FNR_Anaerobic_pool_peaks.narrowPeak \
--input-file-type narrowPeak --output-file FNR_anaerobic_idr_peaks.bed --plot
cd ../..

## Préparation des fichiers pour l'analyse de motif
module add bedtools/2.29.2

mkdir 06-MotifAnalysis
cd 06-MotifAnalysis
samtools faidx ../data/Escherichia_coli_K12.fasta
bedtools getfasta -fi ../data/Escherichia_coli_K12.fasta \
-bed ../05-PeakCalling/idr/FNR_anaerobic_idr_peaks.bed -fo FNR_anaerobic_idr_peaks.fa

cd ..

## Annotation des pics
module add homer/4.10

mkdir 07-PeakAnnotation
cd 07-PeakAnnotation
gunzip ../data/Escherichia_coli_K_12_MG1655.annotation.fixed.gtf.gz
gunzip ../data/Escherichia_coli_K12.fasta.gz

# On met le fichier de pics dans le bon format
cut -f1-5 ../05-PeakCalling/idr/FNR_anaerobic_idr_peaks.bed |   awk -F "\t" '{print $0"\t+"}'  > FNR_anaerobic_idr_peaks.bed

# lancement de l'annotation
annotatePeaks.pl   FNR_anaerobic_idr_peaks.bed   ../data/Escherichia_coli_K12.fasta   -gtf ../data/Escherichia_coli_K_12_MG1655.annotation.fixed.gtf   > FNR_anaerobic_idr_annotated_peaks.tsv

# calcul de statistiques sur les résultats
tail -n 2 FNR_anaerobic_idr_final_peaks_annotation.tsv | awk '{print $11}'
vi FNR_anaerobic_idr_final_peaks_annotation.tsv
tail -n +2 FNR_anaerobic_idr_final_peaks_annotation.tsv | awk '{print $11}'
tail -n +2 FNR_anaerobic_idr_final_peaks_annotation.tsv | awk '{print $8}' | sort | uniq -c
tail -n +2 FNR_anaerobic_idr_final_peaks_annotation.tsv | awk '{if ($8=="promoter-TSS") print $11}'
tail -n +2 FNR_anaerobic_idr_final_peaks_annotation.tsv | awk '{if ($8=="promoter-TSS") print $11}' | wc -l
tail -n +2 FNR_anaerobic_idr_final_peaks_annotation.tsv | awk '{if ($8=="promoter-TSS") print $11}' > FNR_anaerobic_idr_final_peaks_annotation_officialGeneSymbols.tsv

# On compresse de nouveau les fichiers
gzip ../data/Escherichia_coli_K_12_MG1655.annotation.fixed.gtf
gzip ../data/Escherichia_coli_K12.fasta

cd ..

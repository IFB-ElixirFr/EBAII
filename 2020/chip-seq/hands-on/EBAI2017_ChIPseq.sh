#! /bin/env bash

###################################################
################# Variable definition
login=slegras ## To be changed with your login!!
home=/shared/projects/training/${login}/EBA2019_chipseq

###################################################
################# Working environment

## Create a working directory for entire hands-on part
mkdir $home
cd $home

## Create a directory for raw data
# mkdir data
cd data
## copy fastq/fasta files from local computer
## fastq file got downloaded from EBI
# ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR576/SRR576933/SRR576933.fastq.gz
# ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR576/SRR576934/SRR576934.fastq.gz
# ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR576/SRR576938/SRR576938.fastq.gz
## genome fasta file got downloaded from NCBI website
## annotation track (.gtf) was downloaded from UCSC table browser (gtf file)
## annotation track (tsv) was downloaded from UCSC table browser (selecting fields to be output)
# http://microbes.ucsc.edu/cgi-bin/hgTables?org=Escherichia+coli+K12&db=eschColi_K12&hgsid=1465191&hgta_doMainPage=1

## Changing chromosome names
# srun zcat Escherichia_coli_K_12_MG1655.annotation.gtf | perl -pe 's/^chr/gi\|49175990\|ref\|NC_000913.2\|/' | \
#  gzip > Escherichia_coli_K_12_MG1655.annotation.fixed.gtf.gz

## Copy data
srun cp -r /shared/projects/training/slegras/EBA2019_chipseq/data \
/shared/projects/training/slegras/EBA2019_chipseq/scripts/ .
cd $home

# ## Create a directory for scripts
# mkdir scripts
#
# ## Go to the newly created directory
# cd scripts
#
# ## Clone phantompeakqualtools repository
# git clone https://github.com/crazyhottommy/phantompeakqualtools.git

## go back to home working directory
# cd $home

## Loading conda ChIP-Seq environment
source activate eba2017_chipseq

###################################################
################# Quality controls
## Create a directory for quality control
mkdir 01-QualityControl

## Go to quality control directory
cd 01-QualityControl

## Test fastqc tool
## To be run with srun
srun fastqc --help

## Run fastqc on the IP
srun fastqc ../data/SRR576933.fastq.gz -o .

## Run fastqc on the control
srun fastqc ../data/SRR576938.fastq.gz -o .

## List output file
ls

## Go to home working directory
cd $home

###################################################
################# Mapping
## creating output directory for alignment
mkdir 02-Mapping

## Go to newly created directory
cd 02-Mapping

## Create a directory for index files
mkdir index

## Go to index directory
cd index

## testing bowtie-build
srun bowtie-build

## Unzip genome fasta file
srun gunzip ../../data/Escherichia_coli_K12.fasta.gz

## Creating genome index
srun bowtie-build ../../data/Escherichia_coli_K12.fasta Escherichia_coli_K12

## Compress back genome fasta file
srun gzip ../../data/Escherichia_coli_K12.fasta

## Go back to upper directory i.e 02-Mapping
cd ..

## Create a directory for IP alignment
mkdir IP

## Go to IP directory
cd IP

## Unzip fastq IP file
srun gunzip ../../data/SRR576933.fastq.gz

## Run alignment
srun bowtie ../index/Escherichia_coli_K12 ../../data/SRR576933.fastq -v 2 -m 1 -3 1 -S 2> SRR576933.out > SRR576933.sam

## Compress back fastq IP file
srun gzip ../../data/SRR576933.fastq

## Create a sorted bam file
srun samtools sort SRR576933.sam | samtools view -b > SRR576933.bam

## create an index for the bam file
srun samtools index SRR576933.bam

## Compress the .sam file
gzip SRR576933.sam

## Go back to upper directory
cd ..

## Create a directory for IP alignment
mkdir Control

## Go to Control directory
cd Control

## Unzip fastq IP file
srun gunzip ../../data/SRR576938.fastq.gz

## Run alignment
srun bowtie ../index/Escherichia_coli_K12 ../../data/SRR576938.fastq -v 2 -m 1 -3 1 -S 2> SRR576938.out > SRR576938.sam

## Compress back fastq IP file
srun gzip ../../data/SRR576938.fastq

## Create a sorted bam file
srun samtools sort SRR576938.sam | samtools view -b > SRR576938.bam

## create an index for the bam file
srun samtools index SRR576938.bam

## Compress the .sam file
gzip SRR576938.sam

## Go to the IP directory
cd ../IP

## Run picard markDuplicates on the IP sample
srun picard MarkDuplicates \
CREATE_INDEX=true \
INPUT=SRR576933.bam \
OUTPUT=Marked_SRR576933.bam \
METRICS_FILE=metrics \
VALIDATION_STRINGENCY=STRICT

## Go to the Control directory
cd ../Control

## Run picard markDuplicates on the Control sample
srun picard MarkDuplicates \
CREATE_INDEX=true \
INPUT=SRR576938.bam \
OUTPUT=Marked_SRR576938.bam \
METRICS_FILE=metrics \
VALIDATION_STRINGENCY=STRICT

## Go to home working directory
cd $home

###################################################
################# ChIP Quality Controls
## creating output directory for alignment
mkdir 03-ChIPQualityControls

## Go to newly created directory
cd 03-ChIPQualityControls

### Deeptools fingerprint
## Run deeptools fingerprint
srun plotFingerprint -b ../02-Mapping/IP/SRR576933.bam ../02-Mapping/Control/SRR576938.bam -plot fingerprint.png

### phantompeakqualtools
## convert the BAM file into TagAlign format, specific to the program that calculates the quality metrics
srun samtools view -F 0x0204 -o - ../02-Mapping/IP/SRR576933.bam | \
gawk 'BEGIN{OFS="\t"}{if (and($2,16) > 0) {print $3,($4-1),($4-1+length($10)),"N","1000","-"}
else {print $3,($4-1),($4-1+length($10)),"N","1000","+"} }' \
 | gzip -c > SRR576933_experiment.tagAlign.gz

 ## create an R environment and load it
 # conda create -c r --name eba2017_spp r
 # source activate eba2017_spp
 # conda install -c bioconda r-spp
 # conda install -c bioconda samtools
 # conda install -c bioconda gawk

 source activate eba2017_spp

## Run phantompeakqualtools
srun Rscript ../scripts/phantompeakqualtools/run_spp.R -c=SRR576933_experiment.tagAlign.gz  -savp -out=SRR576933_IP_phantompeaks

## Source back the chipseq environment
source activate eba2017_chipseq

## Go to home working directory
cd $home

###################################################
################# visualization
## Create a directory to store visualization files
mkdir 04-Visualization

## Go to the newly created directory
cd 04-Visualization

## Test bamCoverage
srun bamCoverage --help

## Run bamCoverage on IP
srun --mem=3G bamCoverage --bam ../02-Mapping/IP/Marked_SRR576933.bam \
--outFileName SRR576933_nodup.bw --outFileFormat bigwig --normalizeTo1x 4639675 \
--skipNonCoveredRegions --extendReads 200 --ignoreDuplicates

## Run bamCoverage on Control
srun --mem=5G bamCoverage --bam ../02-Mapping/Control/Marked_SRR576938.bam \
--outFileName SRR576938_nodup.bw --outFileFormat bigwig --normalizeTo1x 4639675 \
--skipNonCoveredRegions --extendReads 200 --ignoreDuplicates

###################################################
################# Peak Calling
## Create a directory to store peak calling result files
mkdir 05-PeakCalling

## Go to the newly created directory
cd 05-PeakCalling

## Check macs parameters
srun macs

## Run macs on the IP and the Control file
srun macs -t ../02-Mapping/IP/SRR576933.bam -c ../02-Mapping/Control/SRR576938.bam --format BAM  --gsize 4639675 \
--name "FNR_Anaerobic_A" --bw 400 --diag &> MACS.out

###################################################
################# Peak Annotation
## Create a directory to store peak annotation data
mkdir 06-PeakAnnotation

## Go to the newly created directory
cd 06-PeakAnnotation

### CEAS
srun ceas  --help

### Homer annotatePeaks
## Uncompress annotation file
srun gunzip ../data/Escherichia_coli_K_12_MG1655.annotation.fixed.gtf.gz

## Uncompress genome file
srun gunzip ../data/Escherichia_coli_K12.fasta.gz

## Create a BED file with 6 columns
awk -F "\t" '{print $0"\t+"}' ../05-PeakCalling/FNR_Anaerobic_A_peaks.bed > FNR_Anaerobic_A_peaks.bed

## Try it out
srun annotatePeaks.pl

## Annotation peaks with Homer
srun annotatePeaks.pl \
FNR_Anaerobic_A_peaks.bed \
../data/Escherichia_coli_K12.fasta \
-gtf ../data/Escherichia_coli_K_12_MG1655.annotation.fixed.gtf \
> FNR_Anaerobic_A_annotated_peaks.tsv

## Compress annotation file
srun gzip ../data/Escherichia_coli_K_12_MG1655.annotation.fixed.gtf

## Run srun in an interactive mode
srun --pty bash

## Load the environment with R
source activate eba2017_spp

## Launch R
R

## read the file with peaks annotated with homer
## data are loaded into a data frame
## sep="\t": this is a tab separated file
## h=T: there is a line with headers (ie. column names)
d <- read.table("FNR_Anaerobic_A_annotated_peaks.tsv", sep="\t", h=T)

## Load a 2-columns files which contains in the first column gene IDs
## and in the second column gene symbols
## data are loaded into a data frame
## h=F: there is no header line
gene.symbol <- read.table("../data/Escherichia_coli_K_12_MG1655.annotation.tsv.gz", h=F)

## Merge the 2 data frames based on a common field
## by.x gives the columns name in which the common field is for the d data frame
## by.y gives the columns name in which the common field is for the gene.symbol data frame
## d contains several columns with no information. We select only interesting columns
## -> d[,c(seq(1,6,1),8,10,11)]
d.annot <- merge(d[,c(seq(1,6,1),8,10,11)], gene.symbol, by.x="Nearest.PromoterID", by.y="V1") # to link the two tables by the gene ref

## Change column names of the resulting data frame
colnames(d.annot)[2] <- "PeakID"  # name the 2d column of the new file "PeakID"
colnames(d.annot)[dim(d.annot)[2]] <- "Gene.Symbol" # name the last column of the new file "Gene.Symbol"

## output the merged data frame to a file named "FNR_Anaerobic_A_final_peaks_annotation.tsv"
## col.names=T: output column names
## row.names=F: don't output row names
## sep="\t": table fields are separated by tabs
## quote=F: don't put quote around text.
write.table(d.annot, "FNR_Anaerobic_A_final_peaks_annotation.tsv", col.names=T, row.names=F, sep="\t", quote=F)

## Leave R
quit()

## Do not save the environment
n

## exit node and go back to master node
exit

## Retrieve the list of closest genes
tail -n +2 FNR_Anaerobic_A_final_peaks_annotation.tsv | awk '{print $11}'

## Retrieve only the genes that encode for proteins
tail -n +2 FNR_Anaerobic_A_final_peaks_annotation.tsv | awk '{print $8}' | sort | uniq -c

## How many protein-coding genes are there in the file?
tail -n +2 FNR_Anaerobic_A_final_peaks_annotation.tsv | awk '{if ($8=="promoter-TSS") print $11}'

## Is the number of genes in your file consistent with the previous reply?
tail -n +2 FNR_Anaerobic_A_final_peaks_annotation.tsv | awk '{if ($8=="promoter-TSS") print $11}' | wc -l
tail -n +2 FNR_Anaerobic_A_final_peaks_annotation.tsv | awk '{if ($8=="promoter-TSS") print $11}' \
> FNR_Anaerobic_A_final_peaks_annotation_officialGeneSymbols.tsv

cd $home

###################################################
################# Motif analysis
## Create a directory to store data needed from motif analysis
mkdir 07-MotifAnalysis

## Go to the newly created directory
cd 07-MotifAnalysis

### Extract fasta sequence of peaks

## Create an index of fasta file
srun samtools faidx ../data/Escherichia_coli_K12.fasta

## Extract fasta sequence from genomic coordinate of peaks
srun bedtools getfasta -fi ../data/Escherichia_coli_K12.fasta \
-bed ../05-PeakCalling/FNR_Anaerobic_A_peaks.bed -fo FNR_Anaerobic_A_peaks.fa

### Extract fasta sequence of peaks
## Extract genomic coordinates of peaks summit +/- 100bp
srun bedtools slop -b 100 -i ../05-PeakCalling/FNR_Anaerobic_A_summits.bed \
-g ../data/Escherichia_coli_K12.fasta.fai > FNR_Anaerobic_A_summits+-100.bed

## Extract fasta sequence from genomic coordinate of peaks
srun bedtools getfasta -fi ../data/Escherichia_coli_K12.fasta -bed FNR_Anaerobic_A_summits+-100.bed -fo FNR_Anaerobic_A_summits+-100.fa

## Compress back the genome file
srun gzip ../data/Escherichia_coli_K12.fasta

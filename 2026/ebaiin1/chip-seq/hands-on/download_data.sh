#!/usr/bin/env bash

declare -A geo=(
["RtrC_input_rep1"]=GSM6065735
["RtrC_input_rep2"]=GSM6065736
["RtrC_input_rep3"]=GSM6065737
["RtrC_ChIP_rep1"]=GSM6065738
["RtrC_ChIP_rep2"]=GSM6065739
["RtrC_ChIP_rep3"]=GSM6065740
)

# Download fastq files
mkdir log
echo "declare -A SRX=(" > ./srx.sh
mkdir -p html
cd html
for sample in "${!geo[@]}"
do
  # wget "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=${geo[${sample}]}&targ=self&form=text&view=full" -O ${geo[${sample}]}.txt
  srx=$(grep SRX ${geo[${sample}]}.txt | head -1 | perl -ne '($t) = $_ =~ /(SRX[0-9]+)/; print $t')
  echo "[\"${sample}\"]=\"${srx}\"" >> ../srx.sh
done
cd ..
echo ")" >> ./srx.sh

## Récupération des SRA
source ./srx.sh
echo "declare -A SRA=(" > ./sra.sh
cd html
for sample in "${!SRX[@]}"
do
  # wget "https://www.ncbi.nlm.nih.gov/sra?term=${SRX[${sample}]}" -O ${SRX[${sample}]}
  sra=$(grep SRR ${SRX[${sample}]} | perl -ne '(@t) = $_ =~ /trace.ncbi.nlm.nih.gov\/Traces\?run=(SRR[0-9]+)/g; print join " ", @t;')
  i=0
  for s in ${sra[@]}
  do
    echo "[\"${sample}${i}\"]=\"${s}\"" >> ../sra.sh
    i=$(( i + 1 ))
  done
done
cd ..
echo ")" >> ./sra.sh

module load sra-tools/3.4.1

source sra.sh
for sample in "${!SRA[@]}"
do
sbatch -o log/dump_${SRA[${sample}]}.nohup \
  -J EB3I.${SRA[${sample}]}.download  <<-EOF
#!/bin/bash
fastq-dump --split-files --gzip -A ${SRA[${sample}]}
EOF
done

## on renomme tous les fichiers 
for sample in "${!SRA[@]}"
do
  mv ${SRA[${sample}]}_1.fastq.gz ${sample}_1.fastq.gz
  mv ${SRA[${sample}]}_2.fastq.gz ${sample}_2.fastq.gz
done  

cd ..

# Download FASTA and GFF3 from the SAME GenBank accession (CP001340.1, C. crescentus NA1000).
EUTILS='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=CP001340.1&retmode=text'
curl -fLsS --retry 3 "${EUTILS}&rettype=fasta" -o data/Caulobacter_CP001340.fasta || exit 1
sleep 1 # NCBI limits requests to 3/s without an API key
curl -fLsS --retry 3 "${EUTILS}&rettype=gff3" -o data/Caulobacter_CP001340.gff3 || exit 1
grep -q '^>CP001340.1' data/Caulobacter_CP001340.fasta || { echo "Invalid FASTA" >&2; exit 1; }
grep -q '^##gff-version' data/Caulobacter_CP001340.gff3 || { echo "Invalid GFF3" >&2; exit 1; }
module load samtools 2>/dev/null || true
# Build data/Caulobacter_CP001340.gtf from the GFF3: one record per 'gene' feature,
# written in 9-column (GFF-style) format with feature type 'exon'.
# This format is deliberate: HOMER (parseGTF.pl, called by annotatePeaks.pl -gtf) only
# keeps records whose feature type is 'exon' (by default), and a plain 'gffread -T'
# conversion of this NCBI GFF3 emits only 212 'exon' records (the tRNA/rRNA/ncRNA genes),
# which silently degrades HOMER peak annotation. One 'exon' record per gene (4,097 records)
# is read completely by HOMER and also displays with gene names in IGV.
awk -F'\t' 'BEGIN{OFS="\t"} NF==9 && $3=="gene"{
  id=""; lt=""; nm=""; bt=""
  if(match($9,/ID=[^;]+/))           id=substr($9,RSTART+3,RLENGTH-3)
  if(match($9,/locus_tag=[^;]+/))    lt=substr($9,RSTART+10,RLENGTH-10)
  if(match($9,/Name=[^;]+/))         nm=substr($9,RSTART+5,RLENGTH-5)
  if(match($9,/gene_biotype=[^;]+/)) bt=substr($9,RSTART+13,RLENGTH-13)
  if(id=="") id=lt
  if(nm=="") nm=(lt!=""?lt:id)
  tid=id; sub(/^(gene|ncRNA|pseudogene|protein_coding|tRNA|rRNA)-/,"rna-",tid)
  gid=id; sub(/^rna-/,"gene-",gid)
  print $1,"Genbank","exon",$4,$5,".",$7,".", \
        "transcript_id \"" tid "\"; gene_id \"" gid "\"; gene_name \"" nm "\"; gene_biotype \"" bt "\";"
}' data/Caulobacter_CP001340.gff3 | sort -k1,1 -k4,4n > data/Caulobacter_CP001340.gtf
echo "GTF records (one per gene, expect 4097):"; wc -l < data/Caulobacter_CP001340.gtf
samtools faidx data/Caulobacter_CP001340.fasta
# Inspect the first sequence name in FASTA and GTF before using HOMER/IGV.
echo 'Reference length:'; awk '{s+=$2} END{print s}' data/Caulobacter_CP001340.fasta.fai

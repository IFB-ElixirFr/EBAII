# IGV Single-cell Hands-on Roscoff 2022


1. [Introduction](#introduction)  
2. [Downloading ChIP-seq reads from NCBI](#download)
3. [Connect to the server and set up your environment](#setup)
4. [Quality control of the reads and statistics](#qc)
5. [Mapping the reads with Bowtie](#mapping)
6. [Estimating the number of duplicated reads](#dup)
7. [ChIP quality controls](#cqc)
8. [Visualizing the data in a genome browser](#visualize)
9. [Peak calling with MACS](#macs)
10. [Motif analysis](#motif)
11. [Peak annotation](#annotation)
12. [Bonus: Peak annotation using R](#peakr)
13. [FAQ](#faq)
14. [References](#ref)

## Introduction <a name="introduction"></a>
### Goal
The aim is to :

  * blah balh
View the data in their genomic context, to check whether the IP worked 

### Summary
Some text


## Connect to the server and set up your environment <a name="setup"></a>
### Connect to JupytherHub
1. Sign in to [Jupyterhub](https://jupyterhub.cluster.france-bioinformatique.fr) and open a Terminal
2. select in the reservation field **form_2022_32**
4. In the launcher, click on "Terminal" in "Other" section. You should be in your home directory by default. Check it:
```bash
pwd
```

### 2 - Set up your working environment
1. Go to your project directory
```bash
cd /shared/projects/<your_project>
```
2. Create a directory that will contain all results of the upcoming analyses.
```bash
mkdir ebaii22_igv
```
3. Go to the newly created directory
```bash
cd ebaii22_igv
```
4. Copy the directory containing data

```bash
cp -r /shared/project/form_22_32/SingleCellRNASeq/Visualization .
```

7. Your directory structure should be like this
 ```
/shared/projects/<your_project>/ebaii22_igv
│
└───Visualization
```

You can check your directory structure:
 ```bash
 tree
```
You should see a BAM file (= mapped reads), a BAI file (for technical reasons this file must be present for IGV) and a BED file (contains a summary of the BAM with alignment position. This file is not mandatory, but it can be useful as less heavy than BAM)
 ```
 └── Visualization
    ├── pbmc1k_rdx.bam
    ├── pbmc1k_rdx.bam.bai
    └── pbmc1k_rdx.bed
```



## Visualizing the signal in a genome browser <a name="visualize"></a>


### 1 - Viewing the raw alignment data in IGV
1. Download the following files from the server onto your computer (laptop)
  * pbmc1k_rdx.bam
  * pbmc1k_rdx.bam.bai
  * pbmc1k_rdx.bed
2. Open IGV on your computer
3. Keep the default genome (GRCh38/hg38)
4. Load BAM file : 
  * File / Load from File...
  * Select the BAM file **pbmc1k_rdx.bam**
5. Load the BED file **pbmc1k_rdx.bed**

(the BAI file does not need to be loaded)

**Browse around in the genome. Do you see peaks?** 
**Where are located the peaks, in refence to the annotated genes ?**   



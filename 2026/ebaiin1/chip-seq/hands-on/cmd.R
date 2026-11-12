## - Associate peaks to closest genes
# 5. Add gene symbol annotation using R with Rstudio

setwd("/shared/projects/2621_eb3in1_2026/atelier_chipseq/elo_test/07-PeakAnnotation")

peaks_annotated <- read.table("RtrC_idr_annotated_peaks.tsv", sep="\t", header=TRUE)



## select the peaks associated with a gene feature
peaks_w_genes <- peaks_annotated$Annotation != "Intergenic"
peaks_annotated_genes <- peaks_annotated[peaks_w_genes,]

## extract the column containing the gene name
selected_genes <- unique(peaks_annotated_genes$Gene.Name)
length(selected_genes) # how many genes are selected ?

## Save the list in a text file
writeLines(selected_genes,"genes_associated_rtrc.txt")


## - Performing a first evaluation of peak sets using R
# 1. Go to Rstudio and execute the R code below (show results in the report)

library(RColorBrewer)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm9.knownGene)
library(org.Mm.eg.db)
# define the annotation of the mouse genome
txdb = TxDb.Mmusculus.UCSC.mm9.knownGene
# define colors
col = brewer.pal(9,'Set1')

# read the peaks for each dataset
peaks.forebrain = readPeakFile('GSM348064_p300_peaks.txt.gz')
peaks.midbrain = readPeakFile('GSM348065_p300_peaks.txt.gz')
peaks.limb = readPeakFile('GSM348066_p300_peaks.txt.gz')

# create a list containing all the peak sets
all.peaks = list(forebrain=peaks.forebrain,
                 midbrain=peaks.midbrain,
                 limb=peaks.limb)

# check the number of peaks for the forebrain dataset
length(peaks.forebrain)

# compute the number of peaks for all datasets using the list object
sapply(all.peaks,length)

# display this as a barplot
barplot(sapply(all.peaks,length),col=col)

# statistics on the peak length for forebrain
summary(width(peaks.forebrain))

# size distribution of the peaks
peaks.width = lapply(all.peaks,width)
lapply(peaks.width,summary)

# boxplot of the sizes
boxplot(peaks.width,col=col)

# genome wide distribution
covplot(peaks.forebrain, weightCol="Maximum.Peak.Height")

# define gene promoters
promoter = getPromoters(TxDb=txdb, upstream=5000, downstream=5000)

# compute the density of peaks within the promoter regions
tagMatrix = getTagMatrix(peaks.limb, windows=promoter)

# plot the density
tagHeatmap(tagMatrix, palette = "RdYlBu")

peakAnno.forebrain = annotatePeak(peaks.forebrain, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
peakAnno.midbrain = annotatePeak(peaks.midbrain, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
peakAnno.limb = annotatePeak(peaks.limb, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")

# distribution of genomic compartments for forebrain peaks
plotAnnoPie(peakAnno.forebrain)

# for all the peaks
plotAnnoBar(list(forebrain=peakAnno.forebrain, midbrain=peakAnno.midbrain,limb=peakAnno.limb))


### - functional annotation

# load the library
library(clusterProfiler)

# define the list of all mouse genes as a universe for the enrichment analysis
universe = mappedkeys(org.Mm.egACCNUM)

## extract the gene IDs of the forebrain target genes
genes.forebrain = peakAnno.forebrain@anno$geneId
ego.forebrain = enrichGO(gene          = genes.forebrain,
                         universe      = universe,
                         OrgDb         = org.Mm.eg.db,
                         ont           = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.01,
                         qvalueCutoff  = 0.05,
                         readable      = TRUE)

# display the results as barplots        
barplot(ego.forebrain,showCategory=10)









## - Associate peaks to closest genes
# 5. Add gene symbol annotation using R with Rstudio

setwd("/shared/projects/2538_eb3i_n1_2025/atelier_chipseq/EBAII2025_chipseq/07-PeakAnnotation")

d <- read.table("FNR_anaerobic_idr_annotated_peaks.tsv", sep="\t", header=TRUE)

gene.symbol <- read.table("../data/Escherichia_coli_K_12_MG1655.annotation.tsv.gz", header=FALSE)

d.annot <- merge(d[,c(1,2,3,4,5,6,8,10,11)], gene.symbol, by.x="Nearest.PromoterID", by.y="V1")

colnames(d.annot)[2] <- "PeakID"  # name the 2d column of the new file "PeakID"
colnames(d.annot)[dim(d.annot)[2]] <- "Gene.Symbol"
write.table(d.annot, "FNR_anaerobic_idr_final_peaks_annotation.tsv", col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)

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









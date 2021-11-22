#!/usr/bin/env R

library("clusterProfiler")  # Make enrichment analysis
library("limma")            # A lots of math-related operations
library("DOSE")             # Disease Ontology
library("enrichplot")       # Awesome graphs
library("pathview")         # Nice pathway plot
library("org.At.tair.db")   # A. Thaliana annotation


##################################
# PART 1 : Introduction
# From DESeq2, we have
# a large table with many columns
##################################

deseq_genes <- read.table(
  file="tables/KOvsWT.complete.txt",
  sep="\t",
  header=TRUE
)
print(names(deseq_genes))

print(head(deseq_genes$Id))


#####################################################
# PART 2: Gene Identifiers
# For a computer, gene:AT1G01010 ≠ AT1G01010
# For a human, AT1G01010 is horrible to remember
# We need to clean identifiers,
# and add human-readable names
#####################################################


# Fix identifiers for the computer
head(deseq_genes$Id)
deseq_genes$Id <- sub("gene:", "", deseq_genes$Id)
head(deseq_genes$Id)

# Translate for for humans
annotation <- bitr(
  geneID   = deseq_genes$Id,         # Our gene list
  fromType = "TAIR",                 # We have TAIR ID
  toType   = c("ENTREZID", "SYMBOL"),# Other ID list
  OrgDb    = org.At.tair.db          # Our annotation
)
print(head(annotation))

# Add the translation to the result table
deseq_genes <- merge(
  x = deseq_genes,  y = annotation,
  by.x = "Id",      by.y = "TAIR"
)
print(head(deseq_genes))


#######################################################
# PART 3: Gene sets
# We need to filter differentially expressed genes
# in order to perform ORA.
#######################################################

# Filter DEG
de_genes <- deseq_genes[deseq_genes[, "padj"] <= 0.001, ]
de_genes <- de_genes[!is.na(de_genes[, "log2FoldChange"]), ]
dim(deseq_genes)
dim(de_genes)

# Perform ORA
ego <- enrichGO(
  gene    = de_genes$ENTREZID,  # Ranked gene list
  universe= deseq_genes$ENTREZID,#All genes
  OrgDb   = org.At.tair.db,     # Annotation
  keyType = "ENTREZID",         # The genes ID
  ont     = "CC",               # Cellular Components
  pvalueCutoff  = 1,            # Significance Threshold
  pAdjustMethod = "BH",         # Adjustment method
  readable = TRUE               # For human beings
)

# Plot ORA results
barplot(ego, showCategory=15)
dotplot(object = ego, showCategory=15)

# Search for enriched terms that relates to roots
res_ego <- ego@result
print(head(res_ego, 3))
roots <- res_ego[with(res_ego, grepl("root", Description)), ]
print(head(roots))


# Correct the previous ORA to include biological processes
# instead of cellular components
ego <- enrichGO(
  gene    = de_genes$ENTREZID,  # Ranked gene list
  universe= deseq_genes$ENTREZID,#All genes
  OrgDb   = org.At.tair.db,     # Annotation
  keyType = "ENTREZID",         # The genes ID
  ont     = "BP",               # Biological Process
  pvalueCutoff  = 1,            # Significance Threshold
  pAdjustMethod = "BH",         # Adjustment method
  readable = TRUE               # For human beings
)

# Plot corrected results
barplot(ego, showCategory=15)
dotplot(object = ego, showCategory=15)

# Search for roots
res_ego <- ego@result
roots <- res_ego[with(res_ego, grepl("root", Description)), ]
print(head(roots))


##################################################################
# PART 4: GSEA
# We need to build a named vector which contains sorted numbers
# Then we can run a Gene set enrichment analysis
###################################################################

# Explore results to guess the right column to extract
print(colnames(deseq_genes))

# We choose to use the 'stat' column
geneList <- as.numeric(de_genes$stat)
names(geneList) <- de_genes$ENTREZID
geneList <- sort(geneList, decreasing=TRUE)

# Run GSEA on Gene Onthology
gsea <- gseGO(
  geneList = geneList,       # Ranked gene list
  ont      = "BP",           # Biological Process
  OrgDb    = org.At.tair.db, # Annotation
  keyType  = "ENTREZID",     # Identifiers
  pAdjustMethod = "BH",      # Pvalue Adjustment
  pvalueCutoff = 1           # Significance Threshold
)

# Explore results
columns_of_interest <- c(
  "Description",
  "enrichmentScore",
  "p.adjust"
)
head(
  x = gsea[, columns_of_interest], # Pathway ID
  8                                # lines to display
)


# Search plant organ morphogenesis
gsea_line <- match(
  "plant organ morphogenesis",
  gsea$Description
)
gseaplot2(
  x         = gsea,               # Our analysis
  geneSetID = gsea$ID[gsea_line],  # Pathway ID
  title     = "plant organ morphogenesis" # Its name
)

# Look at the roots pathways
gsea_line <- match(
  "root morphogenesis",
  gsea$Description
)
gseaplot2(
  x         = gsea,               # Our analysis
  geneSetID = gsea$ID[gsea_line],  # Pathway ID
  title     = "root morphogenesis" # Its name
)

# Plot multiple pathways on the same graph
gseaplot2(
  x = gsea,
  geneSetID = 1:3,
  title = "Most enriched terms"
)


##############################################
# PART 5: Compare sets
# We want to compare pathways to each others
##############################################

# Use a Heatmap
heatplot(
  x =  ego,                   # Our enrichment
  showCategory = 15,          # Nb of terms to display
  foldChange = geneList[1:10] # Our fold changes
)


# Use an upsetplot
upsetplot(x = ego) # From our enrichment analysis

# Use an enrichment map
emapplot(ego) # From our enrichment analysis


# Relate enriched terms with each others:
goplot(ego) # From our enrichment analysis

# Plot complete pathways
#names(geneList) <- de_genes$TAIR  # Use TAIR id

pv.out <- pathview(
  gene.data = geneList,     # Our gene list
  pathway.id = "ath00630",  # Our pathway
  species = "ath",          # Our organism
  # The color limits
  limit = list(gene=max(abs(geneList))),
  gene.idtype = "ENTREZID" # The genes identifiers
)

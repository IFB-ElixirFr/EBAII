#!/bin/R

library("clusterProfiler");
library("limma");
library("DOSE");
library("enrichplot");
library("pathview");
library("org.At.tair.db");

de_res <- read.table(
  file   = "tables/KOvsWT.complete.txt",   # Path to results
  header = TRUE,           # There are column names
  sep    = "\t",           # This is a tabulation
  stringsAsFactors = FALSE # Colnames are not factors
);

# Keep selected columns
de_res <- de_res[ ,c("Id", "padj", "log2FoldChange", "stat")];

# List all values that are to be kept
threshold <- de_res[, "padj"] <= 0.001;

# Keep only non-null and significant values
de_res <- de_res[which(threshold), ];

# Get line order
ordered_lines_number <- order(
  de_res$log2FoldChange,     # Select the column
  decreasing = TRUE            # Decreasing sort
);

# Gather sorted lines
de_res <- de_res[ordered_lines_number, ];

# Replace the names in the ID column
de_res[, "Id"] <- sub("gene:", "", de_res[, "Id"]);


# Translate gene identifiers as genes ID
annotation <- bitr(
  geneID   = de_res$Id,            # Our gene list
  fromType = "TAIR",               # We have TAIR ID
  toType   = c("ENTREZID", "SYMBOL"),# Other ID list
  OrgDb    = org.At.tair.db        # Our annotation
);

# Check ENTREZID
entrez <- data.frame(table(annotation$ENTREZID));
head(entrez[entrez$Freq > 1, ]);

# Check genes symbol
symbol <- data.frame(table(annotation$SYMBOL));
head(symbol[symbol$Freq > 1, ]);

# run enrichment analysis
ego <- enrichGO(
  gene    = annotation$ENTREZID,# Ranked gene list
  OrgDb   = org.At.tair.db,     # Annotation
  keyType = "ENTREZID",         # The genes ID
  ont     = "CC",               # Cellular Components
  pvalueCutoff = 0.001,         # Significance Threshold
  pAdjustMethod = "BH"          # Adjustment method
);

# plot results
barplot(ego, showCategory=20);
dotplot(object = ego, showCategory=20);

# Rename the columns
colnames(de_res) <- c("TAIR", "padj", "log2FoldChange");

# Merge the frames
geneFrame <- merge(de_res, annotation, by="TAIR");
geneFrame <- unique(geneFrame); # Drop duplicates

# Build a numeric vector
geneList <- as.numeric(geneFrame$log2FoldChange);
# Get the genes identifiers
names(geneList) <- geneFrame$ENTREZID;
# Sort this list
geneList <- sort(geneList, decreasing = TRUE);

# Run GSEA
gsea <- gseGO(
  geneList = geneList,       # Ranked gene list
  ont      = "BP",           # Biological Process
  OrgDb    = org.At.tair.db, # Annotation
  keyType  = "ENTREZID",     # Identifiers
  pAdjustMethod = "BH",      # Pvalue Adjustment
  pvalueCutoff = 1           # Significance Threshold
);

# Plot GSEA
# We need the number of the line
# Containing our pathway of interest
gsea_line <- match(
  "plant organ morphogenesis",
  gsea$Description
);
gseaplot2(
  x         = gsea,               # Our analysis
  geneSetID = gsea$ID[gsea_line],  # Pathway ID
  title     = "plant organ morphogenesis" # Its name
);

# Plot multiple GSEA results
gseaplot2(
  x = gsea,
  geneSetID = 1:3,
  title = "Most enriched terms"
);

# Heatmaps
heatplot(
  x =  ego,             # Our enrichment
  showCategory = 15,    # Nb of terms to display
  foldChange = geneList # Our fold changes
);

# Compare sets of pathways and genes
upsetplot(x = ego); # From our enrichment analysis
emapplot(ego); # From our enrichment analysis
goplot(ego); # From our enrichment analysis

# Run Kegg erichment and plot-it
names(geneList) <- geneFrame$TAIR;  # Use TAIR id

pv.out <- pathview(
  gene.data = geneList,     # Our gene list
  pathway.id = "ath00630",  # Our pathway
  species = "ath",          # Our organism
  # The color limits
  limit = list(gene=max(abs(geneList))),
  gene.idtype = "TAIR"      # The genes identifiers
);

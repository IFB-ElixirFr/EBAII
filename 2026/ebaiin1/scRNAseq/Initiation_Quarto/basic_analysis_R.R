## Here I will check for Ptcra expression in my data.
## I expect half cells to express it.


data_dir <- "/shared/home/bjob/gits/EBAII/2026/ebaiin1/scRNAseq/Initiation_Quarto"


## Load data
sobj <- readRDS(paste0(data_dir, "/", "sobj.rds"))

## Check
print(sobj)

## Load the Seurat package
# library(Seurat)     # <= I don't do that anymore, now I borrow functions from the packages !!

## Visualize our favorite gene
Seurat::FeaturePlot(sobj, assay = "RNA", slot = "scale.data", reduction = "tsne", features = "Ptcra", cols = c("white", "red")) + Seurat::DarkTheme() + ggplot2::theme(aspect.ratio = 1)

## => Ptcra is more expressed on the left and lower parts of the tSNE space, but how many cells are concerned ?

## Which cells are positive for Ptcra ?

sobj$is_Ptcra_pos = Seurat::FetchData(sobj, "Ptcra")[["Ptcra"]] > 0
Seurat::DimPlot(object = sobj, reduction = "tsne", group.by = "is_Ptcra_pos") +
Seurat::DarkTheme() + ggplot2::theme(aspect.ratio = 1)

## How many cells are positive ?
table(sobj$is_Ptcra_pos)

## Conclusion : not strictly half of cells are Ptcra+, but close !


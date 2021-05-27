library(mixOmics)

# Liver toxicity

data(liver.toxicity)
help(liver.toxicity)

# Recode the factor Dose
liver.toxicity$treatment$Dose.Group <- as.factor(c(rep("low",32),
                                                   rep("high",32)))

# Modification of rownames
rownames(liver.toxicity$clinic) <- paste0("ID", rownames(liver.toxicity$clinic))

## Question 1: Based on clinical data, do we naturally observe clusters of samples
## which correspond to the different dose or exposure treatments?
## PCA

data.lt.clinic <- liver.toxicity$clinic
res.pca.lt.clinic <- pca(data.lt.clinic, scale=TRUE, ncomp=10)

plotIndiv(res.pca.lt.clinic, ind.names = FALSE,
          title = "plotIndiv PCA: clinical data",
          group=liver.toxicity$treatment$Time.Group,
          pch = as.numeric(factor(liver.toxicity$treatment$Dose.Group)),
          pch.levels =liver.toxicity$treatment$Dose.Group,
          legend = TRUE)

plotVar(res.pca.lt.clinic,
        title = "plotVar PCA: clinical")

## Question 2: based on transcriptomics data, do we naturally observe clusters of
## samples which correspond to the different dose or exposure treatments?
## PCA

data.lt.gene <- liver.toxicity$gene

res.pca.lt.gene <- pca(data.lt.gene, scale=TRUE, ncomp=10)

plotIndiv(res.pca.lt.gene, ind.names = FALSE,
          title = "plotIndiv PCA: transcriptomics data",
          group=liver.toxicity$treatment$Time.Group,
          pch = as.numeric(factor(liver.toxicity$treatment$Dose.Group)),
          pch.levels =liver.toxicity$treatment$Dose.Group,
          legend = TRUE)

plotVar(res.pca.lt.gene, var.names = FALSE)


## Question 3: based on transcriptomics data, do we naturally observe clusters of
## samples which correspond to the different doses or exposure treatments when
## we select some genes highly involved in the variability of the data?
## Sparse PCA

res.spca.lt.gene <- spca(data.lt.gene, scale=TRUE,
                         ncomp=3, keepX=c(10,10,10))

plotIndiv(res.spca.lt.gene, ind.names = FALSE,
          group=liver.toxicity$treatment$Time.Group,
          pch = as.numeric(factor(liver.toxicity$treatment$Dose.Group)),
          pch.levels =liver.toxicity$treatment$Dose.Group,
          legend = TRUE)

plotVar(res.spca.lt.gene, var.names = FALSE, pch=16, 
        title = "plotVar SPCA: transcriptomics")


## Question 4: Based on transcriptomics data, can we identify a molecular signature
## that characterizes the different treatment times?
## Sparse PLS-DA

res.splsda.lt.gene.time <- splsda(data.lt.gene,
                                  liver.toxicity$treatment$Time.Group,
                                  ncomp = 3, keepX = c(10,10,10))

plotIndiv(res.splsda.lt.gene.time, ind.names = FALSE, legend = TRUE)
plotVar(res.splsda.lt.gene.time)

plotLoadings(res.splsda.lt.gene.time, comp=1, contrib = "max")
plotLoadings(res.splsda.lt.gene.time, comp=2, contrib = "max")

## Question 5: Do we observe a better discrimination with the clinical data?
## PLS-DA

res.splsda.lt.clinical.time <- splsda(data.lt.clinic,
                                      liver.toxicity$treatment$Time.Group,
                                      ncomp = 3)

plotIndiv(res.splsda.lt.clinical.time, ind.names = FALSE, legend = TRUE)
plotVar(res.splsda.lt.clinical.time)

plotLoadings(res.splsda.lt.clinical.time, comp=1, contrib = "max")
plotLoadings(res.splsda.lt.clinical.time, comp=2, contrib = "max")

## Question 6: Can we unravel relationships between transcriptomics data and
## clinical data?
## PLS

res.pls.lt.gene.clinic <- pls(data.lt.gene,data.lt.clinic)

plotIndiv(res.pls.lt.gene.clinic,rep.space="XY-variate",
          title = "plotIndiv PLS",
          ind.names=FALSE,
          group=liver.toxicity$treatment$Time.Group,
          pch = as.numeric(factor(liver.toxicity$treatment$Dose.Group)),
          pch.levels =liver.toxicity$treatment$Dose.Group,
          legend = TRUE)

plotVar(res.pls.lt.gene.clinic, var.names = c(FALSE, TRUE))

## Question 7: Can we unravel relationships between transcriptomics data and
## clinical data? What are the genes that characterize these relationships?
## Sparse PLS

res.spls.lt.gene.clinic <- spls(data.lt.gene,data.lt.clinic,
                                ncomp=3, keepX = c(10,10,10))

plotIndiv(res.spls.lt.gene.clinic,rep.space="XY-variate",
          title = "plotIndiv Sparse PLS",
          ind.names=FALSE,
          group=liver.toxicity$treatment$Time.Group,
          pch = as.numeric(factor(liver.toxicity$treatment$Dose.Group)),
          pch.levels =liver.toxicity$treatment$Dose.Group,
          legend = TRUE)

plotVar(res.spls.lt.gene.clinic, var.names = c(FALSE, TRUE))

plotLoadings(res.spls.lt.gene.clinic, comp=1, size.title = 1)

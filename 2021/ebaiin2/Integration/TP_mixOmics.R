library(mixOmics)

data(srbct)

help(srbct)

X <- srbct$gene

#### PCA

res.pca.srbct <- pca(X)
?pca

plotVar(res.pca.srbct)
?plotvar
plotVar(res.pca.srbct, var.names = FALSE)

plotIndiv(res.pca.srbct)
?plotIndiv
plotIndiv(res.pca.srbct, group = srbct$class)
plotIndiv(res.pca.srbct, group = srbct$class, ind.names = FALSE)
plotIndiv(res.pca.srbct, group = srbct$class,
          ind.names = FALSE, legend=TRUE)

res.pca.srbct.3pc <- pca(X, ncomp=3)

#### Sparse PCA

res.spca.srbct <- spca(X, ncomp=3, keepX=c(10,10,10))
?spca
plotIndiv(res.spca.srbct, group = srbct$class)
plotVar(res.spca.srbct)

#### PLS-DA

Y <- srbct$class

res.plsda.srbct <- plsda(X, Y)
?plsda

plotIndiv(res.plsda.srbct)
plotIndiv(res.plsda.srbct, ind.names=FALSE, legend=TRUE)
plotIndiv(res.plsda.srbct, ind.names=FALSE, ellipse = TRUE, legend=TRUE)

#### Sparse PLS-DA

res.splsda.srbct <- splsda(X, Y, ncomp = 3, keepX = c(10,10,10))
?splsda

plotIndiv(res.splsda.srbct, ind.names=FALSE, legend=TRUE)
plotIndiv(res.splsda.srbct, ind.names=FALSE, legend=TRUE,
          comp = c(1,3))
plotVar(res.splsda.srbct)

plotLoadings(res.splsda.srbct, comp=1, contrib = 'max')
?plotLoadings
plotLoadings(res.splsda.srbct, comp=2, contrib = 'max')
plotLoadings(res.splsda.srbct, comp=3, contrib = 'max')


#### PLS

data(nutrimouse)

X <- nutrimouse$gene
Y <- nutrimouse$lipid

# PCA first!
## Genes
res.pca.nutrimouse.gene <- pca(X)
plotIndiv(res.pca.nutrimouse.gene)
plotIndiv(res.pca.nutrimouse.gene, ind.names = FALSE, legend=TRUE,
          group=nutrimouse$genotype, pch=nutrimouse$diet)
plotIndiv(res.pca.nutrimouse.gene, ind.names = FALSE, legend=TRUE,
          pch=nutrimouse$genotype, group=nutrimouse$diet)
plotIndiv(res.pca.nutrimouse.gene, ind.names = nutrimouse$diet, legend=TRUE,
          group=nutrimouse$genotype)
plotVar(res.pca.nutrimouse.gene)

## Lipids
res.pca.nutrimouse.lipid <- pca(Y)
plotIndiv(res.pca.nutrimouse.lipid, ind.names = nutrimouse$diet, legend=TRUE,
          group=nutrimouse$genotype)
plotVar(res.pca.nutrimouse.lipid)

# PLS
res.pls.nutrimouse <- pls(X,Y)
plotIndiv(res.pls.nutrimouse)
plotIndiv(res.pls.nutrimouse, rep.space = "XY-variate")
plotIndiv(res.pls.nutrimouse, rep.space = "XY-variate",
          ind.names = nutrimouse$diet, legend=TRUE,
          group=nutrimouse$genotype)
plotVar(res.pls.nutrimouse)

## S-PLS

res.spls.nutrimouse <- spls(X,Y, ncomp = 3, keepX=c(5,5,5))
plotIndiv(res.spls.nutrimouse, rep.space = "XY-variate",
          ind.names = nutrimouse$diet, legend=TRUE,
          group=nutrimouse$genotype)
plotVar(res.spls.nutrimouse)

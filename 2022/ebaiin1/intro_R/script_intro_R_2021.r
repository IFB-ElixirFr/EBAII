########################################################################
#                     Introduction à R... en 1h45                      #
#                                                                      #
# Hugo Varet (Institut Pasteur)                                        #
# Olivier Kirsh (Univ. Paris Diderot)                                  #
# Thomas Denecker (I2BC)                                               #
# Jacques van Helden (Univ. Marseille / IFB)                           #
# Ecole de Bioinformatique AVIESAN-IFB                                 #
# Roscoff 2021                                                         #
########################################################################

#############################################
#        R vu comme une calculatrice        #
#############################################
2 + 3
4 * 5
6 / 4

#############################################
#          notion de variable/objet         #
#############################################
a <- 2
print(a)
b <- 3
resultat <- a + b
print(resultat)

a <- resultat
print(a)

#############################################
#    cas pratique : données d'expression    #
#############################################
# chargement des données
exprs <- read.table(file = "expression.txt", header = TRUE, sep = "\t")

# affichage de l'objet "exprs"
print(exprs)

# accéder à l'aide d'une fonction
?read.table
help(read.table)

# affichage des premières lignes de l'objet
head(exprs)

# un peu plus de lignes
head(exprs, n = 20)

# dimensions du tableau
dim(exprs)

# résumé rapide de l'ensemble des données
summary(exprs)

# valeurs stockées dans la colonne WT1
exprs$WT1
exprs[, "WT1"]

# histogramme de ces valeurs
hist(exprs$WT1)

# histogramme du logarithme de ces valeurs
hist(log(exprs$WT1))

# Exprssions KO1 vs WT1
plot(x = log(exprs$WT1), y = log(exprs$KO1))

# Personnalisation des paramètres graphiques
plot(x = log(exprs$WT1), y = log(exprs$KO1),
     main = "Expression KO1 vs WT1",
     xlab = "WT1", ylab = "KO1", pch = 16, col = "red")

# sélection des lignes 4 et 11 du tableau des expressions
exprs[c(4, 11), ]

# sélection des gènes ENSG00000059728 et ENSG00000140807
exprs[which(exprs$id %in% c("ENSG00000059728", "ENSG00000140807")), ]

# expressions moyennes des WT et des KO
exprs$meanWT <- rowMeans(exprs[ , c("WT1", "WT2")])
exprs$meanKO <- rowMeans(exprs[ , c("KO1", "KO2")])

# fold-change KO vs WT
exprs$FC <- exprs$meanKO / exprs$meanWT

# moyenne de tous les échantillons
exprs$mean <- rowMeans(exprs[ , c("WT1", "WT2", "KO1", "KO2")])

# MA-plot: log2FC vs intensité
plot(x = log(exprs$mean), y = log2(exprs$FC),
     pch = 16, xlab = "Intensity", ylab = "log2FC")

# charger les annotations des gènes
annot <- read.table(file = "annotation.csv", header = TRUE, sep = ";")
print(annot)

# combien de gènes par chromosome
table(annot$chr)

# selectionner les données d'expression seulement pour le chromosome 8
# 1ere étape: merger les deux tableaux exprs et annot
exprs.annot <- merge(exprs, annot, by = "id")
head(exprs.annot)
# 2eme étape: sous-ensemble des lignes pour lesquelles chr vaut 8
exprs8 <- exprs.annot[which(exprs.annot$chr == 8),]
print(exprs8)

# exporter exprs8 dans un fichier
write.table(x = exprs8, file = "exprs8.txt", sep = "\t",
            row.names = TRUE, col.names = NA)

#############################################
#                 packages R                #
#############################################
# CRAN: 16353 packages (sept 2020) avec croissance exponentielle
# Bioconductor: 1903 "software" packages (aout 2020)

# installation depuis le CRAN: exemple du package packHV
install.packages("packHV")

# installation depuis Bioconductor: exemple du package DESeq2
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
BiocManager::install("DESeq2")

# chargement d'un package
library(packHV)


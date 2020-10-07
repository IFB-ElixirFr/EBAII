# TP croisement de données

## Préparation de l'environnement de travail

Dans votre dossier de projet, création d'un répertoire spécifique pour ce TP de croisement de données:
```{bash}
mkdir tp_croisement
cd tp_croisement
```
Chargement du programme bedtools pour le cluster:
```{bash}
module avail bedtools
module load bedtools/2.27.1
```
Astuce pour racourcir l'écriture des commandes en utilisant la commande bash "ln":
création d'un dossier (lienData) qui servira de lien vers le repertoire /shared/projects/ebaii2020/atelier_croisement/data/
```{bash}
ln -s /shared/projects/ebaii2020/atelier_croisement/data/ lienData
```

## Question 1

### pour les gènes différentiellement exprimés

Extraction des positions des gènes différentiellement exprimés:
```{bash}
srun grep -f lienData/differentially_expressed_genes.txt /shared/bank/homo_sapiens/GRCh38/gff3/Homo_sapiens.GRCh38.94.gff3 | srun awk '($3 == "gene")' > differentially_expressed_genes.gff
```

Comparaison (par les positions sur le génome) des gènes différentiellement exprimés avec les régions des sites de fixation de H3K4me3:
```{bash}
srun bedtools intersect -a differentially_expressed_genes.gff -b lienData/h3k4me3_k562.bed -u > differentially_expressed_genes_h3k4me3_k562.gff
```

### pour tous les gènes (background)

Extraction des positions de tous les gènes:
```{bash}
srun awk '($3 == "gene")' /shared/bank/homo_sapiens/GRCh38/gff3/Homo_sapiens.GRCh38.94.gff3 > all_genes.gff
```

Sélection des gènes non différentiellement exprimés:
```{bash}
srun bedtools intersect -v -a all_genes.gff -b differentially_expressed_genes.gff > not_differentially_expressed_genes.gff
```
Comparaison (par les positions sur le génome) des gènes non différentiellement exprimés avec les régions des sites de fixation de H3K4me3:
```{bash}
srun bedtools intersect -a not_differentially_expressed_genes.gff -b lienData/h3k4me3_k562.bed -u > not_differentially_expressed_genes_h3k4me3_k562.gff
```

Résultat des comparaisons (en terme de nombre de lignes des fichiers résultats):
```{bash}
wc -l differentially_expressed_genes.gff differentially_expressed_genes_h3k4me3_k562.gff not_differentially_expressed_genes.gff not_differentially_expressed_genes_h3k4me3_k562.gff
```

## Question 2


Extraction des régions promotrices des gènes différentiellement exprimés:
```{bash}
srun bedtools flank -i differentially_expressed_genes.gff -l 2000 -r 0 -g lienData/chrs.len -s > differentially_expressed_genes_prom2k.gff
```
Récupération des variants présents dans les régions promotrices:
```{bash}
srun bedtools intersect -a lienData/common_all_20180418_div.vcf -b differentially_expressed_genes_prom2k.gff > differentially_expressed_genes_snps_in_prom2k.vcf
```
Dénombrement des variants :
```{bash}
wc -l differentially_expressed_genes_snps_in_prom2k.vcf
```
Dénombrement des variants pour chaque promoteur des gènes différentiellement exprimés (option -c):
```{bash}
srun bedtools intersect -a differentially_expressed_genes_prom2k.gff3 -b lienData/common_all_20180418_div.vcf -c > differentially_expressed_genes_prom2k_n_snps.gff
```
Observation du résultat (aller à la dernière colonne):
```{bash}
less -S differentially_expressed_genes_prom2k_n_snps.gff
```
Contrôle (calcul de la somme de cette dernière colonne):
```{bash}
srun awk -F '\t' 'BEGIN{s=0}{s+=$NF}END{print s}' differentially_expressed_genes_prom2k_n_snps.gff
```
Le compte est bon !

## Question 3

Croisement des régions des sites de fixation de H3K4me3 avec la liste des variants:
```{bash}
srun bedtools intersect -a lienData/h3k4me3_k562.bed -b lienData/common_all_20180418_div.vcf -u > h3k4me3_k562_snps.bed
```
Réorganisation des résultats (sort) par chromosome (alphanumérique) puis par position (numérique):
```{bash}
srun sort -k1,1 -k2,2n h3k4me3_k562_snps.bed > h3k4me3_k562_snps_sorted.bed
```
Recherche des gènes différentiellement exprimés et à proximité des sites de fixation de H3K4me3 contenant un variant:
```{bash}
srun bedtools closest -a h3k4me3_k562_snps_sorted.bed -b differentially_expressed_genes.gff > differentially_expressed_genes_closest_h3k4me3_k562_snps.bed
```
Idem en retrignant à la région amont (upstrem) des gènes, en ajoutant des filtres (awk) pour contraindre la région à 2Kb (et en focalisant les résultats sur les gènes):
```{bash}
srun bash -c "bedtools closest -a h3k4me3_k562_snps_sorted.bed -b differentially_expressed_genes.gff -id | awk -F "\t" '($14 <= 0 && $14 >= -2000)' | cut -f 5- | cut -f -9 | sort -u > differentially_expressed_genes_closest_h3k4me3_k562_snps_selected.bed"
```

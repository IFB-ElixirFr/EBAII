## Exemple de solution du TP de croisement de données, en version "chemins absolus"
## WARNING : ce script est à utiliser dans un environnement contrôlé (Jupyter) sur une machine de faible puissance (1CPU 1GB RAM suffit)
## WARNING : pour l'utiliser directement sur le cluster ifb, il faut recourir à "srun" pour toutes les commandes non-système (et toute commande derrière un pipe "|" )

# PREPARATION
## Aller dans un espace de travail
### Dans ce script, je prends comme exemple mon 'home', mais vous pouvez utiliser votre espace projet
cd /shared/home/${USER}
## Créer le répertoire de travail de ce TP, et y aller
mkdir tp_croisement
cd tp_croisement

## Chargement de bedtools
### Obtenir la liste des modules dont le nom commence par "bedtools"
module avail bedtools
### Charger l'outil bedtools dans sa version la plus récente disponible dans les différents modules
module load bedtools
### Contrôle du chargement effectif
module list

## Environnement prêt pour travailler, on peut y aller !


# QUESTION 1 (gènes différentiellement exprimés et marques d'histone)

## Extraction du GFF des gènes différentiellement exprimés de ma liste TXT (isssue d'un test différentiel)
### On va commencer par observer le contenu des différents fichiers à utiliser : less -S
less -S /shared/projects/form_2021_26/data/atelier_croisement/differentially_expressed_genes.txt
less -S /shared/bank/homo_sapiens/GRCh38/gff3/Homo_sapiens.GRCh38.94.gff3
### On dénombre les lignes
wc -l /shared/bank/homo_sapiens/GRCh38/gff3/Homo_sapiens.GRCh38.94.gff3 # 2814993
wc -l /shared/projects/form_2021_26/data/atelier_croisement/differentially_expressed_genes.txt # 162
### On sélectionne du GFF total les lignes correspondant à nos gènes différentiellement exprimés
grep -f /shared/projects/form_2021_26/data/atelier_croisement/differentially_expressed_genes.txt \
/shared/bank/homo_sapiens/GRCh38/gff3/Homo_sapiens.GRCh38.94.gff3 > differentially_expressed_from_genelist.gff
### On dénombre combien de lignes ont ainsi été sélectionnées
wc -l /shared/projects/form_2021_26/data/atelier_croisement/differentially_expressed_genes.txt \
differentially_expressed_from_genelist.gff # 1106 ?!
## PROBLEME : on a beaucoup de lignes qui ne sont pas des gènes (transcrits, linc, etc) dont il faut les filtrer pur ne conserver que les entrées marquées 'gene' dans la colonne 3
less -S differentially_expressed_from_genelist.gff
awk '($3 == "gene")' differentially_expressed_from_genelist.gff > differentially_expressed_genes.gff
wc -l /shared/projects/form_2021_26/data/atelier_croisement/differentially_expressed_genes.txt \
differentially_expressed_genes.gff # 162 : OK !::
less -S differentially_expressed_genes.gff

## Intersection A=gff_exp_global B=bed_histone_global : donne un gff restreint aux régions communes
## '-u' (unique) ne rapporte qu'une occurrence de A s'il y a plusieurs matches avec B
less -S /shared/projects/form_2021_26/data/atelier_croisement/h3k4me3_k562.bed
wc -l /shared/projects/form_2021_26/data/atelier_croisement/h3k4me3_k562.bed # 52422
bedtools intersect -a differentially_expressed_genes.gff \
-b /shared/projects/form_2021_26/data/atelier_croisement/h3k4me3_k562.bed \
-u > differentially_expressed_genes_h3k4me3_k562.gff
## Sélection des entrées du gff croisé pour lesquelles la 3e colonne a la valeur "gene", de la même façon que tout à l'heure
awk '($3 == "gene")' /shared/bank/homo_sapiens/GRCh38/gff3/Homo_sapiens.GRCh38.94.gff3 > ALL_genes.gff
wc -l ALL_genes.gff # 21492 genes au total
## intersection "inversée" (exclusion) de ce GFF gènes totaux avec le GFF des gènes différentiellement exprimés, pour avoir les gènes non différentiellement exprimés
## -v : Ne donne que les entrées de A qui N'ONT PAS d'overlap avec B (similaire à 'grep -v' dans sa syntaxe)
bedtools intersect -v -a ALL_genes.gff \
-b differentially_expressed_genes.gff > NOT_differentially_expressed_genes.gff
wc -l NOT_differentially_expressed_genes.gff # 21180 genes non différentiels (tous les gènes moins 162)
## intersection du gff des non-différentiels avec le bed de méthylation
bedtools intersect -a NOT_differentially_expressed_genes.gff \
-b /shared/projects/form_2021_26/data/atelier_croisement/h3k4me3_k562.bed \
-u > NOT_differentially_expressed_genes_h3k4me3_k562.gff
## dénombrement pour chaque gff
wc -l differentially_expressed_genes.gff \
differentially_expressed_genes_h3k4me3_k562.gff \
NOT_differentially_expressed_genes.gff \
NOT_differentially_expressed_genes_h3k4me3_k562.gff
## Petite ligne en Perl pour faire le calcul des ratios et afficher cela joliment
perl -e 'print "\n\t\tINACT\t/ ALL\nDIFF:\t\t77\t/ 162\t= ".(77/162)."\nNON-DIFF:\t9362\t/ 21180\t= ".(9362/21180)."\n\n"'
## REPONSE 1 : NA !



# QUESTION 2

## Créer de nouvelles régions de 2Kb en amont de chaque gene différentiel pour symbolyser le promoteur avec 'bedtools flank'
## -l 2000 : 2Kb à gauche de mes régions du GFF
## -r 0 : pas indispensable, mais par sécurité je précise que je ne veux pas de région flanquante à droite du gène
## -g : fichier génome donnant la position de fin de chaque chromosome (pour capper les nouveaux intervales à générer)
## -s : pour définir l'orientation de '-l' non pas sur le sens positif du génome, mais sur le strand du gène (lu dans le GFF)
bedtools flank -i differentially_expressed_genes.gff -l 2000 -r 0 \
-g /shared/projects/form_2021_26/data/atelier_croisement/chrs.len \
-s > differentially_expressed_genes_prom2k.gff
## Je compare les coordonnées que j'obtiens pour les 2 premiers gènes (1er en brin+, 2e en brin-) avec ses coordonnées originelles
head -n 2 differentially_expressed_genes.gff
head -n 2 differentially_expressed_genes_prom2k.gff

## Bonus : Intersection "par la gauche uniquement" avec le paramètre "-loj" : donne le nombre de promoteurs ayant des SNP dedans
wc -l /shared/projects/form_2021_26/data/atelier_croisement/common_all_20180418_div.vcf
bedtools intersect -a differentially_expressed_genes_prom2k.gff \
-b /shared/projects/form_2021_26/data/atelier_croisement/common_all_20180418_div.vcf \
-loj > differentially_expressed_genes_prom2k_snps.gff # 450

## Pour répondre réellement à la question : croisement du vcf des SNP avec les régions promoteur des gènes différentiels
bedtools intersect -a /shared/projects/form_2021_26/data/atelier_croisement/common_all_20180418_div.vcf \
-b differentially_expressed_genes_prom2k.gff  > snps_in_prom2k.vcf # 405
wc -l snps_in_prom2k.vcf

## REPONSE 2 : 405 variants SNP

## Question complémentaire : dénombrer les SNPs matchés par promoteur :
## Option '-c' pour donner le dénombrement
bedtools intersect -a differentially_expressed_genes_prom2k.gff \
-b /shared/projects/form_2021_26/data/atelier_croisement/common_all_20180418_div.vcf \
-c > differentially_expressed_genes_prom2k_n_snps.gff
## Vérification (somme des valeurs de comptage du gff)
awk -F '\t' 'BEGIN{s=0}{s+=$NF}END{print s}' differentially_expressed_genes_prom2k_n_snps.gff # 405

## REPONSE 2prime : valeur par ligne dont la somme fait bien 405


# QUESTION 3

## Croisement pics de méthylations avec les variants du VCF
bedtools intersect -a /shared/projects/form_2021_26/data/atelier_croisement/h3k4me3_k562.bed \
-b /shared/projects/form_2021_26/data/atelier_croisement/common_all_20180418_div.vcf \
-u > h3k4me3_k562_snps.bed

## Tri sur les colonnes 1 (chr) puis 2 (start, en numérique) car 'bedtools closest' requière obligatoirement un fichier trié
sort -k1,1 -k2,2n h3k4me3_k562_snps.bed > h3k4me3_k562_snps_sorted.bed

## Recherche des régions proximales des gènes (mode par défaut : avec ou sans overlap, des deux côtés)
bedtools closest -a h3k4me3_k562_snps_sorted.bed \
-b differentially_expressed_genes.gff3 > differentially_expressed_genes_closest_h3k4me3_k562_snps.bed

##
## closest :  '-D' donne la distance signée entre B et A dans une colonne supplémentaire,
##              la valeur "b" spécifie que le signe est donné en respect de l'orientation de l'entrée de B
## awk :      '-F' séparateur de champ (tabulation)
##            expression : tous les champs dont la colonne 14 est comprise dans les 2 Kb en amont (d'où
##              valeur négative) de B
## cut :      on récupère toutes les colonnes de la 5e jusqu'à la fin de ligne, puis de la 1e à la 9e (à tester en 5-13)
## sort :     on trie en mode unique ('-u')
bedtools closest -a h3k4me3_k562_snps_sorted.bed \
-b differentially_expressed_genes.gff3 -D b \
| awk -F "\t" '($14 <= 0 && $14 >= -2000)' \
| cut -f 5- \
| cut -f -9 \
| sort -u > differentially_expressed_genes_closest_h3k4me3_k562_snps_selected.bed

## REPONSE 3 : contenu du bed final

## J'ai fini, je peux décharger mon usage des modules
module unload bedtools ## ou "module purge" pour faire une ardoise blanche
echo YOUPI!

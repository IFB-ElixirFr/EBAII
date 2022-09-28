
# 1ère école de bioinformatique AVIESAN - IFB - INSERM (EBAII 2022)

## Assemblage et annotation de novo de génomes

**Du 25-09-2022 au 30-09-2022 à la Station biologique de Roscoff**


## Web site

- Web pages: <https://ifb-elixirfr.github.io/EBAII/2022/Assemblage>

From there, you can also find a link to download the whole repository with `git`.

****

## Accès à Galaxy

Tous les TPs se feront sur [UseGalaxy.fr](https://usegalaxy.fr). Chaque participant doit se connecter à son compte personnel (à créer avant la formation).

Pour bénéficier des ressources de calcul réservées pour cette formation, chaque participant doit :

- se connecter sur [UseGalaxy.fr](https://usegalaxy.fr)
- cliquer sur [ce lien](https://usegalaxy.fr/join-training/ebaii_aa/)

## Supports de cours

### Cours généraux

| Cours                             | Intervenants         | Supports                                                                                                                                         |
|-----------------------------------|----------------------|--------------------------------------------------------------------------------------------------------------------------------------------------|
| Introduction                      | A. Cormier; E. Corre | [Cours](https://training.galaxyproject.org/training-material/topics/assembly/tutorials/get-started-genome-assembly/slides.html)                  |
| Outils de séquençage              | J. Kreplak           | [Cours](https://docs.google.com/presentation/d/1rtTCyVF4dz0Trmny5e8r1brzfNab8ZUN/edit?usp=sharing&ouid=109813995176155673766&rtpof=true&sd=true) |
| QC des lectures                   | O. Rué               | [Cours](https://drive.google.com/file/d/1Mv33oQ-_h-ZCxemvlcqYqpQEHd97tJzt/view?usp=sharing)                                                      |
| Validation d'assemblage (Theorie) | E. Corre; A. Cormier | [Cours](https://training.galaxyproject.org/training-material/topics/assembly/tutorials/assembly-quality-control/slides.html)                     |
| Validation d'assemblage (TP)      | A. Cormier; E. Corre | [TP](https://training.galaxyproject.org/training-material/topics/assembly/tutorials/assembly-quality-control/tutorial.html)                      |


### Assemblage

| Cours         | Intervenants                     | Supports                                                                                                      |
|---------------|----------------------------------|---------------------------------------------------------------------------------------------------------------|
| Assemblage    | D. Naquin; C. Klopp              | [Cours](Genome_assembly.pdf) <br> [TP](Genome_assembly_tp.pdf) <br> [Conclusion TP](conclusion_TP.pdf)        |
| Polishing     | J.-M. Aury                       | [Cours](https://docs.google.com/presentation/d/1RAScBkXvWkRCuD2WAbgNLJZ8zJNXz9skkHJ-MGp4VBk/edit?usp=sharing) |
| Scaffolding   | J. Kreplak; C. Klopp; A. Cormier | [Cours](https://drive.google.com/file/d/1SRBBqRPUUTePJ7K1wsqbmaFGqAuvVIt6/view?usp=sharing)                   |

### Annotation

#### Annotation procaryotes

| Cours                         | Intervenants | Supports                                                                                                                                                                             |
|-------------------------------|--------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Introduction                  | G. Gautreau  | [1. Introduction to prokaryotic genomes annotation](https://github.com/IFB-ElixirFr/EBAII/files/9662414/1.Introduction.to.prokaryotic.genomes.annotation.pdf)                        |
| Dataset                       | H. Chiapello | [2.1 Dataset Construction](https://github.com/IFB-ElixirFr/EBAII/files/9659668/2.1-Dataset-Construction-EBAIIA.A.pdf)                                                                |
| Analyse                       | H. Chiapello | [2.2 Dataset Analysis](https://github.com/IFB-ElixirFr/EBAII/files/9659684/2.2-Dataset-Analysis-EBAIIA.A.pdf)                                                                        |
| Déreplication                 | H. Chiapello | [2.3 Dataset Dereplication](https://github.com/IFB-ElixirFr/EBAII/files/9659690/2.3-Dataset-Dereplication-EBAIIA.A.pdf)                                                              |
| CDS prediction                | G. Gautreau  | [3. Prediction of coding genes](https://github.com/IFB-ElixirFr/EBAII/files/9663528/3.Prediction.of.coding.genes.pdf)                                                              |
| Non coding gene prediction    | G. Gautreau  | [4. Prediction of non coding genes](https://github.com/IFB-ElixirFr/EBAII/files/9662493/4.Prediction.of.non.coding.genes.pdf)                                                    |
| Functional annotation (intro) | G. Gautreau  | [5. Introduction to functional annotation](https://github.com/IFB-ElixirFr/EBAII/files/9662561/5.Introduction.to.functional.annotation.pdf)                                      |
| Prokka + Vizualisation        | G. Gautreau  | [6. and 7. Let's play with a prokaryotic annotation tool_ Prokka](https://github.com/IFB-ElixirFr/EBAII/files/9662563/6._7.Let.s.play.with.a.prokaryotic.annotation.tool_.Prokka.pdf)|
| InterScan + EggNOG            | J. Joets     | [8. FunctAnnot.jjoets](https://github.com/IFB-ElixirFr/EBAII/files/9662567/FunctAnnot.jjoets.pdf)                                                                                   |

#### Annotation eucaryotes

| Cours                    | Intervenants       | Supports |
|--------------------------|--------------------|----------|
| Repeat masking           | Jonathan Kreplak   | [Cours](https://drive.google.com/file/d/1rcF9d7ZG4gPMrMYjt3vhwGzCGhdF2BVy/view?usp=sharing) <br> [TP](https://training.galaxyproject.org/training-material/topics/genome-annotation/tutorials/repeatmasker/tutorial.html) |
| Annotation de gènes      | Anthony Bretaudeau | [Cours](https://training.galaxyproject.org/training-material/topics/genome-annotation/slides/introduction.html)<br>[TP](https://training.galaxyproject.org/topics/genome-annotation/tutorials/funannotate/tutorial.html) |
| Visualisation            | Anthony Bretaudeau | [TP](https://training.galaxyproject.org/topics/genome-annotation/tutorials/funannotate/tutorial.html#visualisation-with-a-genome-browser) |
| ARN long non-codants     | Stéphanie Robin    | [Cours](FEELnc_Sept_2022.pdf)<br>[TP](https://training.galaxyproject.org/topics/genome-annotation/tutorials/lncrna/tutorial.html) |
| BUSCO                    | Stéphanie Robin    | [Cours](BUSCO_Sept_2022.pdf)<br>[TP]() |
| Annotation fonctionnelle | Johann Joëts       | [Cours](https://drive.google.com/file/d/1vP2NMW0c0aOWSRHKwU3FNHzHRovg9GMp/view?usp=sharing)<br>TP |
| Orthologie               | Johann Joëts       |  Cours<br>TP |


### Analyses tierces

| Cours      | Intervenants | Supports |
|------------|--------------|----------|
| Procaryote | À venir      |  À venir |
| Eucaryote  | À venir      |  À venir |

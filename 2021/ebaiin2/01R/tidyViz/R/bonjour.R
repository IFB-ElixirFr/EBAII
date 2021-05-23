#' La fonction qui dit bonjour
#' 
#' "Bonjour !" : Cette fonction est un exemple de fonction graphique avec deux paramètres simples : la couleur du texte et la couleur du fond.
#' 
#' PS : cette fonction ne servira jamais pour le cours, elle sert juste à dire bonjour !
#' 
#' @param couleur la couleur du texte
#' @param fond la couleur du fond
#'
#' @return Un graphe simple
#' @export
#'
#' @examples
#' bonjour()
bonjour <- function(couleur = "pink", fond = "limegreen") {
  base::plot(0:1, 0:1, type = "n", axes = F, xlab ="", ylab = "",)
  parusr <- graphics::par("usr")
  graphics::rect(parusr[1], parusr[3], parusr[2], parusr[4], col = fond, border = NA)
  graphics::text(x = 0.5, y = 0.5, labels = "Bonjour !", col = couleur, cex = 8, font = 2)
}

#' Le jeu de données fruits
#' 
#' Composition nutritionnelle de 51 composés alimentaires issus du groupe "fruits" défini dans les données Ciqual de 2020.
#'
#' La table fruits est de classe "tibble" et contient les colonnes suivantes : 
#'\itemize{
#'  \item nom : Nom du fruit
#'  \item groupe : Groupe de fruit
#'  \item Energie : Energie, (kJ/100 g)
#'  \item Eau : Eau (g/100 g)
#'  \item Proteines : Protéines, N x 6.25 (g/100 g)
#'  \item Glucides : Glucides (g/100 g)
#'  \item Lipides : Lipides (g/100 g)
#'  \item Sucres : Sucres (g/100 g)
#'  \item Fructose : Fructose (g/100 g)
#'  \item Fibres : Fibres alimentaires (g/100 g)
#'  \item Calcium : Calcium (mg/100 g)
#'  \item Magnesium : Magnésium (mg/100 g)
#'  \item Phosphore : Phosphore (mg/100 g)
#'  \item Potassium : Potassium (mg/100 g)
#'  \item Zinc : Zinc (mg/100 g)
#'  \item BetaCarotene : Beta-Carotène (µg/100 g)
#'  \item VitamineE : Vitamine E (mg/100 g)
#'  \item VitamineC : Vitamine C (mg/100 g)
#'}
#'  
#' @name fruits
#' @docType data
#' @references \url{https://ciqual.anses.fr/}
#' @keywords data
NULL


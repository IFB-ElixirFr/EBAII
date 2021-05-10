#' La fonction qui dit bonjour
#' 
#' "Bonjour !" : Cette fonction est un exemple de fonction graphique avec deux paramètres simples : la couleur du texte et la couleur du fond.
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
#' Composition nutritionnelle de 57 composés alimentaires issus du groupe "fruits" défini dans les données Ciqual de 2020.
#'
#' @name fruits
#' @docType data
#' @references \url{https://ciqual.anses.fr/}
#' @keywords data
NULL


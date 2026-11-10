assess_covar <- function(mat = NULL, covar.df = NULL, markers = NULL, ndim = 10L, center = TRUE, scale = TRUE, plot = TRUE, title = NULL, return.weights = FALSE, return.plot = FALSE) {
  ## Checks
  if(is.null(mat)) stop('A matrix is required !')
  if(is.null(covar.df)) stop('A data.frame with covariates is required !')
  if(!is.data.frame(covar.df)) stop('covar.df should be a data.frame !')
  if (ncol(mat) != nrow(covar.df)) stop("ncol(mat) != nrow(covar.df) !")
  if(is.null(ndim)) stop('A number of dimensions to use is required !')
  if(!is.integer(ndim)) stop('ndim must be an integer > 0 !')
  if(ndim > ncol(mat)) stop('ndim must be <= ndim(mat) !')
  if(all(c(return.weights, return.plot))) stop("Can't return both weight matrix and plot !")
  
  ## Convert to pure R matrix
  mat <- as.matrix(mat)
  
  ## Check markers
  if(!is.null(markers)) markers <- markers[markers %in% rownames(mat)]
  
  ## Split factors and continuous data
  conti.names <- factor.names <- c()
  for (x in colnames(covar.df)) {
    if(is.numeric(covar.df[1,x])) conti.names <- c(conti.names, x) else factor.names <- c(factor.names, x)
  }
  col.names <- c(factor.names, conti.names, markers)
  col.types <- c(rep('factor', length(factor.names)), rep('continuous', length(conti.names)), rep('marker', length(markers)))
  
  ## Clean factors
  if (length(factor.names) > 0) {
    for (ccn in factor.names) {
      covar.df[[ccn]] <- as.factor(covar.df[[ccn]])
      covar.df[[ccn]] <- droplevels(covar.df[[ccn]])
    }
  }
  
  ## Remove factors with single level
  # covar.df[,vapply(seq_len(ncol(covar.df)), nlevels)]
  fac_sing <- vapply(covar.df, function(x) { if (!is.factor(x)) FALSE else if (nlevels(x) > 1) FALSE else TRUE }, TRUE)
  if (any(fac_sing)) {
    covar.df <- covar.df[,!fac_sing]
    factor.names <- factor.names[!factor.names %in% colnames(covar.df)[fac_sing]]
  }
  
  ## Convert conti to Zscores
  if (length(conti.names) > 0) {
    for (x in conti.names) {
      covar.df[[x]] <- (covar.df[[x]] - mean(covar.df[[x]], na.rm = TRUE)) / sd(covar.df[[x]], na.rm = TRUE)
    }
  }
  
  ## Convert markers to Zscores
  if (length(markers) > 0) {
    covar.df <- cbind(covar.df, t(mat[markers,]))
    for (x in markers) {
      covar.df[[x]] <- (covar.df[[x]] - mean(covar.df[[x]], na.rm = TRUE)) / sd(covar.df[[x]], na.rm = TRUE)
    }
  }
  
  ## Center / scale the matrix ?
  if (any(c(center, scale))) mat <- base::scale(x = mat, center = center, scale = scale)
  
  ## Dimension reduction
  norm.red <- irlba::irlba(A = t(mat), nv = min(ndim+1, ncol(mat)))$u
  dimnames(norm.red) <- list(colnames(mat), paste0("PC_", 1:(ndim+1)))
  
  rm(mat)
  
  ## Setting output matrix
  bc.mat <- matrix(NA, ncol = length(col.names), nrow = ndim, dimnames = list(colnames(norm.red)[1:ndim], col.names))
  ## Filling matrix
  for (cn in seq_along(col.names)) {
    # message(col.names[cn])
    if (col.names[cn] %in% c(conti.names, markers)) {
      cv2cor <- covar.df[[col.names[cn]]]
      nona <- !is.na(cv2cor)
      bc.mat[, cn] <-  abs(cor(x = cv2cor[nona], y = norm.red[nona,seq_len(ndim)], method = 'spearman'))
    } else if (col.names[cn] %in% factor.names) {
      b2kw <- covar.df[[col.names[cn]]]
      nona <- !is.na(b2kw)
      for (si in seq_len(ndim)) {
        k_test <- try(k_res <- kruskal.test(x = norm.red[nona,si], g = as.factor(b2kw[nona])), silent = TRUE)
        if (!is(k_test, class2 = 'try-error')) {
          bc.mat[si,cn] <- k_res$statistic / nrow(norm.red)
        } else bc.mat[si,cn] <- 0
      }
    }
  }
  ## Heatmap
  color.palette = c("white", "orangered3")
  myRamp.col <- circlize::colorRamp2(c(0, 1), color.palette)
  BC.hm <- ComplexHeatmap::Heatmap(matrix = bc.mat,
                                   name = 'Weight',
                                   col = myRamp.col,
                                   na_col = 'grey75',
                                   cluster_rows = FALSE,
                                   cluster_columns = FALSE,
                                   rect_gp = grid::gpar(col = "darkgrey", lwd=0.5),
                                   column_title = paste0(c('Covariates weight', title), collapse = ' ; '),
                                   row_title = 'Dimensions',
                                   column_split = col.types,
                                   top_annotation = ComplexHeatmap::HeatmapAnnotation(
                                     Type = col.types, 
                                     col = list(Type = setNames(
                                       object = c('lightblue','pink', 'green3'), 
                                       nm = c('factor', 'continuous', 'marker')))))
  if (plot) print(BC.hm)
  
  ## Return the weight matrix ?
  if(return.plot) return(BC.hm)
  
  ## Return the weight matrix ?
  if(return.weights) return(bc.mat)
}

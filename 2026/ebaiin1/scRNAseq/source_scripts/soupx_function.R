## Allows to run SoupX with or without the unfiltered matrix (filtered matrix is required). By default, does not return the "unsouped" matrix, just the soup rate estimate (rho). Unsouped matrix is returned when 'return_object' is set to TRUE.
SoupX_auto <- function(scmat_filt = NULL,
                       scmat_raw = NULL,
                       soupQuantile = 0.9,
                       contaminationRange = c(.01, .8),
                       soupRange = c(0,100),
                       return_object = FALSE,
                       doPlot = FALSE) {
  
  ## Checks
  if(is.null(scmat_filt)) {
    stop('A filtered count matrix is required !')
  }
  
  if(is.null(scmat_raw)) {
    message('No unfiltered raw counts matrix provided. Estimation will be based on filtered matrix only.') # <9>
  }
  
  ## If no raw matrix
  if (is.null(scmat_raw)) {
    spChanRaw <- SoupX::SoupChannel(
      tod = scmat_filt, 
      toc = scmat_filt, 
      calcSoupProfile = FALSE)
    sc_rowsum <- sparseMatrixStats::rowSums2(scmat_filt)
    spProf <- data.frame(
      row.names = rownames(scmat_filt), 
      est = sc_rowsum/sum(scmat_filt), 
      counts = sc_rowsum)
    spChan <- SoupX::setSoupProfile(spChanRaw, spProf)
  } else { # <11>
    spChan <- SoupX::SoupChannel( # <11>
      tod = scmat_raw,  # <11>
      toc = scmat_filt,  # <11>
      calcSoupProfile = FALSE) # <11>
    if (min(spChan$nDropUMIs) > max(soupRange)) stop( # <11>
      'Minimum found counts per barcode is : ',  # <11>
      min(spChan$nDropUMIs),  # <11>
      ', which is smaller than the upper bound of soupRange !', # <11>
      ' Please increase soupRange max !') # <11>
    spChan <- SoupX::estimateSoup(sc = spChan,  # <11>
                                  soupRange = soupRange) # <11>
  } # <11>
  ## Display Top 20 contributing genes
  if (!return_object) {
    cat('\nSoup-contributing features (Top 20) :\n')
    print(knitr::kable(head(
      spChan$soupProfile[order(spChan$soupProfile$est, 
                               decreasing = TRUE), ], 
      n = 20)))
  }
  ## Quick clustering needed
  spClust <- scran::quickCluster(scmat_filt, method = "igraph")
  ## Adding clusters to the SoupChannel object
  spChan <- SoupX::setClusters(sc = spChan, clusters = spClust)
  ## Estimating soup
  sX <- SoupX::autoEstCont(sc = spChan, 
                           doPlot = doPlot, 
                           tfidfMin = 1, 
                           soupQuantile = soupQuantile, 
                           maxMarkers = 100, 
                           contaminationRange = contaminationRange, 
                           rhoMaxFDR = .2, 
                           priorRho = .05, 
                           priorRhoStdDev = .1, 
                           forceAccept = FALSE)
  
  ## Removing soup (adjusting counts)
  if(return_object) {
    cat('Counts BEFORE SoupX : ', sum(scmat_filt), '\n')
    scmat_soupx <- SoupX::adjustCounts(
      sX, 
      method = 'subtraction', 
      roundToInt = TRUE, 
      tol = .001, 
      pCut = .01)
    cat('Counts AFTER SoupX : ', sum(scmat_soupx), '\n')
    rm(scmat_filt)
    return(scmat_soupx)
  } else return(sX$fit$rhoEst)
}

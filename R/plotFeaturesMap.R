plotFeaturesMap <- function(ptn,
<<<<<<< HEAD
                           regOnly = TRUE,
                           comparisons = NULL,
                           featSel,
                           remBinary = TRUE,
                           featCol = NULL,
                           scaled = FALSE,
                           remExtreme = NULL,
                           pdfName = NULL
                           ){
  

  check_ptn(ptn)
  
  if(!is.null(ptn_features(ptn))){
    featuresIn <- ptn_features(ptn)
  } else {
    stop("features analysis is null, please run featureIntegration analysis first")
  }
  if(!check_featSel(featSel, features = featuresIn)){
    stop("Please provide character vector of features. The recommedation is to provide these that were selected in the final featureIntegration model.
         These can be obtained by running ptn_selectedFeatures function and providing ptn object, analysis_type (ie lm or rf) and comparison (it would be a consecutive number from the comparisons provided running featureIntegration)")
  }
  
  if(!is.null(featCol)){
    if(!check_featSel(featCol, features = featuresIn)){
      stop("Please provide character vector of features to be plotted. All of them must be already calculated, to check what these are run: colnames(ptn_features(ptn)) ")
    }
  }
  featuresIn <- featuresIn[,colnames(featuresIn) %in% featSel]
  
  if(!check_logical(regOnly)){
=======
                            regOnly = TRUE,
                            comparisons = NULL,
                            featSel,
                            scaled = TRUE,
                            remExtreme = NULL,
                            pdfName = NULL) {
  check_ptn(ptn)

  if (!is.null(ptn_features(ptn))) {
    featuresTmp <- ptn_features(ptn)
  } else {
    stop("features analysis is null, please run featureIntegration analysis first")
  }
  if (!check_featSel(featSel, features = featuresTmp)) {
    stop("Please provide character vector of features to be plotted. The recommedation is to provide these that were selected in the final featureIntegration model.
         These can be obtained by running ptn_selectedFeatures function and providing ptn object, analysis_type (ie lm or rf) and comparison , it would be a consecutivenumber from the comparisons provided running featureIntegration")
  }
  featuresTmp <- featuresTmp[, colnames(featuresTmp) %in% featSel]
  if (!check_logical(regOnly)) {
>>>>>>> 15c7b1170844d1e76dfc836b2cb76526172b0a48
    stop("'regOnly' can only be TRUE or FALSE")
  }
  if (!check_logical(scaled)) {
    stop("'scaled' can only be TRUE or FALSE")
  }
<<<<<<< HEAD
  if(!check_logical(remBinary)){
    stop("'remBinary' can only be TRUE or FALSE")
  }
  if(isTRUE(regOnly)){
    if(!check_comparisons(comparisons)){
=======
  if (isTRUE(regOnly)) {
    if (!check_comparisons(comparisons)) {
>>>>>>> 15c7b1170844d1e76dfc836b2cb76526172b0a48
      stop("'comparisons' must be a list of numeric vector for paired comparisons example: list(c(0,2),c(0,1)). 0 is always a background.")
    }
    if (length(which(unique(unlist(comparisons)) == 0)) > 0 && is.null(ptn_background(ptn))) {
      stop(" 0 is always a background, but no background provided")
    }
    if (length(comparisons) != 1) {
      stop("Although potentially posssible, please run each comparison separetely to fit well to your feature selection")
    }
  }
  if (!is.null(remExtreme)) {
    if (!check_number(remExtreme)) {
      stop("remExtreme should be a number indicating percentile of extreme values to be removed for colouring purposes ")
    } else {
      if (remExtreme <= 0 | remExtreme >= 1) {
        stop("remExtreme can be only 1 number between (0,1). It represents percentile in decimal of extreme values from both sides of distribution to be removed only for colouring scale")
      }
    }
  }
  #
  if (length(which(unique(unlist(comparisons)) == 0)) > 0 && is.null(ptn_background(ptn))) {
    stop(" 0 is always a background, but no background provided")
  }
  #
  if (isTRUE(regOnly)) {
    resOut <- resQuant(qvec = ptn_effect(ptn), ptn = ptn)
    # for (i in 1:length(comparisons)) {
    if (names(resOut)[1] == "background") {
      compTmp <- comparisons[[1]] + 1
    } else {
      compTmp <- comparisons[[1]]
    }
    listSel <- c(names(resOut[[compTmp[1]]]), names(resOut[[compTmp[2]]]))
<<<<<<< HEAD
    featuresIn<- featuresIn[row.names(featuresIn) %in% listSel,]
  } 
  if(isTRUE(remBinary)){
    binaryCols <- sapply(featuresIn, is_binary)
    featuresCluster <- featuresIn[, !binaryCols, drop = FALSE]
=======
    features <- featuresTmp[row.names(featuresTmp) %in% listSel, ]
  }
  if (isTRUE(scaled)) {
    featuresOut <- scale(na.omit(features), center = T)
  } else {
    featuresOut <- na.omit(features)
>>>>>>> 15c7b1170844d1e76dfc836b2cb76526172b0a48
  }
  if(isTRUE(scaled)){
    featuresCluster <- scale(na.omit(featuresCluster), center=T)
  } else {
    featuresCluster <- na.omit(featuresCluster)
  }

  fmapRes <- umap::umap(featuresCluster, n_components = 2)
  fmapRes <- fmapResOut <- as.data.frame(fmapRes$layout)
<<<<<<< HEAD
  colnames(fmapRes) <- colnames(fmapResOut) <- c("UMAP1", "UMAP2")
  fmapRes$Gene <- rownames(featuresCluster)

  #Plot selected features
  
  effTmp <- ptn_effect(ptn)
  eff <- effTmp[match(row.names(fmapRes),names(effTmp))]
  effect_fmap <- plot_fmap(fmapRes, colVec = eff, remExtreme = remExtreme, name='Effect')
  
  if(is.null(featCol)){
    featCol <- featSel
  }
  
  featuresIn <- featuresIn[match(row.names(fmapRes),row.names(featuresIn)),]
  for(feat in featCol){

    featTmp <- featuresIn[,feat]
    feature_fmap <- plot_fmap(fmapRes, colVec = featTmp, remExtreme = remExtreme, name=feat)
    
    pdf(ifelse(is.null(pdfName),paste(feat, '_featureUMAP.pdf', sep=''), paste(pdfName,feat, 'featureUMAP.pdf',sep='_')),width= 16,height=8, useDingbats = F)
    par(mar=c(5,5,5,5),bty='l',font=2, font.axis=2, font.lab=2, cex.axis=0.9, cex.main=0.7,cex.lab=0.9)
    
=======
  colnames(fmapRes) <- colnames(fmapResOut) <- c("fUMAP1", "fUMAP2")
  fmapRes$Gene <- rownames(featuresOut)

  # Plot every feature together with eff. Prepare effect plot first
  effTmp <- ptn_effect(ptn)
  eff <- effTmp[match(row.names(featuresOut), names(effTmp))]

  effect_fmap <- plot_fmap(fmapRes, colVec = eff, remExtreme = remExtreme, name = "Effect")

  for (feat in featSel) {
    featTmp <- features[, feat]
    # is_binary(featTmp)

    feature_fmap <- plot_fmap(fmapRes, colVec = featTmp, remExtreme = remExtreme, name = feat)

    pdf(ifelse(is.null(pdfName), paste(feat, "_fumap.pdf", sep = ""), paste(pdfName, feat, "fumap.pdf", sep = "_")), width = 16, height = 8, useDingbats = F)
    par(mar = c(5, 5, 5, 5), bty = "l", font = 2, font.axis = 2, font.lab = 2, cex.axis = 0.9, cex.main = 0.7, cex.lab = 0.9)

>>>>>>> 15c7b1170844d1e76dfc836b2cb76526172b0a48
    gridExtra::grid.arrange(
      effect_fmap$legend[[1]],
      effect_fmap$mainPlot,
      feature_fmap$mainPlot,
      feature_fmap$legend[[1]],
      ncol = 4,
      widths = c(1, 4, 4, 1)
    )

    dev.off()
  }
  ptn@analysis@featureIntegration@featureMap <- fmapResOut
  return(ptn)
}

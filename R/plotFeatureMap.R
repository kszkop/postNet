plotFeatureMap <- function(ptn,
                           regOnly = TRUE,
                           comparisons = NULL,
                           featSel,
                           scaled = TRUE,
                           remExtreme = NULL,
                           pdfName
                           ){
  
  check_ptn(ptn)
  
  if(!is.null(ptn_features(ptn))){
    featuresTmp <- ptn_features(ptn)
  } else {
    stop("features analysis is null, please run featureIntegration analysis first")
  }
  if(!check_featSel(featSel, features = featuresTmp)){
    stop("Please provide character vector of features to be plotted. The recommedation is to provide these that were selected in the final featureIntegration model.
         These can be obtained by running ptn_selectedFeatures function and providing ptn object, analysis_type (ie lm or rf) and comparison , it would be a consecutivenumber from the comparisons provided running featureIntegration")
  }
  featuresTmp <- featuresTmp[,colnames(featuresTmp) %in% featSel]
  if(!check_logical(regOnly)){
    stop("'regOnly' can only be TRUE or FALSE")
  }
  if(!check_logical(scaled)){
    stop("'scaled' can only be TRUE or FALSE")
  }
  if(isTRUE(regOnly)){
    if(!check_comparisons(comparisons)){
      stop("'comparisons' must be a list of numeric vector for paired comparisons example: list(c(0,2),c(0,1)). 0 is always a background.")
    }
    if(length(which(unique(unlist(comparisons))==0))>0 && is.null(ptn_bg(ptn))){
      stop(" 0 is always a background, but no background provided")
    }
    if(length(comparisons) != 1){
      stop('Although potentially posssible, please run each comparison separetely to fit well to your feature selection')
    }
  }
  if(!is.null(remExtreme)){
    if(!check_number(remExtreme)){
      stop('remExtreme should be a number indicating percentile of extreme values to be removed for colouring purposes ')
    } else {
      if(remExtreme <= 0 | remExtreme >= 1){
        stop('remExtreme can be only 1 number between (0,1). It represents percentile in decimal of extreme values from both sides of distribution to be removed only for colouring scale')
      }
    }
  }
  #
  if(length(which(unique(unlist(comparisons))==0))>0 && is.null(ptn_bg(ptn))){
    stop(" 0 is always a background, but no background provided")
  }
  resOut <- resQuant(qvec = ptn_effect(ptn), ptn = ptn)
  #
  if(isTRUE(regOnly)){
    #for (i in 1:length(comparisons)) {
    if (names(resOut)[1] == 'background') {
      compTmp <- comparisons[[i]] + 1
    } else {
      compTmp <- comparisons[[i]]
    }
    listSel <- c(names(resOut[[compTmp[1]]]), names(resOut[[compTmp[2]]]))
    features <- featuresTmp[row.names(featuresTmp) %in% listSel,]
  } 
  if(isTRUE(scaled)){
    featuresOut <- scale(na.omit(features), center=F)
  } else {
    featuresOut <- na.omit(features)
  }
  fmapRes <- umap::umap(featuresOut, n_components = 2)
  fmapRes <- fmapResOut <- as.data.frame(fmapRes$layout)
  colnames(fmapRes) <- colnames(fmapResOut) <- c("fMAP1", "fMAP2")
  fmapRes$Gene <- rownames(featuresOut)
  
  #Plot every feature together with eff. Prepare effect plot first
  effTmp <- ptn_effect(ptn)
  eff <- effTmp[match(row.names(featuresOut),names(effTmp))]
  
  effect_fmap <- plot_fmap(fmapRes, colVec = eff, remExtreme = remExtreme, name='Effect')
  
  for(feat in featSel){
    featTmp <- features[,feat]
    feature_fmap <- plot_fmap(fmapRes, colVec = featTmp, remExtreme = remExtreme, name=feat)
    
    pdf(ifelse(is.null(pdfName),paste(feat, '_fmap.pdf', sep=''), paste(pdfName,feat, 'fmap.pdf',sep='_')),width= 16,height=8, useDingbats = F)
    par(mar=c(5,5,5,5),bty='l',font=2, font.axis=2, font.lab=2, cex.axis=0.9, cex.main=0.7,cex.lab=0.9)
    
    gridExtra::grid.arrange(
      effect_fmap$legend,
      effect_fmap$mainPlot,
      feature_fmap$mainPlot,
      feature_fmap$legend,
      ncol = 4,
      widths = c(1, 4, 4, 1))
    
    dev.off()
  }
}